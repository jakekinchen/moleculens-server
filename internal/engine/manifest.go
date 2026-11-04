package engine

import (
	"bytes"
	"compress/gzip"
	"crypto/sha256"
	"encoding/base64"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"

	"github.com/jakekinchen/moleculens-server/internal/contract"
)

const MaxResult = 32 << 20
const MaxArtifacts = 256 << 20

type File struct {
	Name   string `json:"name"`
	Size   int64  `json:"size"`
	SHA256 string `json:"sha256"`
}
type Manifest struct {
	Protocol      int               `json:"protocol"`
	Kind          string            `json:"kind"`
	Key           string            `json:"key"`
	Token         string            `json:"token"`
	Compatibility string            `json:"compatibility"`
	InputSHA256   string            `json:"inputSha256"`
	Provenance    map[string]string `json:"provenance"`
	Files         []File            `json:"files"`
	Payload       string            `json:"payload"`
	ElapsedMS     float64           `json:"elapsedMs"`
}

func validName(s string) bool {
	return s != "" && s != "." && s != ".." && filepath.Base(s) == s && !strings.ContainsAny(s, "/\\")
}
func readBounded(path string, max int64) ([]byte, error) {
	info, err := os.Lstat(path)
	if err != nil {
		return nil, err
	}
	if !info.Mode().IsRegular() || info.Size() > max {
		return nil, errors.New("artifact not a bounded regular file")
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	b, err := io.ReadAll(io.LimitReader(f, max+1))
	if int64(len(b)) > max {
		return nil, errors.New("artifact exceeds limit")
	}
	return b, err
}
func Load(root string, m *Manifest) ([]byte, error) {
	if len(m.Files) == 0 || len(m.Files) > 64 || m.Payload != "response.json" {
		return nil, errors.New("invalid artifact manifest")
	}
	seen := map[string]bool{}
	var total int64
	var payload []byte
	for _, f := range m.Files {
		if !validName(f.Name) || seen[f.Name] || f.Size <= 0 {
			return nil, errors.New("invalid artifact entry")
		}
		seen[f.Name] = true
		if f.Size > MaxArtifacts-total {
			return nil, errors.New("artifact budget exceeded")
		}
		total += f.Size
		path := filepath.Join(root, f.Name)
		if f.Name != m.Payload {
			info, e := os.Lstat(path)
			if e != nil {
				return nil, e
			}
			if !info.Mode().IsRegular() || info.Size() != f.Size {
				return nil, errors.New("artifact size/type mismatch")
			}
			stream, e := os.Open(path)
			if e != nil {
				return nil, e
			}
			hash := sha256.New()
			n, e := io.Copy(hash, io.LimitReader(stream, f.Size+1))
			stream.Close()
			if e != nil {
				return nil, e
			}
			if n != f.Size || hex.EncodeToString(hash.Sum(nil)) != f.SHA256 {
				return nil, errors.New("artifact integrity failure")
			}
			continue
		}
		data, err := readBounded(path, MaxResult)
		if err != nil {
			return nil, err
		}
		h := sha256.Sum256(data)
		if int64(len(data)) != f.Size || hex.EncodeToString(h[:]) != f.SHA256 {
			return nil, errors.New("artifact integrity failure")
		}

		if f.Name == m.Payload {
			if len(data) > MaxResult {
				return nil, errors.New("response too large")
			}
			payload = data
		}
	}
	if payload == nil {
		return nil, errors.New("payload absent")
	}
	return payload, nil
}
func Validate(root string, r *contract.Request, token string) (*Manifest, []byte, error) {
	data, err := readBounded(filepath.Join(root, "manifest.json"), 128<<10)
	if err != nil {
		return nil, nil, err
	}
	var m Manifest
	if err = json.Unmarshal(data, &m); err != nil {
		return nil, nil, err
	}
	if m.Protocol != contract.Protocol || m.Kind != r.Kind || m.Key != r.Key || m.Token != token || m.Compatibility != r.Compatibility || m.InputSHA256 != r.InputSHA256 {
		return nil, nil, errors.New("manifest identity mismatch")
	}
	if m.ElapsedMS < 0 || !finite(m.ElapsedMS) || m.Provenance["build"] == "" || m.Provenance["python"] == "" {
		return nil, nil, errors.New("missing build provenance")
	}
	if r.Kind == "orbital" && m.Provenance["psi4"] == "" {
		return nil, nil, errors.New("Psi4 provenance absent")
	}
	payload, err := Load(root, &m)
	if err != nil {
		return nil, nil, err
	}
	if err = ValidatePayload(r, payload); err != nil {
		return nil, nil, err
	}
	return &m, payload, nil
}
func finite(f float64) bool { return !math.IsInf(f, 0) && !math.IsNaN(f) }
func obj(v any) (map[string]any, error) {
	m, ok := v.(map[string]any)
	if !ok {
		return nil, errors.New("object required")
	}
	return m, nil
}
func count(m map[string]any, key string, max int) (int, error) {
	f, ok := m[key].(float64)
	if !ok || !finite(f) || f < 0 || f > float64(max) || f != math.Trunc(f) {
		return 0, fmt.Errorf("invalid %s", key)
	}
	return int(f), nil
}
func buffer(v any, n int, isIndex bool, vertices int) error {
	s, ok := v.(string)
	if !ok || len(s) > MaxResult {
		return errors.New("missing or oversized buffer")
	}
	compressed, err := base64.StdEncoding.DecodeString(s)
	if err != nil {
		return err
	}
	z, err := gzip.NewReader(bytes.NewReader(compressed))
	if err != nil {
		return err
	}
	defer z.Close()
	data, err := io.ReadAll(io.LimitReader(z, int64(n)*4+1))
	if err != nil {
		return err
	}
	if len(data) != n*4 {
		return errors.New("buffer length mismatch")
	}
	for i := 0; i < len(data); i += 4 {
		u := binary.LittleEndian.Uint32(data[i : i+4])
		if isIndex {
			if uint64(u) >= uint64(vertices) {
				return errors.New("mesh index out of range")
			}
		} else if !finite(float64(math.Float32frombits(u))) {
			return errors.New("nonfinite mesh value")
		}
	}
	return nil
}
func mesh(v any, potential bool) error {
	m, err := obj(v)
	if err != nil {
		return err
	}
	nv, err := count(m, "vertexCount", 2_000_000)
	if err != nil {
		return err
	}
	nt, err := count(m, "triangleCount", 4_000_000)
	if err != nil {
		return err
	}
	if nv == 0 || nt == 0 {
		return errors.New("empty mesh must be explicit null")
	}
	for _, k := range []string{"vertices", "normals"} {
		if err = buffer(m[k], nv*3, false, 0); err != nil {
			return fmt.Errorf("%s: %w", k, err)
		}
	}
	if err = buffer(m["indices"], nt*3, true, nv); err != nil {
		return err
	}
	if potential {
		return buffer(m["potential"], nv, false, 0)
	}
	return nil
}
func vector(v any) error {
	a, ok := v.([]any)
	if !ok || len(a) != 3 {
		return errors.New("three vector components required")
	}
	for _, v := range a {
		f, ok := v.(float64)
		if !ok || !finite(f) {
			return errors.New("invalid vector")
		}
	}
	return nil
}
func ValidatePayload(r *contract.Request, data []byte) error {
	var p map[string]any
	if err := json.Unmarshal(data, &p); err != nil {
		return err
	}
	if len(p) == 0 {
		return errors.New("empty result")
	}
	switch r.Kind {
	case "orbital":
		meta, e := obj(p["meta"])
		if e != nil {
			return e
		}
		if meta["method"] != r.Params["method"] || meta["basis"] != r.Params["basis"] || meta["gridSpacingAngstrom"] != r.Params["gridSpacing"] || meta["charge"] != r.Params["charge"] || meta["multiplicity"] != r.Params["multiplicity"] || meta["densityIsovalue"] != r.Params["densityIsovalue"] {
			return errors.New("orbital metadata mismatch")
		}
		orbs, e := obj(p["orbitals"])
		if e != nil {
			return e
		}
		for _, v := range r.Params["orbitals"].([]any) {
			name := v.(string)
			if name == "density" {
				if _, exists := p["density"]; !exists {
					return errors.New("density field missing")
				}
				if p["density"] != nil {
					if e = mesh(p["density"], false); e != nil {
						return e
					}
				}
				continue
			}
			o, e := obj(orbs[name])
			if e != nil {
				return fmt.Errorf("requested orbital missing: %s", name)
			}
			if o["isovalue"] != r.Params["isovalue"] {
				return errors.New("orbital isovalue mismatch")
			}
			energy, ok := o["energyEv"].(float64)
			if !ok || !finite(energy) {
				return errors.New("orbital energy missing")
			}
			for _, sign := range []string{"positive", "negative"} {
				m, exists := o[sign]
				if !exists {
					return errors.New("orbital lobe missing")
				}
				if m != nil {
					if e = mesh(m, false); e != nil {
						return e
					}
				}
			}
		}
	case "electrostatics":
		s, e := obj(p["surface"])
		if e != nil {
			return e
		}
		if e = mesh(s["mesh"], true); e != nil {
			return e
		}
		units, e := obj(s["units"])
		if e != nil {
			return e
		}
		if units["coords"] != "angstrom" || units["potential"] != "kcal/mol/e" {
			return errors.New("unsupported result units")
		}
		d, e := obj(p["dipole"])
		if e != nil {
			return e
		}
		for _, k := range []string{"originAngstrom", "vectorDebye"} {
			if e = vector(d[k]); e != nil {
				return e
			}
		}
		mag, ok := d["magnitudeDebye"].(float64)
		if !ok || mag < 0 || !finite(mag) || d["convention"] != "physics" {
			return errors.New("invalid dipole")
		}
		v := d["vectorDebye"].([]any)
		sum := 0.
		for _, x := range v {
			sum += x.(float64) * x.(float64)
		}
		if math.Abs(math.Sqrt(sum)-mag) > 1e-5 {
			return errors.New("dipole magnitude mismatch")
		}
		meta, e := obj(p["meta"])
		if e != nil {
			return e
		}
		if meta["method"] != r.Params["method"] || meta["molecularCharge"] != r.Params["charge"] || meta["multiplicity"] != r.Params["multiplicity"] {
			return errors.New("electrostatics metadata mismatch")
		}
	case "conformer":
		sdf, ok := p["sdf3d"].(string)
		if !ok {
			return errors.New("conformer SDF missing")
		}
		if e := contract.ValidateSDF(sdf); e != nil {
			return e
		}
		meta, e := obj(p["meta"])
		if e != nil {
			return e
		}
		if meta["geomVersion"] != r.Params["geomVersion"] || meta["method"] != r.Params["method"] {
			return errors.New("conformer metadata mismatch")
		}
		if _, e = obj(p["timingMs"]); e != nil {
			return e
		}
	case "nomenclature":
		name, ok := p["displayName"].(string)
		if !ok || name == "" || p["source"] != "openclatura" {
			return errors.New("invalid naming result")
		}
		if _, ok = p["verified"].(bool); !ok {
			return errors.New("missing naming verification")
		}
		if _, ok = p["trace"].([]any); !ok {
			return errors.New("missing naming trace")
		}
	default:
		return errors.New("unknown result kind")
	}
	return nil
}

// Seal only after the process group was reaped. Files are immutable thereafter.
func Seal(root string, m *Manifest) error {
	for _, f := range m.Files {
		p := filepath.Join(root, f.Name)
		fd, e := os.Open(p)
		if e != nil {
			return e
		}
		e = fd.Sync()
		fd.Close()
		if e != nil {
			return e
		}
		if e = os.Chmod(p, 0444); e != nil {
			return e
		}
	}
	d, e := os.Open(root)
	if e != nil {
		return e
	}
	defer d.Close()
	return d.Sync()
}
