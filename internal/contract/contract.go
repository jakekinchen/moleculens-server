// Package contract defines the bounded HTTP and chemistry adapter protocol.
package contract

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"math"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

const MaxInput = 2_000_000
const Protocol = 1

var Versions = map[string]string{"orbital": "psi4-2", "electrostatics": "esp-2", "conformer": "rdkit-conformer-2", "nomenclature": "openclatura-1"}

type Request struct {
	Protocol      int            `json:"protocol"`
	Kind          string         `json:"kind"`
	Compatibility string         `json:"compatibility"`
	Input         string         `json:"input"`
	InputSHA256   string         `json:"inputSha256"`
	Params        map[string]any `json:"params"`
	Fallback      string         `json:"fallback"`
	Key           string         `json:"key"`
	WaitMS        int            `json:"-"`
}

func digest(b []byte) string { h := sha256.Sum256(b); return hex.EncodeToString(h[:]) }

// Identity encodes typed values and exact IEEE-754 parameters, with sorted map
// keys. No rounded geometry, caller key or dependency build label overrides it.
func canonical(v any) string {
	switch t := v.(type) {
	case string:
		b, _ := json.Marshal(t)
		return "s" + string(b)
	case float64:
		return fmt.Sprintf("f%016x", math.Float64bits(t))
	case bool:
		if t {
			return "b1"
		}
		return "b0"
	case nil:
		return "n"
	case []any:
		var b strings.Builder
		b.WriteString("[")
		for _, x := range t {
			b.WriteString(canonical(x))
			b.WriteByte(',')
		}
		b.WriteString("]")
		return b.String()
	case map[string]any:
		keys := make([]string, 0, len(t))
		for k := range t {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		var b strings.Builder
		b.WriteString("{")
		for _, k := range keys {
			b.WriteString(canonical(k))
			b.WriteByte(':')
			b.WriteString(canonical(t[k]))
			b.WriteByte(',')
		}
		b.WriteString("}")
		return b.String()
	default:
		panic(fmt.Sprintf("non canonical type %T", v))
	}
}
func (r *Request) Identity() string {
	return digest([]byte(canonical(map[string]any{"protocol": float64(Protocol), "kind": r.Kind, "compatibility": r.Compatibility, "inputSha256": r.InputSHA256, "params": r.Params, "fallback": r.Fallback})))
}

var aliases = map[string]string{"sdf_content": "sdfContent", "grid_spacing": "gridSpacing", "density_isovalue": "densityIsovalue", "inchi_key": "inchiKey", "client_cache_key": "clientCacheKey", "probe_radius": "probeRadius", "softening_epsilon": "softeningEpsilon", "clamp_percentiles": "clampPercentiles", "geom_version": "geomVersion", "wait_ms": "waitMs", "max_attempts": "maxAttempts", "max_opt_iters": "maxOptIters", "add_hs": "addHs", "verify_opsin": "verifyOpsin"}

func alias(m map[string]any) error {
	for k, v := range m {
		if target, ok := aliases[k]; ok {
			if old, exists := m[target]; exists && canonical(old) != canonical(v) {
				return fmt.Errorf("conflicting aliases: %s and %s", k, target)
			}
			m[target] = v
			delete(m, k)
		}
	}
	return nil
}
func allowed(m map[string]any, names string) error {
	set := map[string]bool{}
	for _, k := range strings.Fields(names) {
		set[k] = true
	}
	for k := range m {
		if !set[k] {
			return fmt.Errorf("unsupported field %q", k)
		}
	}
	return nil
}
func number(m map[string]any, key string, def, lo, hi float64, integer bool) (float64, error) {
	v, ok := m[key]
	if !ok {
		return def, nil
	}
	n, ok := v.(float64)
	if !ok || math.IsNaN(n) || math.IsInf(n, 0) || n < lo || n > hi || (integer && math.Trunc(n) != n) {
		return 0, fmt.Errorf("invalid %s", key)
	}
	if n == 0 {
		n = 0
	}
	return n, nil
}
func str(m map[string]any, key, def string, max int) (string, error) {
	v, ok := m[key]
	if !ok {
		return def, nil
	}
	s, ok := v.(string)
	if !ok || len(s) > max {
		return "", fmt.Errorf("invalid %s", key)
	}
	return s, nil
}
func boolean(m map[string]any, key string, def bool) (bool, error) {
	v, ok := m[key]
	if !ok {
		return def, nil
	}
	b, ok := v.(bool)
	if !ok {
		return false, fmt.Errorf("invalid %s", key)
	}
	return b, nil
}
func object(m map[string]any, key string) (map[string]any, error) {
	v, ok := m[key]
	if !ok {
		return map[string]any{}, nil
	}
	o, ok := v.(map[string]any)
	if !ok {
		return nil, fmt.Errorf("invalid %s", key)
	}
	return o, alias(o)
}
func Normalize(kind string, data []byte) (*Request, error) {
	var m map[string]any
	d := json.NewDecoder(bytes.NewReader(data))
	if err := d.Decode(&m); err != nil {
		return nil, err
	}
	if m == nil {
		return nil, errors.New("object required")
	}
	if d.Decode(new(any)) != io.EOF {
		return nil, errors.New("one JSON object required")
	}
	if err := alias(m); err != nil {
		return nil, err
	}
	r := &Request{Protocol: Protocol, Kind: kind, Compatibility: Versions[kind], Params: map[string]any{}, Fallback: "error"}
	p := r.Params
	var err error
	setn := func(src map[string]any, key string, def, lo, hi float64, integer bool) {
		if err != nil {
			return
		}
		p[key], err = number(src, key, def, lo, hi, integer)
	}
	sets := func(src map[string]any, key, def string, max int) {
		if err != nil {
			return
		}
		p[key], err = str(src, key, def, max)
	}
	switch kind {
	case "orbital", "electrostatics":
		for _, k := range []string{"no_com", "no_reorient"} {
			if v, exists := m[k]; exists && v != true {
				return nil, fmt.Errorf("%s must be true", k)
			}
		}
		fields := "sdfContent method charge multiplicity inchiKey clientCacheKey smiles"
		if kind == "orbital" {
			fields += " basis orbitals gridSpacing isovalue densityIsovalue"
		} else {
			fields += " surface potential no_com no_reorient"
		}
		if err = allowed(m, fields); err != nil {
			return nil, err
		}
		r.Input, err = str(m, "sdfContent", "", MaxInput)
		if err != nil {
			return nil, err
		}
		if err = ValidateSDF(r.Input); err != nil {
			return nil, fmt.Errorf("SDF: %w", err)
		}
		setn(m, "charge", 0, -10, 10, true)
		maxSpin := float64(10)
		if kind == "orbital" {
			maxSpin = 1
		}
		setn(m, "multiplicity", 1, 1, maxSpin, true)
		if kind == "orbital" {
			sets(m, "method", "scf", 40)
			sets(m, "basis", "sto-3g", 100)
			setn(m, "gridSpacing", .25, .1, 1, false)
			setn(m, "isovalue", .05, .001, .5, false)
			setn(m, "densityIsovalue", .01, .001, .5, false)
			v, ok := m["orbitals"]
			if !ok {
				v = []any{"homo", "lumo"}
			}
			a, ok := v.([]any)
			if !ok || len(a) == 0 || len(a) > 5 {
				return nil, errors.New("invalid orbitals")
			}
			set := map[string]bool{}
			for _, x := range a {
				s, ok := x.(string)
				s = strings.ToLower(s)
				if !ok || !strings.Contains("|homo|homo-1|lumo|lumo+1|density|", "|"+s+"|") {
					return nil, errors.New("unsupported orbital")
				}
				set[s] = true
			}
			names := []string{}
			for s := range set {
				names = append(names, s)
			}
			sort.Strings(names)
			a = []any{}
			for _, s := range names {
				a = append(a, s)
			}
			p["orbitals"] = a
		} else {
			sets(m, "method", "xtb-gfn2", 40)
			surface, e := object(m, "surface")
			if e != nil {
				return nil, e
			}
			if e = allowed(surface, "gridSpacing probeRadius padding"); e != nil {
				return nil, e
			}
			setn(surface, "gridSpacing", .25, .1, 1, false)
			setn(surface, "probeRadius", 1.4, 0, 5, false)
			setn(surface, "padding", 3, 1, 10, false)
			potential, e := object(m, "potential")
			if e != nil {
				return nil, e
			}
			if e = allowed(potential, "units softeningEpsilon clampPercentiles"); e != nil {
				return nil, e
			}
			sets(potential, "units", "kcal/mol/e", 20)
			setn(potential, "softeningEpsilon", .3, .1, 2, false)
			v, ok := potential["clampPercentiles"]
			if !ok {
				v = []any{float64(5), float64(95)}
			}
			a, ok := v.([]any)
			if !ok || len(a) != 2 {
				return nil, errors.New("two clamp percentiles required")
			}
			lo, lok := a[0].(float64)
			hi, hok := a[1].(float64)
			if !lok || !hok || lo < 0 || hi > 100 || lo >= hi {
				return nil, errors.New("invalid clamp percentiles")
			}
			p["clampPercentiles"] = a
		}
	case "conformer":
		if err = allowed(m, "molblock2d params waitMs geomVersion seed forceField optimize numConformers clientCacheKey smiles inchiKey"); err != nil {
			return nil, err
		}
		r.Input, err = str(m, "molblock2d", "", MaxInput)
		if err != nil {
			return nil, err
		}
		if err = ValidateSDF(r.Input); err != nil {
			return nil, fmt.Errorf("SDF: %w", err)
		}
		opts, e := object(m, "params")
		if e != nil {
			return nil, e
		}
		if e = allowed(opts, "method opt maxAttempts maxOptIters addHs seed"); e != nil {
			return nil, e
		}
		n, e := number(m, "numConformers", 1, 1, 1, true)
		if e != nil || n != 1 {
			return nil, errors.New("only one conformer supported")
		}
		if seed, ok := m["seed"]; ok {
			if old, exists := opts["seed"]; exists && canonical(seed) != canonical(old) {
				return nil, errors.New("conflicting seed")
			}
			opts["seed"] = seed
		}
		_, hasOpt := m["optimize"]
		_, hasFF := m["forceField"]
		if hasOpt || hasFF {
			b, e := boolean(m, "optimize", true)
			if e != nil {
				return nil, e
			}
			ff, e := str(m, "forceField", "mmff", 10)
			if e != nil {
				return nil, e
			}
			ff = strings.ToLower(ff)
			if !b {
				ff = "none"
			}
			if old, ok := opts["opt"]; ok && old != ff {
				return nil, errors.New("conflicting force field")
			}
			opts["opt"] = ff
		}
		if method, ok := opts["method"].(string); ok {
			switch strings.ToLower(method) {
			case "etkdgv3", "etkdg3":
				opts["method"] = "etkdg_v3"
			case "etkdgv2", "etkdg2":
				opts["method"] = "etkdg_v2"
			}
		}
		if opt, ok := opts["opt"].(string); ok {
			opt = strings.ToLower(opt)
			if opt == "off" || opt == "false" || opt == "no" {
				opt = "none"
			}
			opts["opt"] = opt
		}
		sets(opts, "method", "etkdg_v3", 30)
		sets(opts, "opt", "uff", 10)
		sets(m, "geomVersion", "", 100)
		if p["geomVersion"] == "" {
			return nil, errors.New("geomVersion required")
		}
		setn(opts, "seed", 0, 0, 2147483647, true)
		setn(opts, "maxAttempts", 10, 1, 100, true)
		setn(opts, "maxOptIters", 200, 0, 2000, true)
		b, e := boolean(opts, "addHs", false)
		if e != nil {
			return nil, e
		}
		p["addHs"] = b
		wait, e := number(m, "waitMs", 0, 0, 10000, true)
		if e != nil {
			return nil, e
		}
		r.WaitMS = int(wait)
		r.Fallback = "preserve-original-planar"
	case "nomenclature":
		if err = allowed(m, "smiles molblock verifyOpsin"); err != nil {
			return nil, err
		}
		r.Input, err = str(m, "smiles", "", 4096)
		if err != nil || len(r.Input) == 0 {
			return nil, errors.New("smiles required")
		}
		if m["molblock"] == nil {
			delete(m, "molblock")
		}
		sets(m, "molblock", "", MaxInput)
		b, e := boolean(m, "verifyOpsin", true)
		if e != nil {
			return nil, e
		}
		p["verifyOpsin"] = b
	default:
		return nil, errors.New("unsupported operation")
	}
	if err != nil {
		return nil, err
	}
	if method, ok := p["method"].(string); ok {
		method = strings.ToLower(method)
		p["method"] = method
		allowedMethods := map[string]string{"orbital": "|scf|hf|b3lyp|pbe|pbe0|mp2|", "electrostatics": "|xtb-gfn2|gasteiger|", "conformer": "|etkdg_v3|etkdg_v2|etkdg|"}
		if (kind == "orbital" && !regexp.MustCompile(`^[a-z0-9_-]{1,40}$`).MatchString(method)) || (kind != "orbital" && !strings.Contains(allowedMethods[kind], "|"+method+"|")) {
			return nil, errors.New("unsupported method")
		}
	}
	if basis, ok := p["basis"].(string); ok {
		basis = strings.ToLower(basis)
		if len(basis) == 0 {
			return nil, errors.New("basis required")
		}
		for _, c := range basis {
			if !(c >= 'a' && c <= 'z' || c >= '0' && c <= '9' || strings.ContainsRune("-+*(),_", c)) {
				return nil, errors.New("invalid basis")
			}
		}
		p["basis"] = basis
	}
	if kind == "electrostatics" {
		if p["units"] != "kcal/mol/e" {
			return nil, errors.New("unsupported potential units")
		}
		if p["method"] == "gasteiger" && p["multiplicity"] != float64(1) {
			return nil, errors.New("Gasteiger open shell unsupported")
		}
	}
	if kind == "conformer" && !strings.Contains("|none|uff|mmff|", "|"+p["opt"].(string)+"|") {
		return nil, errors.New("unsupported optimization")
	}
	r.InputSHA256 = digest([]byte(r.Input))
	r.Key = r.Identity()
	return r, nil
}

const Elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"

// ValidateSDF is structural only. RDKit/Psi4 remain chemical authorities.
func ValidateSDF(s string) error {
	if len(s) == 0 || len(s) > MaxInput {
		return errors.New("invalid SDF size")
	}
	lines := strings.Split(strings.ReplaceAll(strings.ReplaceAll(s, "\r\n", "\n"), "\r", "\n"), "\n")
	if len(lines) < 5 || len(lines[3]) < 6 {
		return errors.New("SDF counts line missing")
	}
	if !strings.Contains(lines[3], "V2000") {
		return errors.New("only V2000 is supported")
	}
	n, err := strconv.Atoi(strings.TrimSpace(lines[3][:3]))
	if err != nil || n < 1 || n > 999 {
		return errors.New("invalid atom count")
	}
	b, err := strconv.Atoi(strings.TrimSpace(lines[3][3:6]))
	if err != nil || b < 0 || b > 999 || len(lines) < 4+n+b+1 {
		return errors.New("invalid bond count or truncated SDF")
	}
	for _, line := range lines[4 : 4+n] {
		if len(line) < 34 {
			return errors.New("truncated atom record")
		}
		for i := 0; i < 3; i++ {
			x, e := strconv.ParseFloat(strings.TrimSpace(line[i*10:i*10+10]), 64)
			if e != nil || math.IsNaN(x) || math.IsInf(x, 0) || math.Abs(x) > 1e6 {
				return errors.New("invalid atom coordinate")
			}
		}
		element := strings.TrimSpace(line[31:34])
		if !strings.Contains(" "+Elements+" ", " "+element+" ") {
			return errors.New("unknown element")
		}
	}
	for _, line := range lines[4+n : 4+n+b] {
		if len(line) < 9 {
			return errors.New("truncated bond record")
		}
		for i := 0; i < 2; i++ {
			x, e := strconv.Atoi(strings.TrimSpace(line[i*3 : i*3+3]))
			if e != nil || x < 1 || x > n {
				return errors.New("invalid bond atom index")
			}
		}
	}
	end := false
	for _, line := range lines[4+n+b:] {
		if strings.TrimSpace(line) == "M  END" {
			end = true
			break
		}
	}
	if !end {
		return errors.New("missing M END")
	}
	return nil
}
