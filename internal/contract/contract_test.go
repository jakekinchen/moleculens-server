package contract

import (
	"encoding/json"
	"os"
	"strings"
	"testing"
)

func water(t *testing.T) string {
	t.Helper()
	b, e := os.ReadFile("../../tests/golden/water.sdf")
	if e != nil {
		t.Fatal(e)
	}
	return string(b)
}
func request(t *testing.T, kind string, m map[string]any) *Request {
	t.Helper()
	b, e := json.Marshal(m)
	if e != nil {
		t.Fatal(e)
	}
	r, e := Normalize(kind, b)
	if e != nil {
		t.Fatal(e)
	}
	return r
}
func TestIdentityAndAliases(t *testing.T) {
	s := water(t)
	a := request(t, "orbital", map[string]any{"sdfContent": s, "clientCacheKey": "hint"})
	b := request(t, "orbital", map[string]any{"sdf_content": s, "client_cache_key": "different", "grid_spacing": .25})
	if a.Key != b.Key {
		t.Fatal("hints or aliases changed identity")
	}
	changed := []map[string]any{{"sdfContent": strings.Replace(s, "0.7586", "0.7587", 1), "clientCacheKey": "hint"}, {"sdfContent": s, "isovalue": .050001}, {"sdfContent": s, "charge": 2}, {"sdfContent": s, "density_isovalue": .011}, {"sdfContent": s, "basis": "6-31g*"}}
	for _, m := range changed {
		if request(t, "orbital", m).Key == a.Key {
			t.Fatal("changed scientific input collided")
		}
	}
	a.Compatibility += "-new"
	if a.Identity() == b.Key {
		t.Fatal("compatibility not included")
	}
	bad, _ := json.Marshal(map[string]any{"sdfContent": "bad", "clientCacheKey": "trusted"})
	if _, e := Normalize("orbital", bad); e == nil {
		t.Fatal("client hint bypassed validation")
	}
	bad, _ = json.Marshal(map[string]any{"sdfContent": s, "gridSpacing": .25, "grid_spacing": .3})
	if _, e := Normalize("orbital", bad); e == nil {
		t.Fatal("ambiguous aliases accepted")
	}
}
func TestConsumerMappings(t *testing.T) {
	s := water(t)
	r := request(t, "conformer", map[string]any{"molblock2d": s, "geomVersion": "etkdgv3", "optimize": true, "forceField": "mmff", "seed": 42, "numConformers": 1, "clientCacheKey": "hint", "smiles": "O"})
	if r.Params["opt"] != "mmff" || r.Params["seed"] != float64(42) {
		t.Fatal(r.Params)
	}
	e := request(t, "electrostatics", map[string]any{"sdfContent": s, "method": "gasteiger", "no_com": true, "no_reorient": true})
	if e.Params["probeRadius"] != 1.4 || e.Params["padding"] != float64(3) {
		t.Fatal("legacy defaults changed")
	}
}
func TestStructuralSDF(t *testing.T) {
	s := water(t)
	blank := s[strings.Index(s, "\n"):]
	if e := ValidateSDF(blank); e != nil {
		t.Fatal("valid blank title:", e)
	}
	for _, bad := range []string{strings.Replace(s, "V2000", "V3000", 1), strings.Replace(s, " O ", " Xx ", 1), strings.Replace(s, "    0.0000", "       nan", 1), strings.Replace(s, "M  END", "M BAD", 1)} {
		if e := ValidateSDF(bad); e == nil {
			t.Fatal("invalid SDF accepted")
		}
	}
}

func TestIndependentFrozenDescriptors(t *testing.T) {
	b, e := os.ReadFile("../../tests/contracts/requests.json")
	if e != nil {
		t.Fatal(e)
	}
	var fixtures []struct {
		Kind     string          `json:"kind"`
		Request  json.RawMessage `json:"request"`
		Params   map[string]any  `json:"params"`
		Key      string          `json:"key"`
		InputSHA string          `json:"inputSha256"`
	}
	if e = json.Unmarshal(b, &fixtures); e != nil {
		t.Fatal(e)
	}
	for _, f := range fixtures {
		r, e := Normalize(f.Kind, f.Request)
		if e != nil {
			t.Fatal(e)
		}
		if r.Key != f.Key || r.InputSHA256 != f.InputSHA || canonical(r.Params) != canonical(f.Params) {
			t.Fatalf("%s independent descriptor disagrees: key=%s expected=%s params=%v", f.Kind, r.Key, f.Key, r.Params)
		}
	}
}
