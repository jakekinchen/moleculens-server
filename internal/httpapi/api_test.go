package httpapi

import (
	"bytes"
	"context"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/jackc/pgx/v5"
	"github.com/jackc/pgx/v5/pgxpool"
	"github.com/jakekinchen/moleculens-server/internal/contract"
	"github.com/jakekinchen/moleculens-server/internal/engine"
	"github.com/jakekinchen/moleculens-server/internal/store"
)

func testAPI(t *testing.T) *API {
	t.Helper()
	dsn := os.Getenv("TEST_DATABASE_URL")
	if dsn == "" {
		t.Skip("real PostgreSQL requires TEST_DATABASE_URL")
	}
	ctx := context.Background()
	admin, e := store.Open(ctx, dsn)
	if e != nil {
		t.Fatal(e)
	}
	name := "mlhttp_" + store.ID()
	if _, e = admin.Pool.Exec(ctx, "CREATE DATABASE "+pgx.Identifier{name}.Sanitize()+" WITH TEMPLATE template0 ENCODING 'UTF8'"); e != nil {
		t.Fatal(e)
	}
	cfg, e := pgxpool.ParseConfig(dsn)
	if e != nil {
		t.Fatal(e)
	}
	cfg.ConnConfig.Database = name
	pool, e := pgxpool.NewWithConfig(ctx, cfg)
	if e != nil {
		t.Fatal(e)
	}
	s := &store.Store{Pool: pool}
	if e = s.Migrate(ctx); e != nil {
		t.Fatal(e)
	}
	t.Cleanup(func() {
		s.Pool.Close()
		admin.Pool.Exec(ctx, "DROP DATABASE "+pgx.Identifier{name}.Sanitize()+" WITH(FORCE)")
		admin.Pool.Close()
	})
	return &API{Store: s, Root: t.TempDir(), LegacyRoot: t.TempDir(), Capacity: 10}
}
func httpCall(t *testing.T, h http.Handler, method, path string, body any) (int, map[string]any) {
	t.Helper()
	b, _ := json.Marshal(body)
	rec := httptest.NewRecorder()
	h.ServeHTTP(rec, httptest.NewRequest(method, path, bytes.NewReader(b)))
	var result map[string]any
	if e := json.Unmarshal(rec.Body.Bytes(), &result); e != nil {
		t.Fatal(e, rec.Body.String())
	}
	return rec.Code, result
}
func water(t *testing.T) string {
	b, e := os.ReadFile("../../tests/contracts/water-v2000.mol")
	if e != nil {
		t.Fatal(e)
	}
	return string(b)
}
func TestImmediateValidationHealthAndKindSafety(t *testing.T) {
	a := testAPI(t)
	h := a.Handler()
	code, _ := httpCall(t, h, "POST", "/v1/orbitals/compute", map[string]any{"sdfContent": "bad", "clientCacheKey": "override"})
	if code != 400 {
		t.Fatalf("SDF status %d", code)
	}
	code, _ = httpCall(t, h, "GET", "/ready", nil)
	if code != 503 {
		t.Fatalf("no worker ready: %d", code)
	}
	code, r := httpCall(t, h, "POST", "/v1/conformers/compute", map[string]any{"molblock2d": water(t), "geomVersion": "v1"})
	if code != 200 {
		t.Fatal(code, r)
	}
	id := r["jobId"].(string)
	code, _ = httpCall(t, h, "GET", "/v1/orbitals/jobs/"+id, nil)
	if code != 404 {
		t.Fatal("cross-kind", code)
	}
	j, e := a.Store.Get(context.Background(), id, "conformer")
	if e != nil || j.State != "pending" {
		t.Fatal(j, e)
	}
}
func TestCutoverRollbackRetainsIDsArtifactsAndQuiesces(t *testing.T) {
	a := testAPI(t)
	ctx := context.Background()
	req, e := contract.Normalize("conformer", mustJSON(map[string]any{"molblock2d": water(t), "geomVersion": "v1"}))
	if e != nil {
		t.Fatal(e)
	}
	j, _, e := a.Store.Admit(ctx, req.Kind, req.Key, mustJSON(req), 10)
	if e != nil {
		t.Fatal(e)
	}
	attempt, e := a.Store.Claim(ctx, "worker", time.Minute)
	if e != nil {
		t.Fatal(e)
	}
	root := filepath.Join(a.Root, j.Key, *attempt.Token)
	os.MkdirAll(root, 0700)
	payload := mustJSON(map[string]any{"sdf3d": water(t), "meta": map[string]any{"geomVersion": "v1", "method": "etkdg_v3"}, "timingMs": map[string]any{"total": 1}})
	sum := sha256.Sum256(payload)
	m := engine.Manifest{Protocol: 1, Kind: j.Kind, Key: j.Key, Token: *attempt.Token, Payload: "response.json", Files: []engine.File{{Name: "response.json", Size: int64(len(payload)), SHA256: hex.EncodeToString(sum[:])}}}
	os.WriteFile(filepath.Join(root, "response.json"), payload, 0600)
	if e = a.Store.Complete(ctx, attempt, mustJSON(m)); e != nil {
		t.Fatal(e)
	}
	// Saturated readers must not report lost artifacts or invalidate a published result.
	func() {
		resultReads <- struct{}{}
		resultReads <- struct{}{}
		defer func() { <-resultReads; <-resultReads }()
		for _, path := range []string{"/v1/conformers/jobs/" + j.ID, "/v1/conformers/artifacts/" + j.Key + "/response.json"} {
			code, response := httpCall(t, a.Handler(), "GET", path, nil)
			if code != 503 || response["detail"].(map[string]any)["code"] != "busy" {
				t.Fatalf("capacity was misreported: %d %v", code, response)
			}
		}
		historical, err := a.Store.Get(ctx, j.ID, j.Kind)
		if err != nil || historical.State != "done" || !historical.Available {
			t.Fatalf("capacity changed publication: %v %v", historical, err)
		}
	}()
	// A fresh API process/handler represents both cutover and fenced rollback.
	for _, quiesced := range []bool{false, true} {
		next := &API{Store: a.Store, Root: a.Root, LegacyRoot: a.LegacyRoot, Capacity: 10, Quiesced: quiesced}
		code, r := httpCall(t, next.Handler(), "GET", "/v1/conformers/jobs/"+j.ID, nil)
		if code != 200 || r["status"] != "done" || r["sdf3d"] != water(t) {
			t.Fatal(code, r)
		}
		if quiesced {
			code, _ = httpCall(t, next.Handler(), "POST", "/v1/conformers/compute", map[string]any{"molblock2d": water(t), "geomVersion": "v1"})
			if code != 503 {
				t.Fatal("rollback admits new work", code)
			}
		}
	}
	os.Remove(filepath.Join(root, "response.json"))
	code, _ := httpCall(t, a.Handler(), "GET", "/v1/conformers/jobs/"+j.ID, nil)
	if code != 503 {
		t.Fatal("missing result", code)
	}
	historical, e := a.Store.Get(ctx, j.ID, j.Kind)
	if e != nil || historical.State != "done" {
		t.Fatal("history changed", historical, e)
	}
}
func TestLegacyNamespaceReadOnlyAndWorkerIsolation(t *testing.T) {
	a := testAPI(t)
	ctx := context.Background()
	meta := map[string]any{"quality": "embedded", "method": "etkdg_v3", "geom_version": "legacy"}
	_, e := a.Store.Pool.Exec(ctx, `INSERT INTO public.jobs(id,cache_key,status,request_params,result_meta) VALUES('legacy-id','legacy-key','done',$1,$2)`, mustJSON(map[string]any{"job_type": "conformer"}), mustJSON(meta))
	if e != nil {
		t.Fatal(e)
	}
	root := filepath.Join(a.LegacyRoot, "legacy-key")
	os.Mkdir(root, 0700)
	os.WriteFile(filepath.Join(root, "structure.sdf"), []byte(water(t)), 0600)
	code, r := httpCall(t, a.Handler(), "GET", "/v1/conformers/jobs/legacy-id", nil)
	if code != 200 || r["status"] != "done" {
		t.Fatal(code, r)
	}
	code, _ = httpCall(t, a.Handler(), "GET", "/v1/orbitals/jobs/legacy-id", nil)
	if code != 404 {
		t.Fatal(code)
	}
	_, _, e = a.Store.Admit(ctx, "orbital", "new-key", []byte(`{}`), 10)
	if e != nil {
		t.Fatal(e)
	}
	var n int
	if e = a.Store.Pool.QueryRow(ctx, `SELECT count(*) FROM public.jobs WHERE status='pending'`).Scan(&n); e != nil || n != 0 {
		t.Fatal("legacy broad claim can see Go queue", n, e)
	}
	old, e := a.Store.Legacy(ctx, "legacy-id", "conformer")
	if e != nil || old.State != "done" {
		t.Fatal(old, e)
	}
}
func mustJSON(v any) []byte {
	b, e := json.Marshal(v)
	if e != nil {
		panic(e)
	}
	return b
}

func TestOrphanCollectionPreservesPublishedAndAmbiguous(t *testing.T) {
	a := testAPI(t)
	ctx := context.Background()
	key := strings.Repeat("a", 64)
	published := strings.Repeat("1", 32)
	retired := strings.Repeat("2", 32)
	ambiguous := strings.Repeat("3", 32)
	fresh := strings.Repeat("4", 32)
	for _, token := range []string{published, retired, ambiguous, fresh} {
		path := filepath.Join(a.Root, key, token)
		if e := os.MkdirAll(path, 0700); e != nil {
			t.Fatal(e)
		}
		os.WriteFile(filepath.Join(path, "bytes"), []byte("retain unless eligible"), 0600)
		if token == ambiguous {
			os.WriteFile(filepath.Join(path, "publication-uncertain"), []byte("job"), 0600)
		}
		if token != fresh {
			old := time.Now().Add(-25 * time.Hour)
			os.Chtimes(path, old, old)
		}
	}
	_, e := a.Store.Pool.Exec(ctx, `INSERT INTO moleculens_v2.jobs(id,kind,cache_key,request,state,token,manifest) VALUES('published','orbital',$1,'{}','done',$2,'{}')`, key, published)
	if e != nil {
		t.Fatal(e)
	}
	candidates, e := engine.CollectOrphans(ctx, a.Store, a.Root, false)
	if e != nil || len(candidates) != 1 || candidates[0] != key+"/"+retired {
		t.Fatal(candidates, e)
	}
	if _, e = os.Stat(filepath.Join(a.Root, key, retired)); e != nil {
		t.Fatal("dry run deleted bytes")
	}
	if _, e = engine.CollectOrphans(ctx, a.Store, a.Root, true); e != nil {
		t.Fatal(e)
	}
	if _, e = os.Stat(filepath.Join(a.Root, key, retired)); !os.IsNotExist(e) {
		t.Fatal("eligible orphan retained", e)
	}
	for _, token := range []string{published, ambiguous, fresh} {
		if _, e = os.Stat(filepath.Join(a.Root, key, token)); e != nil {
			t.Fatal("protected directory removed", token, e)
		}
	}
}
