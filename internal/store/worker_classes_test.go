package store

import (
	"context"
	"encoding/json"
	"testing"
	"time"

	"github.com/jakekinchen/moleculens-server/internal/contract"
)

func TestPostgresWorkerClassesPreserveReadinessAndClaimOwnership(t *testing.T) {
	s := testStore(t)
	ctx := context.Background()
	all, _ := json.Marshal(contract.Versions)
	long := map[string]string{"orbital": contract.Versions["orbital"], "electrostatics": contract.Versions["electrostatics"]}
	short := map[string]string{"conformer": contract.Versions["conformer"], "nomenclature": contract.Versions["nomenclature"]}
	heartbeat := func(id string, selected map[string]string) {
		t.Helper()
		versions, _ := json.Marshal(selected)
		caps, _ := json.Marshal(map[string]string{"psi4": "1.11", "compatibility": string(versions)})
		if err := s.Heartbeat(ctx, id, caps); err != nil {
			t.Fatal(err)
		}
	}
	heartbeat("long", long)
	if s.Ready(ctx, string(all)) == nil {
		t.Fatal("missing short worker was considered ready")
	}
	heartbeat("short", short)
	if err := s.Ready(ctx, string(all)); err != nil {
		t.Fatal(err)
	}
	for _, kind := range []string{"orbital", "nomenclature"} {
		request, _ := json.Marshal(map[string]string{"compatibility": contract.Versions[kind]})
		if _, _, err := s.Admit(ctx, kind, kind, request, 10); err != nil {
			t.Fatal(err)
		}
	}
	versions, _ := json.Marshal(short)
	job, err := s.ClaimCompatible(ctx, "short", time.Minute, versions)
	if err != nil || job.Kind != "nomenclature" {
		t.Fatalf("short worker consumed long work or missed naming: %v %v", job, err)
	}
	if _, err := s.Pool.Exec(ctx, `UPDATE moleculens_v2.workers SET heartbeat=now()-interval '31 seconds' WHERE id='short'`); err != nil {
		t.Fatal(err)
	}
	if s.Ready(ctx, string(all)) == nil {
		t.Fatal("stale short worker was considered ready")
	}
	heartbeat("short", map[string]string{"conformer": "wrong-version", "nomenclature": contract.Versions["nomenclature"]})
	if s.Ready(ctx, string(all)) == nil {
		t.Fatal("wrong compatibility was considered ready")
	}
	heartbeat("normal", contract.Versions)
	if err := s.Ready(ctx, string(all)); err != nil {
		t.Fatalf("normal all-kind worker regressed: %v", err)
	}
}
