package engine

import (
	"testing"

	"github.com/jakekinchen/moleculens-server/internal/contract"
)

func TestWorkerKindSelection(t *testing.T) {
	all, err := workerVersions(nil)
	if err != nil || len(all) != len(contract.Versions) {
		t.Fatalf("default worker lost operations: %v %v", all, err)
	}
	short, err := workerVersions([]string{"conformer", " nomenclature "})
	if err != nil || len(short) != 2 || short["nomenclature"] != contract.Versions["nomenclature"] || short["orbital"] != "" {
		t.Fatalf("short worker must not claim orbital work: %v %v", short, err)
	}
	for _, kinds := range [][]string{{""}, {"orbital", "typo"}} {
		if _, err := workerVersions(kinds); err == nil {
			t.Fatalf("invalid selection accepted: %v", kinds)
		}
	}
	short["conformer"] = "changed"
	if contract.Versions["conformer"] == "changed" {
		t.Fatal("worker selection mutated global compatibility")
	}
}
