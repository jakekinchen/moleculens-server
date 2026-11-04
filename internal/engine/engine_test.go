package engine

import (
	"bytes"
	"compress/gzip"
	"context"
	"crypto/sha256"
	"encoding/base64"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/jakekinchen/moleculens-server/internal/contract"
)

func TestManifestIntegrityAndRequiredPayload(t *testing.T) {
	root := t.TempDir()
	request := &contract.Request{Protocol: 1, Kind: "nomenclature", Key: "key", Compatibility: "name-1", InputSHA256: "input"}
	payload := []byte(`{"displayName":"ethanol","source":"openclatura","verified":true,"trace":[]}`)
	hash := sha256.Sum256(payload)
	m := Manifest{Protocol: 1, Kind: request.Kind, Key: request.Key, Token: "attempt", Compatibility: request.Compatibility, InputSHA256: request.InputSHA256, Payload: "response.json", Files: []File{{Name: "response.json", Size: int64(len(payload)), SHA256: hex.EncodeToString(hash[:])}}, Provenance: map[string]string{"build": "test", "python": "3"}}
	b, _ := json.Marshal(m)
	os.WriteFile(filepath.Join(root, "manifest.json"), b, 0600)
	os.WriteFile(filepath.Join(root, "response.json"), payload, 0600)
	if _, _, e := Validate(root, request, "attempt"); e != nil {
		t.Fatal(e)
	}
	if _, _, e := Validate(root, request, "stale"); e == nil {
		t.Fatal("wrong attempt accepted")
	}
	os.WriteFile(filepath.Join(root, "response.json"), []byte(`{}`), 0600)
	if _, e := Load(root, &m); e == nil {
		t.Fatal("corruption accepted")
	}
	os.Remove(filepath.Join(root, "response.json"))
	os.Symlink(filepath.Join(root, "manifest.json"), filepath.Join(root, "response.json"))
	if _, e := Load(root, &m); e == nil {
		t.Fatal("symlink accepted")
	}
	if e := ValidatePayload(request, []byte(`{}`)); e == nil {
		t.Fatal("empty JSON accepted")
	}
}
func zipped(data []byte) string {
	var b bytes.Buffer
	z := gzip.NewWriter(&b)
	z.Write(data)
	z.Close()
	return base64.StdEncoding.EncodeToString(b.Bytes())
}
func TestMeshRejectsBadBuffersAndBoundedExpansion(t *testing.T) {
	if e := buffer(zipped(make([]byte, 4<<20)), 3, false, 0); e == nil {
		t.Fatal("oversized expansion accepted")
	}
	if e := buffer(zipped([]byte{0, 0, 128, 127}), 1, false, 0); e == nil {
		t.Fatal("infinity accepted")
	}
	if e := buffer(zipped([]byte{3, 0, 0, 0}), 1, true, 3); e == nil {
		t.Fatal("out of range index accepted")
	}
}
func TestCancellationReapsDescendantsAndStripsCredentials(t *testing.T) {
	python, e := exec.LookPath("python3")
	if e != nil {
		t.Skip("python3 unavailable")
	}
	src := t.TempDir()
	pkg := filepath.Join(src, "moleculens", "adapter")
	os.MkdirAll(pkg, 0700)
	script := `import os, subprocess, sys, time
assert 'DATABASE_URL' not in os.environ
subprocess.Popen([sys.executable, '-c', "import time; time.sleep(1); open('late-write','w').write('bad')"])
open('started','w').write('yes')
time.sleep(100)
`
	os.WriteFile(filepath.Join(pkg, "__main__.py"), []byte(script), 0600)
	root := t.TempDir()
	r := Runner{Python: python, PythonPath: src, Root: root, Path: os.Getenv("PATH"), Build: "test"}
	t.Setenv("DATABASE_URL", "secret-should-not-propagate")
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()
	done := make(chan error, 1)
	go func() { done <- r.command(ctx, root) }()
	deadline := time.Now().Add(3 * time.Second)
	for {
		if _, e = os.Stat(filepath.Join(root, "started")); e == nil {
			break
		}
		if time.Now().After(deadline) {
			t.Fatal("helper never started")
		}
		time.Sleep(10 * time.Millisecond)
	}
	cancel()
	select {
	case e = <-done:
		if !errors.Is(e, context.Canceled) {
			t.Fatal(e)
		}
	case <-time.After(3 * time.Second):
		t.Fatal("process not reaped")
	}
	time.Sleep(1200 * time.Millisecond)
	if _, e = os.Stat(filepath.Join(root, "late-write")); e == nil {
		t.Fatal("descendant survived cancellation")
	}
}

func TestAttemptScratchDoesNotExpandPublicationBudget(t *testing.T) {
	python, err := exec.LookPath("python3")
	if err != nil {
		t.Skip("python3 unavailable")
	}
	for _, tc := range []struct {
		name      string
		size      int64
		budget    int64
		wantError bool
	}{{"512-MiB temporary scratch", 512 << 20, 0, false}, {"scratch over default 1 GiB", MaxAttemptDisk + 1, 0, true}, {"configured large workspace", 2 << 30, MaxConfiguredAttemptDisk, false}, {"scratch over configured limit", MaxConfiguredAttemptDisk + 1, MaxConfiguredAttemptDisk, true}} {
		t.Run(tc.name, func(t *testing.T) {
			src := t.TempDir()
			pkg := filepath.Join(src, "moleculens", "adapter")
			if err := os.MkdirAll(pkg, 0700); err != nil {
				t.Fatal(err)
			}
			script := fmt.Sprintf("import time\nwith open('scratch', 'wb') as f:\n    f.truncate(%d)\ntime.sleep(1.2)\n", tc.size)
			if err := os.WriteFile(filepath.Join(pkg, "__main__.py"), []byte(script), 0600); err != nil {
				t.Fatal(err)
			}
			root := t.TempDir()
			r := Runner{Python: python, PythonPath: src, Root: root, Path: os.Getenv("PATH"), Build: "test", AttemptDisk: tc.budget}
			ctx, cancel := context.WithTimeout(context.Background(), 4*time.Second)
			defer cancel()
			err := r.command(ctx, root)
			if tc.wantError {
				if err == nil || !strings.Contains(err.Error(), "disk budget exceeded") {
					t.Fatalf("wanted scratch-budget rejection, got %v", err)
				}
			} else if err != nil {
				t.Fatalf("temporary scratch within its own budget was rejected: %v", err)
			}
		})
	}
	_, err = Load(t.TempDir(), &Manifest{Payload: "response.json", Files: []File{{Name: "response.json", Size: MaxArtifacts + 1}}})
	if err == nil || !strings.Contains(err.Error(), "artifact budget exceeded") {
		t.Fatalf("published artifacts must retain their smaller budget: %v", err)
	}
}
