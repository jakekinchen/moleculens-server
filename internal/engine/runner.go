package engine

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io/fs"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"syscall"
	"time"

	"github.com/jackc/pgx/v5"
	"github.com/jakekinchen/moleculens-server/internal/contract"
	"github.com/jakekinchen/moleculens-server/internal/store"
)

type ProcessFailure struct {
	Err       error
	Retryable bool
}

// Temporary engine scratch has its own bound; published artifacts stay smaller.
const MaxAttemptDisk = 1 << 30
const MaxConfiguredAttemptDisk = 16 << 30

func (e *ProcessFailure) Error() string { return "adapter process: " + e.Err.Error() }

type Runner struct {
	Python      string
	PythonPath  string
	Root        string
	Build       string
	Path        string
	Timeout     time.Duration
	Kinds       []string
	AttemptDisk int64
}

func (r Runner) attemptDisk() int64 {
	if r.AttemptDisk == 0 {
		return MaxAttemptDisk
	}
	return r.AttemptDisk
}

func workerVersions(kinds []string) (map[string]string, error) {
	selected := make(map[string]string)
	if len(kinds) == 0 {
		for kind, version := range contract.Versions {
			selected[kind] = version
		}
		return selected, nil
	}
	for _, value := range kinds {
		kind := strings.TrimSpace(value)
		version, ok := contract.Versions[kind]
		if !ok {
			return nil, fmt.Errorf("unknown worker kind %q", value)
		}
		selected[kind] = version
	}
	return selected, nil
}

type cappedLog struct {
	mu sync.Mutex
	f  *os.File
	n  int
}

func (w *cappedLog) Write(p []byte) (int, error) {
	w.mu.Lock()
	defer w.mu.Unlock()
	n := len(p)
	if w.n < 1<<20 {
		take := n
		if take > (1<<20)-w.n {
			take = (1 << 20) - w.n
		}
		_, e := w.f.Write(p[:take])
		if e != nil {
			return 0, e
		}
		w.n += take
	}
	return n, nil
}
func (r Runner) environment(root string) []string {
	return []string{"PATH=" + r.Path, "PYTHONPATH=" + r.PythonPath, "HOME=" + root, "TMPDIR=" + root, "LANG=C.UTF-8", "PYTHONNOUSERSITE=1", "PYTHONDONTWRITEBYTECODE=1", "ENGINE_BUILD=" + r.Build, "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "PSI4_NUM_THREADS=1", "PSI4_MEMORY_MB=1024"}
}
func (r Runner) command(ctx context.Context, dir string, args ...string) error {
	cmd := exec.Command(r.Python, append([]string{"-m", "moleculens.adapter"}, args...)...)
	cmd.Env = r.environment(dir)
	cmd.Dir = dir
	cmd.SysProcAttr = processAttributes()
	log, err := os.OpenFile(filepath.Join(dir, "adapter.log"), os.O_CREATE|os.O_EXCL|os.O_WRONLY, 0600)
	if err != nil {
		return err
	}
	defer log.Close()
	out := &cappedLog{f: log}
	cmd.Stdout = out
	cmd.Stderr = out
	if err = cmd.Start(); err != nil {
		return err
	}
	done := make(chan error, 1)
	go func() { done <- cmd.Wait() }()
	ticker := time.NewTicker(500 * time.Millisecond)
	defer ticker.Stop()
	// Always kill descendants before publishing, including after a successful leader exit.
	defer syscall.Kill(-cmd.Process.Pid, syscall.SIGKILL)
	for {
		select {
		case err = <-done:
			return err
		case <-ctx.Done():
			syscall.Kill(-cmd.Process.Pid, syscall.SIGKILL)
			<-done
			return ctx.Err()
		case <-ticker.C:
			var size int64
			walkErr := filepath.WalkDir(dir, func(path string, d fs.DirEntry, e error) error {
				if e != nil {
					return e
				}
				if d.Type()&os.ModeSymlink != 0 {
					return errors.New("adapter symlink prohibited")
				}
				if d.IsDir() {
					return nil
				}
				info, e := d.Info()
				if e != nil {
					return e
				}
				size += info.Size()
				if size > r.attemptDisk() {
					return errors.New("attempt disk budget exceeded")
				}
				return nil
			})
			if walkErr != nil {
				syscall.Kill(-cmd.Process.Pid, syscall.SIGKILL)
				<-done
				return walkErr
			}
		}
	}
}
func (r Runner) Probe(ctx context.Context) ([]byte, error) {
	dir, err := os.MkdirTemp(r.Root, "probe-")
	if err != nil {
		return nil, err
	}
	defer os.RemoveAll(dir)
	err = r.command(ctx, dir, "--probe", filepath.Join(dir, "probe.json"))
	if err != nil {
		return nil, err
	}
	return readBounded(filepath.Join(dir, "probe.json"), 128<<10)
}
func (r Runner) Execute(ctx context.Context, req *contract.Request, token string) (*Manifest, error) {
	dir := filepath.Join(r.Root, req.Key, token)
	if err := os.MkdirAll(filepath.Dir(dir), 0700); err != nil {
		return nil, err
	}
	if err := os.Mkdir(dir, 0700); err != nil {
		return nil, err
	}
	b, err := json.Marshal(req)
	if err != nil {
		return nil, err
	}
	input := filepath.Join(dir, "request.json")
	if err = os.WriteFile(input, b, 0400); err != nil {
		return nil, err
	}
	deadline := r.Timeout
	if deadline <= 0 {
		deadline = 5 * time.Minute
	}
	if req.Kind == "nomenclature" {
		deadline = 8 * time.Second
	}
	ctx, cancel := context.WithTimeout(ctx, deadline)
	defer cancel()
	if err = r.command(ctx, dir, "--request", input, "--output", dir, "--token", token); err != nil {
		detail, e := readBounded(filepath.Join(dir, "error.json"), 4096)
		if e == nil {
			return nil, fmt.Errorf("adapter: %s", detail)
		}
		var exitError *exec.ExitError
		retryable := errors.Is(err, context.Canceled) || errors.As(err, &exitError)
		return nil, &ProcessFailure{Err: err, Retryable: retryable}
	}
	m, _, err := Validate(dir, req, token)
	if err != nil {
		return nil, err
	}
	if err = Seal(dir, m); err != nil {
		return nil, err
	}
	return m, nil
}
func (r Runner) Worker(ctx context.Context, s *store.Store, lease time.Duration) error {
	if r.attemptDisk() < MaxAttemptDisk || r.attemptDisk() > MaxConfiguredAttemptDisk {
		return errors.New("attempt disk budget must be 1..16 GiB")
	}
	selected, err := workerVersions(r.Kinds)
	if err != nil {
		return err
	}
	if lease < 15*time.Second {
		return errors.New("lease must be at least 15 seconds")
	}
	if err := os.MkdirAll(r.Root, 0700); err != nil {
		return err
	}
	c, cancel := context.WithTimeout(ctx, 30*time.Second)
	caps, err := r.Probe(c)
	cancel()
	if err != nil {
		return fmt.Errorf("engine capability probe: %w", err)
	}
	var capmap map[string]string
	if err = json.Unmarshal(caps, &capmap); err != nil {
		return err
	}
	for _, key := range []string{"psi4", "rdkit", "openclatura", "java", "xtb"} {
		if capmap[key] == "" || capmap[key] == "unavailable" {
			return fmt.Errorf("required chemistry capability missing: %s", key)
		}
	}
	versions, _ := json.Marshal(selected)
	capmap["compatibility"] = string(versions)
	caps, _ = json.Marshal(capmap)
	owner := "worker-" + store.ID()
	heartbeatCtx, stop := context.WithCancel(ctx)
	defer stop()
	if err = s.Heartbeat(ctx, owner, caps); err != nil {
		return err
	}
	go func() {
		t := time.NewTicker(5 * time.Second)
		defer t.Stop()
		for {
			select {
			case <-heartbeatCtx.Done():
				return
			case <-t.C:
				c, done := context.WithTimeout(heartbeatCtx, 5*time.Second)
				if err := s.Heartbeat(c, owner, caps); err != nil {
					slog.Warn("worker heartbeat failed", "error", err)
				}
				done()
			}
		}
	}()
	for ctx.Err() == nil {
		if err = checkSpace(r.Root, uint64(r.attemptDisk())+StorageReserve); err != nil {
			select {
			case <-ctx.Done():
				return nil
			case <-time.After(time.Second):
				continue
			}
		}
		if err = s.Reap(ctx, 3); err != nil {
			return err
		}
		j, e := s.ClaimCompatible(ctx, owner, lease, versions)
		if errors.Is(e, pgx.ErrNoRows) {
			select {
			case <-ctx.Done():
				return nil
			case <-time.After(250 * time.Millisecond):
				continue
			}
		}
		if e != nil {
			return e
		}
		var req contract.Request
		if e = json.Unmarshal(j.Request, &req); e != nil || req.Protocol != contract.Protocol || req.Key != j.Key || req.Kind != j.Kind || req.Key != req.Identity() {
			s.Fail(ctx, j, "invalid stored request")
			continue
		}
		attempt, cancelAttempt := context.WithCancel(ctx)
		stopRenew := make(chan struct{})
		renewDone := make(chan struct{})
		go func() {
			defer close(renewDone)
			t := time.NewTicker(lease / 3)
			defer t.Stop()
			for {
				select {
				case <-stopRenew:
					return
				case <-attempt.Done():
					return
				case <-t.C:
					c, done := context.WithTimeout(attempt, 5*time.Second)
					e := s.Renew(c, j, lease)
					done()
					if e != nil {
						cancelAttempt()
						return
					}
				}
			}
		}()
		m, runErr := r.Execute(attempt, &req, *j.Token)
		// Keep renewing through validation/sealing above; publication is immediate.
		if runErr == nil {
			b, _ := json.Marshal(m)
			e = s.Complete(attempt, j, b)
			if e != nil {
				if !errors.Is(e, store.ErrFenced) {
					os.WriteFile(filepath.Join(r.Root, j.Key, *j.Token, "publication-uncertain"), []byte(j.ID), 0400)
				}
				slog.Error("publication unresolved; retain attempt", "job", j.ID, "token", *j.Token, "error", e)
			}
		} else {
			c, done := context.WithTimeout(context.Background(), 5*time.Second)
			var processFailure *ProcessFailure
			if attempt.Err() != nil || (errors.As(runErr, &processFailure) && processFailure.Retryable) {
				e = s.Retry(c, j, runErr.Error(), 3)
			} else {
				e = s.Fail(c, j, runErr.Error())
			}
			done()
			slog.Error("attempt failed", "job", j.ID, "error", runErr, "transitionError", e)
		}
		close(stopRenew)
		cancelAttempt()
		<-renewDone
	}
	return nil
}
