package main

import (
	"context"
	"flag"
	"fmt"
	"log/slog"
	"net/http"
	"os"
	"os/signal"
	"path/filepath"
	"strings"
	"syscall"
	"time"

	"github.com/jakekinchen/moleculens-server/internal/engine"
	"github.com/jakekinchen/moleculens-server/internal/httpapi"
	"github.com/jakekinchen/moleculens-server/internal/store"
)

func env(k, def string) string {
	if s := os.Getenv(k); s != "" {
		return s
	}
	return def
}
func run() error {
	if len(os.Args) < 2 {
		return fmt.Errorf("usage: moleculens migrate|api|worker|gc")
	}
	mode := os.Args[1]
	fs := flag.NewFlagSet(mode, flag.ContinueOnError)
	addr := fs.String("listen", env("LISTEN_ADDR", ":8000"), "HTTP listen address")
	root := fs.String("artifacts", env("ARTIFACT_ROOT", "var/go-artifacts"), "new artifact root; separate from legacy")
	legacy := fs.String("legacy-artifacts", env("LEGACY_ARTIFACT_ROOT", "var/cache"), "read-only legacy artifact root")
	python := fs.String("python", env("ENGINE_PYTHON", "python"), "chemistry Python executable")
	pythonPath := fs.String("python-path", env("ENGINE_PYTHONPATH", "src"), "chemistry module source path")
	capacity := fs.Int("capacity", 100, "maximum pending plus running jobs")
	lease := fs.Duration("lease", 30*time.Second, "attempt lease")
	timeout := fs.Duration("timeout", 5*time.Minute, "maximum attempt duration")
	attemptDiskMiB := fs.Int64("attempt-disk-mib", 1024, "private attempt disk budget in MiB (1024..16384)")
	applyGC := fs.Bool("apply", false, "delete eligible orphan attempts; gc defaults to dry run")
	quiesced := fs.Bool("quiesced", false, "serve existing jobs but refuse new admissions")
	kinds := fs.String("kinds", "", "worker operation kinds, comma-separated; empty serves all")
	if err := fs.Parse(os.Args[2:]); err != nil {
		return err
	}
	if *capacity < 1 || *capacity > 10000 {
		return fmt.Errorf("capacity must be 1..10000")
	}
	if *timeout < time.Second || *timeout > 2*time.Hour {
		return fmt.Errorf("timeout must be 1s..2h")
	}
	if *attemptDiskMiB < 1024 || *attemptDiskMiB > 16384 {
		return fmt.Errorf("attempt-disk-mib must be 1024..16384")
	}
	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer stop()
	dsn := os.Getenv("DATABASE_URL")
	if dsn == "" {
		return fmt.Errorf("DATABASE_URL required")
	}
	s, err := store.Open(ctx, dsn)
	if err != nil {
		return err
	}
	defer s.Pool.Close()
	if mode == "migrate" {
		return s.Migrate(ctx)
	}
	if err = s.Check(ctx); err != nil {
		return fmt.Errorf("run explicit migrate first: %w", err)
	}
	*root, err = filepath.Abs(*root)
	if err != nil {
		return err
	}
	*legacy, err = filepath.Abs(*legacy)
	if err != nil {
		return err
	}
	if *root == *legacy {
		return fmt.Errorf("new and legacy artifact roots must differ")
	}
	*pythonPath, err = filepath.Abs(*pythonPath)
	if err != nil {
		return err
	}
	switch mode {
	case "gc":
		paths, e := engine.CollectOrphans(ctx, s, *root, *applyGC)
		for _, p := range paths {
			fmt.Println(p)
		}
		return e
	case "worker":
		var selected []string
		if *kinds != "" {
			selected = strings.Split(*kinds, ",")
		}
		return (engine.Runner{Python: *python, PythonPath: *pythonPath, Root: *root, Build: env("ENGINE_BUILD", "local-unqualified"), Path: os.Getenv("PATH"), Timeout: *timeout, Kinds: selected, AttemptDisk: *attemptDiskMiB << 20}).Worker(ctx, s, *lease)
	case "api":
		a := &httpapi.API{Store: s, Root: *root, LegacyRoot: *legacy, Capacity: *capacity, Origins: strings.Split(env("ALLOWED_ORIGINS", "http://localhost:3000"), ","), Quiesced: *quiesced}
		server := &http.Server{Addr: *addr, Handler: a.Handler(), ReadHeaderTimeout: 5 * time.Second, ReadTimeout: 15 * time.Second, WriteTimeout: 20 * time.Second, IdleTimeout: 30 * time.Second, MaxHeaderBytes: 16 << 10}
		done := make(chan error, 1)
		go func() { slog.Info("HTTP service listening", "address", *addr); done <- server.ListenAndServe() }()
		select {
		case err = <-done:
			if err == http.ErrServerClosed {
				return nil
			}
			return err
		case <-ctx.Done():
			shutdown, cancel := context.WithTimeout(context.Background(), 20*time.Second)
			defer cancel()
			return server.Shutdown(shutdown)
		}
	default:
		return fmt.Errorf("unknown command %s", mode)
	}
}
func main() {
	if err := run(); err != nil {
		slog.Error("moleculens stopped", "error", err)
		os.Exit(1)
	}
}
