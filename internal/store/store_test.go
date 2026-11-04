package store

import (
	"context"
	"encoding/json"
	"errors"
	"os"
	"sync"
	"testing"
	"time"

	"github.com/jackc/pgx/v5"
	"github.com/jackc/pgx/v5/pgxpool"
)

func testStore(t *testing.T) *Store {
	t.Helper()
	dsn := os.Getenv("TEST_DATABASE_URL")
	if dsn == "" {
		t.Skip("real PostgreSQL requires TEST_DATABASE_URL with disposable-database creation permission")
	}
	admin, e := Open(context.Background(), dsn)
	if e != nil {
		t.Fatal(e)
	}
	name := "mltest_" + ID()
	if _, e = admin.Pool.Exec(context.Background(), "CREATE DATABASE "+pgx.Identifier{name}.Sanitize()+" WITH TEMPLATE template0 ENCODING 'UTF8'"); e != nil {
		admin.Pool.Close()
		t.Fatal(e)
	}
	cfg, e := pgxpool.ParseConfig(dsn)
	if e != nil {
		t.Fatal(e)
	}
	cfg.ConnConfig.Database = name
	pool, e := pgxpool.NewWithConfig(context.Background(), cfg)
	s := &Store{Pool: pool}
	if e != nil {
		t.Fatal(e)
	}
	t.Cleanup(func() {
		s.Pool.Close()
		admin.Pool.Exec(context.Background(), "DROP DATABASE "+pgx.Identifier{name}.Sanitize()+" WITH (FORCE)")
		admin.Pool.Close()
	})
	if e = s.Migrate(context.Background()); e != nil {
		t.Fatal(e)
	}
	return s
}

func TestPostgresAdmissionContentionAndCapacity(t *testing.T) {
	s := testStore(t)
	ctx := context.Background()
	var wg sync.WaitGroup
	ids := make(chan string, 32)
	errs := make(chan error, 32)
	for range 32 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			j, _, e := s.Admit(ctx, "orbital", "same", []byte(`{}`), 1)
			if e != nil {
				errs <- e
			} else {
				ids <- j.ID
			}
		}()
	}
	wg.Wait()
	close(errs)
	close(ids)
	for e := range errs {
		t.Error(e)
	}
	unique := map[string]bool{}
	for id := range ids {
		unique[id] = true
	}
	if len(unique) != 1 {
		t.Fatalf("created %d active IDs", len(unique))
	}
	if _, _, e := s.Admit(ctx, "orbital", "other", []byte(`{}`), 1); !errors.Is(e, ErrFull) {
		t.Fatalf("capacity: %v", e)
	}
	var n int
	if e := s.Pool.QueryRow(ctx, `SELECT count(*) FROM moleculens_v2.jobs`).Scan(&n); e != nil || n != 1 {
		t.Fatalf("rows=%d error=%v", n, e)
	}
}
func TestPostgresClaimsFencingAndAbsorbingTerminal(t *testing.T) {
	s := testStore(t)
	ctx := context.Background()
	job, _, e := s.Admit(ctx, "conformer", "key", []byte(`{}`), 20)
	if e != nil {
		t.Fatal(e)
	}
	var wg sync.WaitGroup
	claims := make(chan *Job, 16)
	for range 16 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			j, e := s.Claim(ctx, ID(), time.Minute)
			if e == nil {
				claims <- j
			} else if !errors.Is(e, pgx.ErrNoRows) {
				t.Error(e)
			}
		}()
	}
	wg.Wait()
	close(claims)
	var first *Job
	count := 0
	for j := range claims {
		first = j
		count++
	}
	if count != 1 {
		t.Fatalf("claims=%d", count)
	}
	_, e = s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET lease_until=now()-interval '1 second' WHERE id=$1`, job.ID)
	if e != nil {
		t.Fatal(e)
	}
	if e = s.Reap(ctx, 3); e != nil {
		t.Fatal(e)
	}
	second, e := s.Claim(ctx, "next", time.Minute)
	if e != nil {
		t.Fatal(e)
	}
	if *first.Token == *second.Token {
		t.Fatal("token reused")
	}
	if e = s.Complete(ctx, first, []byte(`{"old":true}`)); !errors.Is(e, ErrFenced) {
		t.Fatalf("late completion: %v", e)
	}
	if e = s.Fail(ctx, first, "late"); !errors.Is(e, ErrFenced) {
		t.Fatalf("late failure: %v", e)
	}
	manifest := []byte(`{"new":true}`)
	if e = s.Complete(ctx, second, manifest); e != nil {
		t.Fatal(e)
	}
	if e = s.Complete(ctx, second, manifest); e != nil {
		t.Fatalf("ambiguous acknowledgement readback: %v", e)
	}
	if e = s.Fail(ctx, second, "after done"); !errors.Is(e, ErrFenced) {
		t.Fatalf("terminal changed: %v", e)
	}
	if e = s.Unavailable(ctx, second); e != nil {
		t.Fatal(e)
	}
	done, e := s.Get(ctx, job.ID, "conformer")
	if e != nil || done.State != "done" || done.Available {
		t.Fatalf("availability changed execution: %+v %v", done, e)
	}
	if _, e = s.Get(ctx, job.ID, "orbital"); !errors.Is(e, pgx.ErrNoRows) {
		t.Fatalf("cross kind read: %v", e)
	}
	replacement, _, e := s.Admit(ctx, "conformer", "key", []byte(`{}`), 20)
	if e != nil || replacement.ID == job.ID {
		t.Fatalf("replacement: %+v %v", replacement, e)
	}
}
func TestPendingDoesNotExpireAndRetryExhausts(t *testing.T) {
	s := testStore(t)
	ctx := context.Background()
	j, _, e := s.Admit(ctx, "orbital", "aged", []byte(`{}`), 20)
	if e != nil {
		t.Fatal(e)
	}
	_, e = s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET created_at=now()-interval '1 year' WHERE id=$1`, j.ID)
	if e != nil {
		t.Fatal(e)
	}
	if e = s.Reap(ctx, 3); e != nil {
		t.Fatal(e)
	}
	p, e := s.Get(ctx, j.ID, j.Kind)
	if e != nil || p.State != "pending" {
		t.Fatal(p, e)
	}
	for range 3 {
		attempt, e := s.Claim(ctx, "crashed", time.Minute)
		if e != nil {
			t.Fatal(e)
		}
		_, e = s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET lease_until=now()-interval '1 second' WHERE id=$1`, attempt.ID)
		if e != nil {
			t.Fatal(e)
		}
		if e = s.Reap(ctx, 3); e != nil {
			t.Fatal(e)
		}
	}
	p, e = s.Get(ctx, j.ID, j.Kind)
	if e != nil || p.State != "error" || p.Attempt != 3 {
		b, _ := json.Marshal(p)
		t.Fatalf("%s %v", b, e)
	}
}
func TestLegacyUniqueIndexRejectsRace(t *testing.T) {
	s := testStore(t)
	ctx := context.Background()
	_, e := s.Pool.Exec(ctx, `INSERT INTO public.jobs(id,cache_key,status,request_params) VALUES('a','key','pending','{}')`)
	if e != nil {
		t.Fatal(e)
	}
	_, e = s.Pool.Exec(ctx, `INSERT INTO public.jobs(id,cache_key,status,request_params) VALUES('b','key','running','{}')`)
	if e == nil {
		t.Fatal("duplicate active legacy key accepted")
	}
}

func TestInterruptedProcessRetriesWithNewAuthority(t *testing.T) {
	s := testStore(t)
	ctx := context.Background()
	j, _, e := s.Admit(ctx, "orbital", "interrupted", []byte(`{}`), 10)
	if e != nil {
		t.Fatal(e)
	}
	first, e := s.Claim(ctx, "old", time.Minute)
	if e != nil {
		t.Fatal(e)
	}
	if e = s.Retry(ctx, first, "supervisor stopping", 3); e != nil {
		t.Fatal(e)
	}
	if e = s.Complete(ctx, first, []byte(`{}`)); !errors.Is(e, ErrFenced) {
		t.Fatal("requeued token published", e)
	}
	next, e := s.Claim(ctx, "next", time.Minute)
	if e != nil || *first.Token == *next.Token || next.Error != nil {
		t.Fatal(next, e)
	}
	if next.ID != j.ID {
		t.Fatal("retry changed public job ID")
	}
}
