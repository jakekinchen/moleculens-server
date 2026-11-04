// Package store owns durable transitions. Engine subprocesses never import it.
package store

import (
	"context"
	"crypto/rand"
	"crypto/sha256"
	"embed"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"sort"
	"time"

	"github.com/jackc/pgx/v5"
	"github.com/jackc/pgx/v5/pgxpool"
)

//go:embed migrations/*.sql
var migrations embed.FS
var ErrFull = errors.New("queue capacity reached")
var ErrFenced = errors.New("attempt no longer owns a live lease")

type Store struct{ Pool *pgxpool.Pool }
type Job struct {
	ID        string          `json:"jobId"`
	Kind      string          `json:"kind"`
	Key       string          `json:"cacheKey"`
	Request   json.RawMessage `json:"-"`
	State     string          `json:"status"`
	Created   time.Time       `json:"createdAt"`
	Started   *time.Time      `json:"-"`
	Completed *time.Time      `json:"-"`
	Attempt   int             `json:"attempt"`
	Token     *string         `json:"-"`
	Owner     *string         `json:"-"`
	Lease     *time.Time      `json:"-"`
	Manifest  json.RawMessage `json:"-"`
	Available bool            `json:"-"`
	Error     *string         `json:"errorMessage"`
}

const columns = `id, kind, cache_key, request, state, created_at, started_at, completed_at, attempt, token, owner, lease_until, manifest, available, error_message`

func scan(row pgx.Row) (*Job, error) {
	j := new(Job)
	err := row.Scan(&j.ID, &j.Kind, &j.Key, &j.Request, &j.State, &j.Created, &j.Started, &j.Completed, &j.Attempt, &j.Token, &j.Owner, &j.Lease, &j.Manifest, &j.Available, &j.Error)
	return j, err
}
func ID() string {
	b := make([]byte, 16)
	if _, err := rand.Read(b); err != nil {
		panic(err)
	}
	return hex.EncodeToString(b)
}
func Open(ctx context.Context, dsn string) (*Store, error) {
	cfg, err := pgxpool.ParseConfig(dsn)
	if err != nil {
		return nil, err
	}
	cfg.MaxConns = 12
	cfg.ConnConfig.RuntimeParams["statement_timeout"] = "5000"
	p, err := pgxpool.NewWithConfig(ctx, cfg)
	if err != nil {
		return nil, err
	}
	if err = p.Ping(ctx); err != nil {
		p.Close()
		return nil, err
	}
	return &Store{p}, nil
}
func (s *Store) Migrate(ctx context.Context) error {
	var encoding string
	if err := s.Pool.QueryRow(ctx, "SHOW server_encoding").Scan(&encoding); err != nil {
		return err
	}
	if encoding != "UTF8" {
		return fmt.Errorf("UTF8 database required, found %s", encoding)
	}
	tx, err := s.Pool.Begin(ctx)
	if err != nil {
		return err
	}
	defer tx.Rollback(ctx)
	if _, err = tx.Exec(ctx, `SELECT pg_advisory_xact_lock(782134092)`); err != nil {
		return err
	}
	if _, err = tx.Exec(ctx, `CREATE TABLE IF NOT EXISTS public.moleculens_migrations(name text PRIMARY KEY, sha256 text NOT NULL, applied_at timestamptz NOT NULL DEFAULT now())`); err != nil {
		return err
	}
	entries, err := migrations.ReadDir("migrations")
	if err != nil {
		return err
	}
	sort.Slice(entries, func(i, j int) bool { return entries[i].Name() < entries[j].Name() })
	for _, e := range entries {
		data, err := migrations.ReadFile("migrations/" + e.Name())
		if err != nil {
			return err
		}
		sum := sha256.Sum256(data)
		digest := hex.EncodeToString(sum[:])
		var old string
		err = tx.QueryRow(ctx, `SELECT sha256 FROM public.moleculens_migrations WHERE name=$1`, e.Name()).Scan(&old)
		if err == nil {
			if old != digest {
				return fmt.Errorf("migration %s checksum changed", e.Name())
			}
			continue
		}
		if !errors.Is(err, pgx.ErrNoRows) {
			return err
		}
		if _, err = tx.Exec(ctx, string(data)); err != nil {
			return fmt.Errorf("migration %s: %w", e.Name(), err)
		}
		if _, err = tx.Exec(ctx, `INSERT INTO public.moleculens_migrations(name,sha256) VALUES($1,$2)`, e.Name(), digest); err != nil {
			return err
		}
	}
	return tx.Commit(ctx)
}
func (s *Store) Check(ctx context.Context) error {
	_, err := s.Pool.Exec(ctx, `SELECT token FROM moleculens_v2.jobs LIMIT 0`)
	return err
}
func (s *Store) Get(ctx context.Context, id, kind string) (*Job, error) {
	return scan(s.Pool.QueryRow(ctx, `SELECT `+columns+` FROM moleculens_v2.jobs WHERE id=$1 AND kind=$2`, id, kind))
}
func (s *Store) Cached(ctx context.Context, key, kind string) (*Job, error) {
	return scan(s.Pool.QueryRow(ctx, `SELECT `+columns+` FROM moleculens_v2.jobs WHERE cache_key=$1 AND kind=$2 AND state='done' AND available ORDER BY completed_at DESC LIMIT 1`, key, kind))
}

// Admit serializes the bounded count and active-key decision, not computation.
func (s *Store) Admit(ctx context.Context, kind, key string, request []byte, limit int) (*Job, bool, error) {
	tx, err := s.Pool.Begin(ctx)
	if err != nil {
		return nil, false, err
	}
	defer tx.Rollback(ctx)
	if _, err = tx.Exec(ctx, `SELECT pg_advisory_xact_lock(782134093)`); err != nil {
		return nil, false, err
	}
	j, err := scan(tx.QueryRow(ctx, `SELECT `+columns+` FROM moleculens_v2.jobs WHERE cache_key=$1 AND kind=$2 AND (state IN ('pending','running') OR (state='done' AND available)) ORDER BY (state='done') DESC,created_at DESC LIMIT 1`, key, kind))
	if err == nil {
		return j, j.State == "done", tx.Commit(ctx)
	}
	if !errors.Is(err, pgx.ErrNoRows) {
		return nil, false, err
	}
	var n int
	if err = tx.QueryRow(ctx, `SELECT count(*) FROM moleculens_v2.jobs WHERE state IN ('pending','running')`).Scan(&n); err != nil {
		return nil, false, err
	}
	if n >= limit {
		return nil, false, ErrFull
	}
	j, err = scan(tx.QueryRow(ctx, `INSERT INTO moleculens_v2.jobs(id,kind,cache_key,request) VALUES($1,$2,$3,$4) RETURNING `+columns, "g_"+ID(), kind, key, request))
	if err != nil {
		return nil, false, err
	}
	return j, false, tx.Commit(ctx)
}
func (s *Store) Claim(ctx context.Context, owner string, lease time.Duration) (*Job, error) {
	return s.ClaimCompatible(ctx, owner, lease, nil)
}
func (s *Store) ClaimCompatible(ctx context.Context, owner string, lease time.Duration, versions []byte) (*Job, error) {
	return scan(s.Pool.QueryRow(ctx, `UPDATE moleculens_v2.jobs SET state='running',started_at=now(),error_message=NULL,attempt=attempt+1,token=$1,owner=$2,lease_until=now()+$3*interval '1 second' WHERE id=(SELECT id FROM moleculens_v2.jobs WHERE state='pending' AND ($4::jsonb IS NULL OR request->>'compatibility'=$4::jsonb->>kind) ORDER BY created_at FOR UPDATE SKIP LOCKED LIMIT 1) RETURNING `+columns, ID(), owner, lease.Seconds(), versions))
}
func (s *Store) Renew(ctx context.Context, j *Job, lease time.Duration) error {
	tag, err := s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET lease_until=now()+$3*interval '1 second' WHERE id=$1 AND token=$2 AND state='running' AND lease_until>now()`, j.ID, j.Token, lease.Seconds())
	if err != nil {
		return err
	}
	if tag.RowsAffected() != 1 {
		return ErrFenced
	}
	return nil
}

// Complete is conditional. A lost commit acknowledgement is resolved by reading
// the same token and exact manifest back; ambiguous files are retained.
func (s *Store) Complete(ctx context.Context, j *Job, manifest []byte) error {
	tag, err := s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET state='done',completed_at=now(),lease_until=NULL,manifest=$3,available=true,error_message=NULL WHERE id=$1 AND token=$2 AND state='running' AND lease_until>now()`, j.ID, j.Token, manifest)
	if err == nil && tag.RowsAffected() == 1 {
		return nil
	}
	c, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()
	var same bool
	readErr := s.Pool.QueryRow(c, `SELECT state='done' AND token=$2 AND manifest=$3::jsonb FROM moleculens_v2.jobs WHERE id=$1`, j.ID, j.Token, manifest).Scan(&same)
	if readErr == nil && same {
		return nil
	}
	if err != nil {
		return fmt.Errorf("completion acknowledgement unresolved; retain attempt: %w", err)
	}
	return ErrFenced
}
func (s *Store) Fail(ctx context.Context, j *Job, message string) error {
	if len(message) > 2048 {
		message = message[:2048]
	}
	tag, err := s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET state='error',completed_at=now(),lease_until=NULL,error_message=$3 WHERE id=$1 AND token=$2 AND state='running' AND lease_until>now()`, j.ID, j.Token, message)
	if err != nil {
		return err
	}
	if tag.RowsAffected() != 1 {
		return ErrFenced
	}
	return nil
}
func (s *Store) Unavailable(ctx context.Context, j *Job) error {
	_, err := s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET available=false WHERE id=$1 AND state='done'`, j.ID)
	return err
}

// A lease grants publication authority, not proof that a process stopped. New
// attempts get private paths and tokens; stale processes cannot publish.
func (s *Store) Reap(ctx context.Context, maxAttempts int) error {
	_, err := s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET state=CASE WHEN attempt<$1 THEN 'pending' ELSE 'error' END,completed_at=CASE WHEN attempt<$1 THEN NULL ELSE now() END,lease_until=NULL,error_message='worker lease expired' WHERE state='running' AND lease_until<=now()`, maxAttempts)
	return err
}
func (s *Store) Heartbeat(ctx context.Context, owner string, caps []byte) error {
	_, err := s.Pool.Exec(ctx, `INSERT INTO moleculens_v2.workers(id,heartbeat,capabilities) VALUES($1,now(),$2) ON CONFLICT(id) DO UPDATE SET heartbeat=now(),capabilities=EXCLUDED.capabilities`, owner, caps)
	return err
}
func (s *Store) Ready(ctx context.Context, versions string) error {
	var n int
	err := s.Pool.QueryRow(ctx, `SELECT count(*) FROM jsonb_each_text($1::jsonb) expected WHERE NOT EXISTS (
 SELECT 1 FROM moleculens_v2.workers w WHERE heartbeat>now()-interval '30 seconds'
 AND capabilities->>'psi4' NOT IN ('','unavailable')
 AND (capabilities->>'compatibility')::jsonb->>expected.key=expected.value
)`, versions).Scan(&n)
	if err != nil {
		return err
	}
	if n != 0 {
		return errors.New("no recently qualified worker for every required operation")
	}
	return nil
}

// Retry is used only after the supervised process group is gone. Scientific
// rejection/deadline failures use Fail; infrastructure interruption is bounded.
func (s *Store) Retry(ctx context.Context, j *Job, reason string, maxAttempts int) error {
	if len(reason) > 2048 {
		reason = reason[:2048]
	}
	tag, e := s.Pool.Exec(ctx, `UPDATE moleculens_v2.jobs SET state=CASE WHEN attempt<$3 THEN 'pending' ELSE 'error' END,completed_at=CASE WHEN attempt<$3 THEN NULL ELSE now() END,lease_until=NULL,error_message=$4 WHERE id=$1 AND token=$2 AND state='running' AND lease_until>now()`, j.ID, j.Token, maxAttempts, reason)
	if e != nil {
		return e
	}
	if tag.RowsAffected() != 1 {
		return ErrFenced
	}
	return nil
}
