-- Sole schema owner: the explicit `moleculens migrate` command.
-- Run after quiescing/draining legacy workers. Duplicate active legacy keys
-- deliberately fail the unique index instead of silently deleting user jobs.
CREATE TABLE IF NOT EXISTS public.jobs (
 id varchar(64) PRIMARY KEY, cache_key varchar(64) NOT NULL, status text NOT NULL DEFAULT 'pending',
 request_params json NOT NULL, result_meta json, error_message text,
 created_at timestamptz NOT NULL DEFAULT now(), started_at timestamptz, completed_at timestamptz,
 compute_time_ms double precision, worker_id varchar(64), attempt_count integer NOT NULL DEFAULT 0,
 attempt_token varchar(64)
);
ALTER TABLE public.jobs ADD COLUMN IF NOT EXISTS attempt_token varchar(64);
CREATE UNIQUE INDEX IF NOT EXISTS uq_jobs_active_key ON public.jobs(cache_key)
 WHERE status IN ('pending','queued','running');
CREATE TABLE IF NOT EXISTS public.cache_entries (
 cache_key varchar(64) PRIMARY KEY, params_hash varchar(64) NOT NULL, result_meta json NOT NULL,
 artifact_paths json NOT NULL, created_at timestamptz NOT NULL DEFAULT now(),
 hit_count integer NOT NULL DEFAULT 0, last_accessed_at timestamptz
);
CREATE SCHEMA IF NOT EXISTS moleculens_v2;
CREATE TABLE moleculens_v2.jobs (
 id text PRIMARY KEY, kind text NOT NULL CHECK (kind IN ('orbital','electrostatics','conformer','nomenclature')),
 cache_key text NOT NULL, request jsonb NOT NULL,
 state text NOT NULL DEFAULT 'pending' CHECK (state IN ('pending','running','done','error')),
 created_at timestamptz NOT NULL DEFAULT now(), started_at timestamptz, completed_at timestamptz,
 attempt integer NOT NULL DEFAULT 0, token text, owner text, lease_until timestamptz,
 manifest jsonb, available boolean NOT NULL DEFAULT false, error_message text,
 CHECK ((state = 'running') = (lease_until IS NOT NULL)),
 CHECK (state <> 'done' OR manifest IS NOT NULL)
);
CREATE UNIQUE INDEX active_key ON moleculens_v2.jobs(cache_key) WHERE state IN ('pending','running');
CREATE INDEX claim_order ON moleculens_v2.jobs(created_at) WHERE state = 'pending';
CREATE INDEX cached_results ON moleculens_v2.jobs(cache_key, completed_at DESC) WHERE state = 'done' AND available;
CREATE TABLE moleculens_v2.workers (id text PRIMARY KEY, heartbeat timestamptz NOT NULL, capabilities jsonb NOT NULL);
