package store

import (
	"context"
	"encoding/json"
	"errors"
	"time"

	"github.com/jackc/pgx/v5"
)

type LegacyJob struct {
	ID, Key, State, Kind string
	Params, Meta         map[string]any
	Created              time.Time
	Compute              *float64
	Error                *string
}

func (s *Store) Legacy(ctx context.Context, id, kind string) (*LegacyJob, error) {
	j := new(LegacyJob)
	var params, meta []byte
	err := s.Pool.QueryRow(ctx, `SELECT id,cache_key,status,request_params,result_meta,created_at,compute_time_ms,error_message FROM public.jobs WHERE id=$1`, id).Scan(&j.ID, &j.Key, &j.State, &params, &meta, &j.Created, &j.Compute, &j.Error)
	if err != nil {
		return nil, err
	}
	if err = json.Unmarshal(params, &j.Params); err != nil {
		return nil, err
	}
	if len(meta) > 0 {
		if err = json.Unmarshal(meta, &j.Meta); err != nil {
			return nil, err
		}
	}
	j.Kind, _ = j.Params["job_type"].(string)
	if j.Kind == "" {
		j.Kind = "orbital"
		if j.Params["basis"] == "electrostatics" {
			j.Kind = "electrostatics"
		}
	}
	if j.Kind != kind {
		return nil, pgx.ErrNoRows
	}
	return j, nil
}
func (s *Store) LegacyArtifacts(ctx context.Context, key, kind string) ([]string, map[string]any, error) {
	var paths, meta []byte
	var id string
	err := s.Pool.QueryRow(ctx, `SELECT id FROM public.jobs WHERE cache_key=$1 AND status='done' ORDER BY created_at DESC LIMIT 1`, key).Scan(&id)
	if err != nil {
		return nil, nil, err
	}
	if _, err = s.Legacy(ctx, id, kind); err != nil {
		return nil, nil, err
	}
	err = s.Pool.QueryRow(ctx, `SELECT artifact_paths,result_meta FROM public.cache_entries WHERE cache_key=$1`, key).Scan(&paths, &meta)
	if err != nil {
		return nil, nil, err
	}
	var names []string
	var m map[string]any
	if json.Unmarshal(paths, &names) != nil || json.Unmarshal(meta, &m) != nil {
		return nil, nil, errors.New("invalid legacy metadata")
	}
	return names, m, nil
}
