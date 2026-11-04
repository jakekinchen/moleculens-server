package httpapi

import (
	"context"
	"encoding/json"
	"errors"
	"io"
	"log/slog"
	"mime"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/jackc/pgx/v5"
	"github.com/jakekinchen/moleculens-server/internal/contract"
	"github.com/jakekinchen/moleculens-server/internal/engine"
	"github.com/jakekinchen/moleculens-server/internal/store"
)

var resultReads = make(chan struct{}, 2)
var errReadCapacity = errors.New("result read capacity reached")

type API struct {
	Store            *store.Store
	Root, LegacyRoot string
	Capacity         int
	Origins          []string
	Quiesced         bool
}

func write(w http.ResponseWriter, status int, v any) {
	w.Header().Set("Content-Type", "application/json")
	w.Header().Set("Cache-Control", "no-store")
	w.WriteHeader(status)
	json.NewEncoder(w).Encode(v)
}
func failure(w http.ResponseWriter, status int, code, msg string) {
	write(w, status, map[string]any{"detail": map[string]string{"code": code, "message": msg}})
}
func resultFailure(w http.ResponseWriter, err error, unavailable string) {
	if errors.Is(err, errReadCapacity) {
		w.Header().Set("Retry-After", "1")
		failure(w, 503, "busy", "Result read capacity reached")
		return
	}
	failure(w, 503, "result_unavailable", unavailable)
}
func (a *API) Handler() http.Handler {
	mux := http.NewServeMux()
	mux.HandleFunc("GET /live", func(w http.ResponseWriter, r *http.Request) { write(w, 200, map[string]any{"status": "alive"}) })
	health := func(w http.ResponseWriter, r *http.Request) {
		c, cancel := context.WithTimeout(r.Context(), 3*time.Second)
		defer cancel()
		versions, _ := json.Marshal(contract.Versions)
		err := a.Store.Ready(c, string(versions))
		if err == nil {
			err = engine.CheckSpace(a.Root)
		}
		if err != nil {
			failure(w, 503, "not_ready", "No ready chemistry worker or database unavailable")
			return
		}
		write(w, 200, map[string]any{"status": "healthy", "version": "go-2", "database": true, "cacheDir": true})
	}
	mux.HandleFunc("GET /health", health)
	mux.HandleFunc("GET /ready", health)
	for path, kind := range map[string]string{"orbitals": "orbital", "electrostatics": "electrostatics", "conformers": "conformer"} {
		mux.HandleFunc("POST /v1/"+path+"/compute", func(w http.ResponseWriter, r *http.Request) { a.submit(w, r, kind) })
		mux.HandleFunc("GET /v1/"+path+"/jobs/{id}", func(w http.ResponseWriter, r *http.Request) { a.get(w, r, kind) })
		mux.HandleFunc("GET /v1/"+path+"/artifacts/{key}/{file}", func(w http.ResponseWriter, r *http.Request) { a.artifact(w, r, kind) })
	}
	mux.HandleFunc("POST /v1/nomenclature/explain", func(w http.ResponseWriter, r *http.Request) { a.submit(w, r, "nomenclature") })
	// This metric describes durable current jobs; it does not fabricate legacy
	// in-process cache-hit counters across service restarts.
	mux.HandleFunc("GET /metrics/cache", func(w http.ResponseWriter, r *http.Request) {
		var pending, running, done, failed int
		err := a.Store.Pool.QueryRow(r.Context(), `SELECT count(*) FILTER(WHERE state='pending'),count(*) FILTER(WHERE state='running'),count(*) FILTER(WHERE state='done'),count(*) FILTER(WHERE state='error') FROM moleculens_v2.jobs`).Scan(&pending, &running, &done, &failed)
		if err != nil {
			failure(w, 503, "database_unavailable", "Metrics unavailable")
			return
		}
		write(w, 200, map[string]any{"scope": "durable_jobs", "pending": pending, "running": running, "done": done, "error": failed})
	})
	slots := make(chan struct{}, 4)
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		origin := r.Header.Get("Origin")
		for _, o := range a.Origins {
			if o == origin || o == "*" {
				w.Header().Set("Access-Control-Allow-Origin", o)
				w.Header().Add("Vary", "Origin")
				w.Header().Set("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
				w.Header().Set("Access-Control-Allow-Headers", "Content-Type")
				break
			}
		}
		if r.Method == http.MethodOptions {
			w.WriteHeader(204)
			return
		}
		select {
		case slots <- struct{}{}:
			defer func() { <-slots }()
		default:
			failure(w, 429, "busy", "HTTP capacity reached")
			return
		}
		mux.ServeHTTP(w, r)
	})
}
func (a *API) submit(w http.ResponseWriter, r *http.Request, kind string) {
	if a.Quiesced {
		failure(w, 503, "quiesced", "Admission is stopped for maintenance")
		return
	}
	if err := engine.CheckSpace(a.Root); err != nil {
		failure(w, 503, "storage_unavailable", "Artifact storage unavailable")
		return
	}
	data, err := io.ReadAll(http.MaxBytesReader(w, r.Body, 4_000_000))
	if err != nil {
		failure(w, 413, "body_too_large", "Request exceeds limit")
		return
	}
	req, err := contract.Normalize(kind, data)
	if err != nil {
		status := 422
		if strings.HasPrefix(err.Error(), "SDF:") {
			status = 400
		}
		failure(w, status, "invalid_request", err.Error())
		return
	}
	body, err := json.Marshal(req)
	if err != nil {
		failure(w, 500, "encoding_error", "Could not encode request")
		return
	}
	var j *store.Job
	var cached bool
	// A corrupt cache pointer is invalidated independently, then one bounded
	// admission retry creates/reuses a fresh job. Historical done is unchanged.
	for attempt := 0; attempt < 2; attempt++ {
		j, cached, err = a.Store.Admit(r.Context(), kind, req.Key, body, a.Capacity)
		if err != nil {
			if errors.Is(err, store.ErrFull) {
				w.Header().Set("Retry-After", "2")
				failure(w, 429, "queue_full", err.Error())
			} else {
				failure(w, 503, "database_unavailable", "Admission unavailable")
			}
			return
		}
		if cached {
			if _, err = a.payload(j); err != nil {
				if errors.Is(err, errReadCapacity) {
					resultFailure(w, err, "Result unavailable")
					return
				}
				if e := a.Store.Unavailable(r.Context(), j); e != nil {
					failure(w, 503, "result_unavailable", "Result unavailable")
					return
				}
				continue
			}
		}
		break
	}
	wait := time.Duration(req.WaitMS) * time.Millisecond
	if kind == "nomenclature" {
		wait = 9 * time.Second
	}
	deadline := time.NewTimer(wait)
	defer deadline.Stop()
	ticker := time.NewTicker(50 * time.Millisecond)
	defer ticker.Stop()
waitLoop:
	for wait > 0 && j.State != "done" && j.State != "error" {
		select {
		case <-r.Context().Done():
			return
		case <-deadline.C:
			break waitLoop
		case <-ticker.C:
			j, err = a.Store.Get(r.Context(), j.ID, kind)
			if err != nil {
				failure(w, 503, "database_unavailable", "Polling unavailable")
				return
			}
		}
	}
	if kind == "nomenclature" {
		if j.State == "done" {
			p, e := a.payload(j)
			if e != nil {
				resultFailure(w, e, "Naming result unavailable")
				return
			}
			write(w, 200, json.RawMessage(p))
			return
		}
		if j.State == "error" {
			status, code, msg := 422, "naming_failed", message(j.Error)
			if strings.HasPrefix(msg, "adapter: ") {
				var detail struct {
					Code       string `json:"code"`
					Message    string `json:"message"`
					HTTPStatus int    `json:"httpStatus"`
				}
				if json.Unmarshal([]byte(strings.TrimPrefix(msg, "adapter: ")), &detail) == nil {
					if detail.HTTPStatus == 500 {
						status = 500
					}
					code = detail.Code
					msg = detail.Message
				}
			}
			failure(w, status, code, msg)
			return
		}
		failure(w, 503, "naming_busy", "Naming job is still pending; retry the same request")
		return
	}
	a.respond(w, r, j, cached)
}
func message(p *string) string {
	if p == nil {
		return ""
	}
	return *p
}
func (a *API) payload(j *store.Job) ([]byte, error) {
	select {
	case resultReads <- struct{}{}:
		defer func() { <-resultReads }()
	default:
		return nil, errReadCapacity
	}
	var m engine.Manifest
	if err := json.Unmarshal(j.Manifest, &m); err != nil {
		return nil, err
	}
	if j.Token == nil || m.Key != j.Key || m.Token != *j.Token {
		return nil, errors.New("manifest pointer mismatch")
	}
	return engine.Load(filepath.Join(a.Root, j.Key, *j.Token), &m)
}
func (a *API) respond(w http.ResponseWriter, r *http.Request, j *store.Job, cached bool) {
	status := j.State
	if j.Kind == "conformer" {
		if status == "error" {
			status = "failed"
		} else if status != "done" {
			status = "pending"
		}
	}
	var publicError *string
	if j.State == "error" {
		publicError = j.Error
	}
	response := map[string]any{"jobId": j.ID, "cacheKey": j.Key, "status": status, "cached": cached, "createdAt": j.Created, "errorMessage": publicError}
	if j.State == "done" {
		p, e := a.payload(j)
		if e != nil {
			if !errors.Is(e, errReadCapacity) {
				a.Store.Unavailable(r.Context(), j)
			}
			resultFailure(w, e, "Job completed but its artifacts are unavailable")
			return
		}
		var m engine.Manifest
		json.Unmarshal(j.Manifest, &m)
		response["computeTimeMs"] = m.ElapsedMS
		if j.Kind == "conformer" {
			var fields map[string]any
			if json.Unmarshal(p, &fields) != nil {
				failure(w, 503, "result_unavailable", "Invalid result")
				return
			}
			for k, v := range fields {
				response[k] = v
			}
		} else {
			response["result"] = json.RawMessage(p)
		}
	}
	write(w, 200, response)
}
func (a *API) get(w http.ResponseWriter, r *http.Request, kind string) {
	id := r.PathValue("id")
	if !strings.HasPrefix(id, "g_") {
		a.legacy(w, r, id, kind)
		return
	}
	j, err := a.Store.Get(r.Context(), id, kind)
	if err != nil {
		if errors.Is(err, pgx.ErrNoRows) {
			failure(w, 404, "not_found", "Job not found")
		} else {
			failure(w, 503, "database_unavailable", "Job unavailable")
		}
		return
	}
	a.respond(w, r, j, false)
}
func safeSegment(s string) bool {
	return s != "" && s != "." && s != ".." && filepath.Base(s) == s && !strings.ContainsAny(s, "/\\")
}
func (a *API) artifact(w http.ResponseWriter, r *http.Request, kind string) {
	key, name := r.PathValue("key"), r.PathValue("file")
	if !safeSegment(key) || !safeSegment(name) {
		failure(w, 400, "invalid_path", "Invalid artifact path")
		return
	}
	j, e := a.Store.Cached(r.Context(), key, kind)
	if e == nil {
		var m engine.Manifest
		if json.Unmarshal(j.Manifest, &m) == nil {
			if _, e = a.payload(j); e == nil {
				for _, f := range m.Files {
					if f.Name == name {
						serveArtifact(w, r, filepath.Join(a.Root, key, m.Token, name))
						return
					}
				}
			}
		}
		resultFailure(w, e, "Artifact unavailable")
		return
	}
	if !errors.Is(e, pgx.ErrNoRows) {
		failure(w, 503, "database_unavailable", "Artifact lookup unavailable")
		return
	}
	names, meta, e := a.Store.LegacyArtifacts(r.Context(), key, kind)
	if e == nil {
		dir, _ := meta["_artifact_dir"].(string)
		if dir != "" && !safeSegment(dir) {
			failure(w, 503, "result_unavailable", "Invalid artifact pointer")
			return
		}
		for _, f := range names {
			if f == name {
				path := filepath.Join(a.LegacyRoot, key, dir, name)
				info, err := os.Lstat(path)
				if err == nil && info.Mode().IsRegular() {
					serveArtifact(w, r, path)
					return
				}
			}
		}
	}
	failure(w, 404, "not_found", "Artifact not found")
}
func serveArtifact(w http.ResponseWriter, r *http.Request, path string) {
	contentType := "application/octet-stream"
	switch filepath.Ext(path) {
	case ".json":
		contentType = "application/json"
	case ".cube":
		contentType = "chemical/x-gaussian-cube"
	}
	w.Header().Set("Content-Type", contentType)
	w.Header().Set("Content-Disposition", mime.FormatMediaType("attachment", map[string]string{"filename": filepath.Base(path)}))
	http.ServeFile(w, r, path)
}
func rename(m map[string]any, keys map[string]string) map[string]any {
	out := map[string]any{}
	for k, v := range m {
		if n, ok := keys[k]; ok {
			k = n
		}
		out[k] = v
	}
	return out
}
func (a *API) legacy(w http.ResponseWriter, r *http.Request, id, kind string) {
	j, e := a.Store.Legacy(r.Context(), id, kind)
	if e != nil {
		if errors.Is(e, pgx.ErrNoRows) {
			failure(w, 404, "not_found", "Job not found")
		} else {
			failure(w, 503, "database_unavailable", "Legacy job unavailable")
		}
		return
	}
	status := j.State
	if kind == "conformer" {
		if status == "error" {
			status = "failed"
		} else if status != "done" {
			status = "pending"
		}
	}
	response := map[string]any{"jobId": j.ID, "cacheKey": j.Key, "status": status, "cached": false, "createdAt": j.Created, "computeTimeMs": j.Compute, "errorMessage": j.Error}
	if j.State == "done" {
		p, e := a.legacyPayload(j)
		if e != nil {
			slog.Warn("legacy artifacts unavailable", "job", j.ID, "error", e)
			failure(w, 503, "result_unavailable", "Historical job completed but artifacts are unavailable")
			return
		}
		if kind == "conformer" {
			for k, v := range p {
				response[k] = v
			}
		} else {
			response["result"] = p
		}
	}
	write(w, 200, response)
}
func (a *API) legacyPayload(j *store.LegacyJob) (map[string]any, error) {
	dir, _ := j.Meta["_artifact_dir"].(string)
	if !safeSegment(j.Key) || (dir != "" && !safeSegment(dir)) {
		return nil, errors.New("invalid legacy pointer")
	}
	root := filepath.Join(a.LegacyRoot, j.Key, dir)
	filename := "result.json"
	if j.Kind == "electrostatics" {
		filename = "electrostatics_result.json"
	}
	if j.Kind == "conformer" {
		filename = "structure.sdf"
	}
	info, e := os.Lstat(filepath.Join(root, filename))
	if e != nil {
		return nil, e
	}
	if !info.Mode().IsRegular() || info.Size() > engine.MaxResult {
		return nil, errors.New("invalid legacy file")
	}
	b, e := os.ReadFile(filepath.Join(root, filename))
	if e != nil {
		return nil, e
	}
	if j.Kind == "conformer" {
		if e = contract.ValidateSDF(string(b)); e != nil {
			return nil, e
		}
		if _, ok := j.Meta["quality"]; !ok {
			return nil, errors.New("legacy conformer metadata absent")
		}
		return map[string]any{"sdf3d": string(b), "timingMs": j.Meta["timings_ms"], "meta": rename(j.Meta, map[string]string{"max_attempts": "maxAttempts", "max_opt_iters": "maxOptIters", "add_hs": "addHs", "geom_version": "geomVersion", "has_metal": "hasMetal", "embed_status": "embedStatus", "opt_status": "optStatus"})}, nil
	}
	var p map[string]any
	if e = json.Unmarshal(b, &p); e != nil {
		return nil, e
	}
	if len(p) == 0 || p["meta"] == nil {
		return nil, errors.New("empty legacy result")
	}
	if j.Kind == "electrostatics" {
		d, ok := p["dipole"].(map[string]any)
		if !ok || p["surface"] == nil {
			return nil, errors.New("missing electrostatics fields")
		}
		p["dipole"] = rename(d, map[string]string{"origin_angstrom": "originAngstrom", "vector_debye": "vectorDebye", "magnitude_debye": "magnitudeDebye"})
		m, ok := p["meta"].(map[string]any)
		if !ok {
			return nil, errors.New("missing ESP metadata")
		}
		p["meta"] = rename(m, map[string]string{"charge_model": "chargeModel", "timings_ms": "timingsMs", "sdf_hash": "sdfHash", "molecular_charge": "molecularCharge"})
	} else {
		if p["orbitals"] == nil {
			return nil, errors.New("missing orbital data")
		}
	}
	return p, nil
}
