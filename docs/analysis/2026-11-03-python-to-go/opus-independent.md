# Moleculens server: independent Python-to-Go migration assessment (Opus)

Date: 2026-11-03. This is analysis only. No implementation, commit, push, deploy, live-service or live-data access happened.

## 0. Identity, scope, and how to read this

| Item | Verified value |
| --- | --- |
| Server repo | `/Users/kelly/Developer/moleculens-server`, branch `main`, HEAD `c9fce3c5d9f69212cc65d39e06ae3c18ef74eafb` |
| Remote | `git ls-remote origin refs/heads/main` returned the same hash, `c9fce3c…`. The dirty edits below are local only and are not on the remote. |
| Dirty tracked files | `README.md`, `docker/Dockerfile.api`, `pyproject.toml`, `src/moleculens/api/main.py`, `src/moleculens/api/schemas.py`, `uv.lock` (6 files, +156/−7) |
| Untracked files | `benchmarks/`, `coverage.xml`, `scripts/build_nomenclature_corpus.py`, `src/moleculens/api/nomenclature_routes.py`, `src/moleculens/compute/nomenclature.py`, `src/moleculens/compute/nomenclature_benchmark.py`, `tests/unit/test_nomenclature.py`, `var/dev.sqlite3` |
| Appeared during my session | `docs/analysis/2026-11-03-python-to-go/` was created at 11:12, after this session started. I did not create it and did not read its contents, to keep this assessment independent. |
| Consumer repo | `/Users/kelly/.openclaw/workspace/projects/moleculens`, HEAD `35dd5f4b4bb174b2005074b8e1e0f16a2702dd5b`. The four named files (`src/app/api/{orbitals,conformers,electrostatics}/schema.ts`, `src/shared/lib/overlays/cacheKey.ts`) are clean. `src/features/nomenclature/data/nomenclature.server.ts` is untracked. |
| Post-analysis state | Both repos are unchanged by me. Status and diffstat were re-checked at the end of the session. |

**Scope covered:** API (all four routers), queue and models, worker, caching and artifact publication, all compute modules, Docker, compose and Caddy, CI and deploy, tests, the uncommitted nomenclature work, and the consumer request builders named in the follow-up request.

**Method:** I read the code under the six code-lenses, did a state-machine extraction, ran disposable probes (SQLite or no DB, temp dirs in my scratchpad, run from a cwd with no `.env`), and checked primary sources for Psi4 and xtb-python behavior.

**Not read:** `.env`, credentials, `var/` contents, `docs/analysis/…`.

**Evidence labels used below:**
- **[RUN]**: I executed a check and observed the result.
- **[SRC]**: complete source trace, not executed.
- **[PRIMARY]**: confirmed against upstream library source or docs.
- **[KNOW]**: background knowledge, not verified in this session.

---

## 1. Short answer

1. **The defects that matter here do not depend on the language.** The highest-impact problems are an unguarded job lifecycle, non-atomic artifact publication, a cache identity taken on trust from an unauthenticated client, contract drift with the actual consumer, and several chemistry-semantics bugs. A Go rewrite fixes none of these automatically. A straight port would carry them over, and it adds a parity risk for the cache keys.
2. **What can move to Go (control plane):**
   - HTTP edge: strict request decoding, readiness that returns correct status codes, metrics.
   - Job submission and dedupe, job status reads, artifact serving.
   - Queue claim, lease, heartbeat and reaper.
   - An engine *supervisor* that runs chemistry as subprocesses with hard timeouts.
3. **What stays behind a chemistry-engine boundary (Python):**
   - Psi4 and cubeprop, cube parsing, marching cubes and mesh encoding.
   - xTB/Gasteiger charges and the ESP surface.
   - RDKit conformers.
   - Openclatura, RDKit and OPSIN nomenclature.
   - The worker-side SDF/chemistry parsing.

   There is no native Go Psi4 or xTB. I found no existing Go RDKit binding: RDKit exposes a C FFI ("MinimalLib"), wrapped today for Rust, not Go [PRIMARY search]. Porting numeric mesh code would be parity work with no measured benefit.
4. **One gray zone: the API-side cache identity** (`parse_sdf` plus `geometry_hash`). Either port it to Go behind a byte-for-byte parity corpus, or move identity computation to the engine side. Do not port it casually. It is a cache contract, and it has Python-specific float, negative-zero and sort semantics.
5. **Recommended order (detailed in §8):**
   - **0.** Evidence baseline.
   - **1.** Two-way contract agreement with the consumer.
   - **2.** Make Postgres the authoritative, guarded state machine (still Python).
   - **3.** Fix engine semantics, with cache-key versioning.
   - **Decision point:** Go or not, based on measured data.
   - **4–8.** Go strangler: read paths, then conformer submits, then the other submits once identity parity is proven, then the supervisor. Nomenclature stays a Python sidecar. FastAPI is retired last.

   Stopping after stage 3 is a legitimate outcome if stage 0 shows no API-tier problem.

---

## 2. Checks actually run (separate from source-only findings)

All probes live under `/private/tmp/claude-501/-Users-kelly-Developer-moleculens-server/c4a0daf6-1b5c-4ef9-8c7b-20967be17449/scratchpad/probes/`.

| # | Check | Result |
| --- | --- | --- |
| R0 | `git rev-parse`, `git status`, `git ls-remote` on both repos | As in §0 |
| R1 | Existing unit tests (`pytest tests/unit`, isolated cwd, `-p no:cacheprovider`) | **46 passed**, with 1 warning: Java is missing for py2opsin |
| R2 | `probe_queue.py`, SQLite queue lifecycle | P1 zombie completion reproduced. P2 terminal-state mutation reproduced. P3 duplicate active jobs reproduced by replaying submit's check-then-insert across two transactions. P5 stale predicate is true for an old PENDING job. |
| R3 | `probe_fix.py` | P4: a CacheEntry primary-key race leaves the second successful job stuck RUNNING, and the worker would then mark it failed. C8: RDKit clears the 2D conformer on embed failure, so the repo's "planar fallback" puts every atom at (0,0,0). C6b: Gasteiger respects M CHG when present. |
| R4 | `probe_chem.py` | C1: basis `6-31G*` becomes `6-31gstar`. C2: the parser rejects RDKit's own molblock with a blank title. C3: the parser rejects ≥100 atoms with ≥100 bonds. C4: Python's negative-zero normalization. C5: Pt, Au, Hg and Gd become carbon. C6: an `xtb-gfn2` request yields `rdkit-gasteiger`, and the request charge is ignored. C7: surface spacing distortion. C9: key fragmentation. |
| R5 | `probe_loop.py` (ASGI in-process) | Sync I/O in `async def` stalls unrelated requests. `/metrics/cache` completed at +51 ms alone and **+506 ms** next to a 0.5 s blocking poll. |
| R6 | `probe_compat.py` | Consumer-shaped bodies (transcribed from consumer source, not executed TS) validated against the server's Pydantic models. Lists of dropped fields, 422 cases and wire key casing are in §5. |
| R7 | `moleculens-nomenclature-benchmark` on `benchmarks/nomenclature_corpus.json`, no OPSIN | 398 cases in 7.13 s (about 18 ms/case). supportedNameRate 1.0, guardrailRecall 1.0. **opsinMatchRate not measured** (no Java), so the release gate cannot pass locally. |
| R8 | Primary sources | Psi4 `libmintsbasisset.make_filename`: lowercases, `*`→`s`, `+`→`p`, appends `.gbs`. Psi4 `export_mints.cc` / `vector.h`: `Vector.get(m)` reads flat irrep-blocked storage, and `dim(h=0)` returns the **irrep-0** dimension. Psi4 cubeprop manual: `CUBEPROP_ORBITALS` indices are in energy order ignoring symmetry. xtb-python `Calculator(param, numbers, positions, charge=None, uhf=None, …)`. RDKit CFFI blog plus search: no Go binding found. |

**Not run** (blocked locally or out of scope): Psi4, xTB, Java/OPSIN, Docker, Postgres, the Go toolchain (not installed), the consumer TypeScript code, anything touching production.

---

## 3. Ranked confirmed defects

Ranking is impact × reach × evidence strength. "Migration relevance" says whether a Go port would inherit the defect, block on it, or be a natural place to fix it.

### D1. The job lifecycle has no transition guards; terminal states are mutable (High) [RUN + SRC]
- **Rule:** only the current owner of a RUNNING job may complete or fail it, and terminal states are absorbing.
- **Where:**
  - `complete_job` sets DONE without checking status or worker: `src/moleculens/db/queue.py:608-617`.
  - `fail_job` sets ERROR unconditionally: `src/moleculens/db/queue.py:647-655`.
  - Stale auto-clear sets ERROR on a job a live worker still owns: `src/moleculens/db/queue.py:202-217` and `:238-241`.
  - Invalidation rewrites DONE to ERROR: `src/moleculens/db/queue.py:276-289`.
- **Counterexample (P1, reproduced):**
  1. Worker w1 runs job j1 for more than 900 s.
  2. A resubmit auto-clears j1 to ERROR and creates j2 (PENDING).
  3. w1 finishes, and `complete_job(j1)` flips **ERROR→DONE** and writes the cache entry.
  4. j2 stays PENDING and is later claimed. It recomputes into the same directory and overwrites the cache entry again.
- **Second counterexample (P2):** DONE followed by `fail_job` gives ERROR while a valid CacheEntry remains. The next submit recomputes, but `hit_count` was still incremented.
- **Migration relevance:** Go must not port this. The guards belong in SQL (§7), so both languages obey them during a strangler period.

### D2. Duplicate active jobs per cache key, and a CacheEntry race marks successful work as failed (High) [RUN replay + SRC]
- **Where:**
  - Check-then-insert with no lock or unique constraint: `src/moleculens/db/queue.py:419-464` and `:481-516`.
  - Indexes are non-unique: `src/moleculens/db/models.py:69-72`.
  - CacheEntry is inserted after a non-locking read: `src/moleculens/db/queue.py:620-628`.
- **Counterexample:**
  1. Two concurrent submits for the same key both see no active job. Two PENDING rows result (P3).
  2. When both complete, the second CacheEntry insert raises IntegrityError (P4).
  3. The worker's generic handler calls `fail_job` on a successful compute: `src/moleculens/worker/main.py:167-169`.
- Under Postgres READ COMMITTED the same interleaving applies: there is no row to lock and no constraint. The reproduction used SQLite with explicit interleaving.
- **Migration relevance:** this needs a partial unique index plus `ON CONFLICT`. It is language-independent.

### D3. Artifact publication is not atomic and shares one final directory (High) [SRC]
- **Where:**
  - The worker computes straight into the final `cache/<cache_key>/`: `src/moleculens/worker/main.py:88-89`.
  - Psi4 writes its output, cubes and `result.json` in place: `src/moleculens/compute/psi4_runner.py:95`, `:176`, `:292-293`.
  - Electrostatics writes in place: `src/moleculens/compute/electrostatics.py:249-251`.
  - Conformers write in place: `src/moleculens/compute/conformers.py:193-194`.
  - Readers check only that the file exists: `src/moleculens/db/queue.py:261-268`, `src/moleculens/api/routes.py:240-252`.
- **Counterexample:** in the D1/D2 duplicate scenario, j2 rewrites `result.json` while the API serves j1. A read during that write fails `json.load`, which triggers invalidation. Invalidation deletes the cache entry and turns j1 into ERROR (`routes.py:32-45`), and the cascade enqueues yet another job.
- The torn-read cascade itself is plausible (R5 in §4). The non-atomic write, and the absence of any staging/rename/manifest step, are confirmed.
- **Migration relevance:** an artifact publisher (staging dir, checksums, manifest, rename, then DB reference) is a natural control-plane responsibility. It is a good Go candidate, but only together with D1/D2.

### D4. No lease, heartbeat, timeout or retry; queue backlog becomes failures (High) [RUN predicate + SRC]
- **Unused settings:** `worker_job_timeout_seconds` and `max_retries` are defined but never used (`src/moleculens/core/config.py:35`, `:42`; confirmed by grep). `attempt_count` is incremented and never read.
- **Crashed worker:** SIGTERM only sets a flag (`src/moleculens/worker/main.py:31-35`), and Psi4 runs inline (`:133-142`). A crash, OOM or deploy recreation leaves the job RUNNING forever. It is only cleared when *another submit for the same key* arrives after 900 s (`queue.py:219-241`).
- **Stale pending jobs:** a job PENDING for longer than `job_stale_seconds` (900 s) is auto-errored *at claim time*, even though it never ran (`queue.py:194`, `:572-575`; predicate confirmed in P5). With one worker and Psi4 jobs taking minutes, a modest backlog turns queued work into "Stale pending job auto-cleared" errors.
- **Dead state:** `QUEUED` is never written (`models.py:21`).
- **Migration relevance:** the lease/reaper and engine supervisor are the strongest technical reasons for a separate control-plane component, whether Go or Python.

### D5. The unauthenticated `client_cache_key` replaces server-derived identity, enabling cache poisoning (High) [SRC + existing test]
- **Where:**
  - `src/moleculens/api/routes.py:100` and `src/moleculens/db/queue.py:395-397` (`molecule_identity = cache_identity or geometry_hash`).
  - `src/moleculens/api/electrostatics_routes.py:98-102`.
  - `/v1/*` is publicly proxied: `docker/Caddyfile.single:22-24`.
  - The existing test `tests/unit/test_job_queue.py:217-242` *asserts* that two different geometries sharing a client key reuse the same job.
- **Counterexample:**
  1. The consumer key is deterministic: SHA-256 of `overlay-v3`, the SDF text and the parameters (`cacheKey.ts:54-85`). Anyone holding a public molecule's SDF can compute it.
  2. An attacker POSTs a different SDF with that `clientCacheKey` and the same method/basis/grid/isovalue/orbitals.
  3. The server's `compute_cache_key` matches, so every consumer request for the real molecule receives the attacker's mesh.
- **Fallback amplifier:** if `crypto.subtle` is unavailable (for example an insecure context), `fallbackHash` produces a 32-bit key (`cacheKey.ts:7-14`, `:37-47`). Collisions then become plausible without any attacker (R4 in §4).
- **Migration relevance:** this is a trust-model decision (U3 in §10), not a language one. It must be settled before porting identity.

### D6. `6-31G*`-style basis names are rewritten to a name Psi4 cannot resolve (High for any non-STO basis) [RUN + PRIMARY]
- **Where:** `src/moleculens/api/schemas.py:65-69` maps `*` to `star`. The value is passed verbatim to `psi4.energy(..., basis=basis)` at `src/moleculens/compute/psi4_runner.py:119-122`. Psi4 `make_filename` maps `*` to `s`, giving `6-31gs.gbs` for `6-31g*`, but `6-31gstar` becomes `6-31gstar.gbs` [PRIMARY].
- That Psi4 ships no `6-31gstar.gbs` is [KNOW]. Not executed: Psi4 is not installed.
- **Counterexample:** `{"basis":"6-31G*"}` → 200 at submit, then the job fails in the worker.
- The golden test encodes the `star` naming (`tests/integration/test_golden.py:162`) and only exercises `sto-3g`.
- **Fix shape:** keep normalization in the cache key only, and pass Psi4 the canonical name.

### D7. The SDF parser rejects valid V2000 molfiles (High reach) [RUN]
- **Blank title line:** `content.strip()` drops the empty title line and shifts the header (`src/moleculens/compute/sdf_parser.py:91`). RDKit's default `MolToMolBlock` output, which has an empty first line, is rejected with "Invalid counts line" (C2). RDKit itself reads the same block.
- **≥100 atoms and ≥100 bonds:** the fixed-width counts field `aaabbb` run together (e.g. `104103`) is parsed as one integer (`sdf_parser.py:115-117`), so these molecules are rejected (C3).
- **Where it bites:**
  - Orbital and electrostatics submits without a client key fail with 400 (`queue.py:388-393`, `electrostatics_routes.py:100-101`).
  - With a client key, the parse happens in the worker and the job ends in ERROR.
  - The server's own conformer output inherits the input's title line, so a blank-title 2D input produces 3D output that the other endpoints reject.
- **Migration relevance:** a Go port must decide whether to replicate or fix this. Fixing changes behavior only for inputs that currently fail. Keys for accepted inputs must stay bit-identical (see S4 in §6.5).

### D8. Consumer/server contract drift is silently ignored (High) [RUN against server models + SRC consumer]
Pydantic's default `extra="ignore"` drops fields the consumer believes it controls. The full matrix is in §5. Highlights:
- **Orbitals:** `charge`, `multiplicity` and `density_isovalue` are dropped.
  - Consumer: `orbitals/schema.ts:18-25`, `:36-61`; `presentations/cache/route.ts:164-181`.
  - Server: `schemas.py:8-54`.
  - Result: ions and radicals are computed as neutral singlets, and the density mesh uses `isovalue` (0.02), not `densityIsovalue` (0.01). The consumer's key includes charge, multiplicity and densityIsovalue, so it believes these are honored.
- **Conformers:** the consumer never sends `params`. It sends `optimize`, `forceField` (default `mmff`), `seed` and `numConformers`, and all of them are dropped (`conformers/schema.ts:25-38`; `conformerClient.ts:40-46`).
  - The server always runs ETKDGv3 plus **UFF**, even when the consumer asked for `optimize:false` or MMFF.
  - The consumer's `method` is mapped to `geomVersion`, so it only changes the cache key, not the method.
- **Electrostatics:** `surface` is optional in the consumer (`electrostatics/schema.ts:20-26`). When it is omitted, the server uses probe 1.4 Å / padding 3.0 Å, while the consumer key assumes 0.9 / 1.2 (`cacheKey.ts:111-113`). Consumer-valid inputs also hit server 422s: method other than `xtb|xtb-gfn2|gasteiger`, padding below 1.0, |charge| above 10.
- **Response casing:** the consumer reads `meta.grid_spacing_angstrom` (`orbitals/schema.ts:115`), but the server emits `gridSpacingAngstrom` (FastAPI serializes by alias). The consumer's value is therefore always undefined.
- **Migration relevance:** this is the best argument for strict decoding at the edge, in Go (`DisallowUnknownFields`) or Pydantic (`extra="forbid"`). It needs a two-way decision first (stage 1).

### D9. No charge or multiplicity path to Psi4; charged species are computed as neutral (Medium-High) [SRC + KNOW]
- **Where:** `parse_sdf` ignores the atom-block charge field and `M  CHG` lines (`sdf_parser.py:129-139`, `:142-178`). `to_psi4_geometry` emits no charge/multiplicity line (`sdf_parser.py:31-47`), and the request schema has no charge (`schemas.py:8-54`).
- Psi4 defaults to charge 0 and the lowest consistent multiplicity [KNOW].
- **Counterexamples:**
  - Sulfate SDF with two `M CHG -1` lines: computed as neutral SO₄, wrong orbitals, cached as success.
  - NH₄⁺: computed as an odd-electron neutral species, so the RHF job fails.

### D10. Electrostatics degrades silently and caches the degraded result under the non-degraded key (Medium) [RUN + PRIMARY]
- **Silent fallback:** if xTB is unavailable or fails, the code falls back to Gasteiger (`charges.py:52-63`, `electrostatics.py:120-131`). The result is cached under the `xtb-gfn2` key (`electrostatics_routes.py:104-111`). In C6, `meta.method=xtb-gfn2` came with `charge_model=rdkit-gasteiger`.
  - The worker image deliberately skips xTB (`docker/Dockerfile.worker:29-31`; commit `5b33b80`), so in production every xTB request is Gasteiger. That is repo evidence; the image contents are unverified.
  - If xTB is later added, the old entries are still served.
- **Request charge ignored by Gasteiger:** `calculate_charges_gasteiger` takes no charge argument (`charges.py:190-270`). `_compute_simple_charges` forces neutrality (`charges.py:313`). In C6 the same pose without `M CHG` and with request charge +1 gave a total charge of −0.34. With `M CHG` present (C6b) it gave +1.00.
- **Same key, different inputs:** `geometry_hash` ignores `M CHG` and bonds, which Gasteiger consumes (`sdf_parser.py:49-70`, `charges.py:216`). In C6 the two blocks had the same key but different engine inputs, so whichever runs first is served for both.
- **Multiplicity ignored by xtb-python:** `Calculator(…, charge=charge)` passes no `uhf` (`charges.py:84`), confirmed against the [PRIMARY] signature. Multiplicity is still in the key and echoed in meta.
- **Unknown elements become carbon:** atomic number 6 (`charges.py:435`) and radius 1.70 Å (`surface_mesh.py:41`, `:57`). C5 showed this for Pt, Au, Hg and Gd.

### D11. Conformer embed failure caches a degenerate geometry (Medium) [RUN]
- **Where:** RDKit `EmbedMolecule` clears existing conformers. On failure, `_ensure_planar_coords` sees zero conformers and writes (0,0,0) for every atom (`conformers.py:65-71`, `:164-168`).
- **Result (C8):** status DONE, `quality=planar_fallback`, 1 distinct position out of 3 atoms. The label says "planar", but the geometry is a single point, and it is cached permanently.
- The failure was forced with an impossible coordinate map. Behavior after any real embed failure is identical.

### D12. Partial engine results are recorded as success (Medium) [SRC]
- **Where:** a missing cube file is skipped with a warning (`psi4_runner.py:216-222`). Per-orbital mesh exceptions are swallowed (`:258-259`). A missing density, or a density failure, is silently omitted (`:262-280`). The job is then DONE and cached.
- **Test gap:** the golden test tolerates all of this. Energies are optional (`test_golden.py:139-151`) and density is optional (`:153-157`).

### D13. Orbital energies are read in irrep-blocked order while cubes are chosen in energy order (Medium) [SRC + PRIMARY]
- **Where:** `epsilon = wfn.epsilon_a()` and `epsilon.get(index-1)`, bounded by `epsilon.dim()` (`psi4_runner.py:135-146`, `:162`). [PRIMARY] establishes that `get(m)` is a flat Pitzer/irrep-blocked offset and `dim()` is the irrep-0 size. Cubeprop indices are energy-ordered [PRIMARY].
- **Effect:** for an input frame Psi4 recognizes as symmetric, the code reports `None` or the wrong orbital's energy, and may skip `lumo+1`.
- **Likely case:** the golden water (O at the origin, H at ±x in the xz-plane) is C2v in its input frame. HOMO index 5 exceeds A1's dimension of 4, so `energyEv` would be omitted, and the test would still pass (D12). This manifestation is not executed (R1 in §4).

### D14. Blocking I/O inside `async def` handlers stalls the event loop (Medium) [RUN]
- **Where:** every DB and file call in `routes.py`, `electrostatics_routes.py` and `conformers_routes.py` is synchronous inside `async def`. `_wait_for_job` polls the DB every 100 ms for up to `waitMs`, and the consumer default is 2000 ms (`conformers_routes.py:170-184`).
- **Probe (R5):** an unrelated request was delayed by about 455 ms.
- **Migration relevance:** this is one of the few problems where Go's runtime helps by default. Python fixes it with plain `def` handlers or an async driver.

### D15. Build, deploy and health integrity (Medium) [SRC]
- **Deploy is not gated on CI:** it runs on every push to `main` with no dependency on CI (`.github/workflows/deploy.yml:3-6`). Committing the nomenclature work would deploy it immediately.
- **Images ignore `uv.lock`:**
  - The API image uses `uv pip install --system -e .` (`docker/Dockerfile.api:33`).
  - The worker uses unpinned `pip install -e . rdkit` on `psi4/psi4:latest` (`docker/Dockerfile.worker:4`, `:31`).
  - The lockfile (rdkit 2026.3.4, py2opsin 1.2.0, openclatura 0.1.2) is therefore not what ships.
  - The API and worker can run different RDKit versions.
- **Rollback is `git reset --hard origin/main`** on the host (`deploy.yml:43-44`): there are no image tags or digests. `up -d` recreates the worker mid-job (feeds D4).
- **Health is always 200:** `/health` returns 200 even when the DB is down. Only the body says `unhealthy` (`routes.py:48-75`). Both `curl -f` health checks and the deploy gate (`scripts/deploy.sh:24-31`) treat it as healthy.
- **No migrations:** `alembic` is a dependency but no migrations exist. Both processes run `create_all` at startup (`queue.py:64-68`), so no schema change can be rolled out safely. This blocks §7.

### D16. Surface mesh spacing distortion (Low) [RUN]
- **Where:** `marching_cubes` is given the requested spacing, but the grid is `linspace` with `ceil(L/g)+1` points, so the actual spacing is ≤ g and differs per axis. Clamping recomputes spacing from the x-axis only (`surface_mesh.py:92-103`, `:132-138`).
- **Measured (C7):**
  - Single carbon: mean radius 3.112 Å against 3.100 expected; centroid offset 0.043 Å.
  - Clamped two-atom case: centroid offset 0.080 Å.
- The effect is small, but any Go or numeric port must decide whether to replicate it.

### Lower-order confirmed facts (no ranking)
- **Duplicate orbital entries fragment the key:** `orbitals:["homo","homo"]` and `["homo"]` produce different keys (C9).
- **`homo-1` implies HOMO:** requesting `homo-1` also computes and returns `homo` (`psi4_runner.py:156-157`).
- **Ignored unit:** `potential.units` is in the key but ignored. Output is always kcal/mol/e (`electrostatics_routes.py:91-95`, `electrostatics.py:208-211`).
- **Writes on read paths:** conformer `GET` writes to the job row (`conformers_routes.py:82`). Every cache hit writes `hit_count` (`queue.py:333-334`).
- **Mixed DB failures reported as 400:** orbital submit maps any exception, DB outage included, to 400 (`routes.py:102-104`).
- **Weak deploy host-key check:** the SSH deploy uses `StrictHostKeyChecking=no` even though it provides `known_hosts`, and it ships the production `.env` over that channel (`deploy.yml:33-35`). See R6 in §4.

---

## 4. Plausible risks (needs a check to confirm or dismiss)

| # | Risk | What would confirm or dismiss it |
| --- | --- | --- |
| R1 | Golden water is detected as C2v, so HOMO/LUMO `energyEv` are silently missing (D13) | Run the pinned worker image on `tests/golden/water.sdf` and inspect `energyEv` and `wfn.molecule().schoenflies_symbol()` |
| R2 | `psi4/psi4:latest` and unpinned pip drift between deploys changes numerics and the RDKit API | Compare running image digests and `pip freeze` across the last deploys |
| R3 | Production worker has no xTB, so every `xtb-gfn2` request has been served Gasteiger | `which xtb` / `import xtb` inside the running worker; `charge_model` histogram in cache entries (read-only) |
| R4 | Consumer key computation sometimes runs without `crypto.subtle`, producing 32-bit keys and cross-molecule collisions | Find where `prefetchQuantumOverlays.ts` and `electrostaticsClient.ts` execute (browser, secure context or not; Node version); look for 1–8 hex-char keys in `jobs.request_params->>'cache_identity'` |
| R5 | Torn `result.json` read triggers the invalidation cascade (D3 + D1) | Disposable Postgres fault injection: a slow writer and a concurrent poller |
| R6 | MITM on deploy SSH exposes the production `.env` | Threat-model review; drop `StrictHostKeyChecking=no` |
| R7 | Real backlogs exceed 900 s, so jobs that never ran are failed (D4) | Production histogram of queue wait times and of "Stale pending job auto-cleared" errors |
| R8 | OPSIN verification (two JVM launches per request, no timeout) exceeds the consumer's 10 s timeout under concurrency, leaving orphaned threadpool work | Benchmark with Java at concurrency 1/8/40; check `subprocess.run` without a timeout in `py2opsin.py:137` |
| R9 | The nomenclature tempfile monkeypatch (`nomenclature.py:60-88`) breaks silently if the image installs a different py2opsin than locked (D15) | Pin installs to the lockfile; add an image-level test of `_safe_py2opsin` |
| R10 | Event-loop stalls (D14) materially hurt p99 in production | p99 of `/metrics/cache` or `/health` during conformer `waitMs` bursts |
| R11 | Unbounded growth: no `max_length` on `sdfContent`/`molblock2d` (`schemas.py:11-16`, `:264-269`, `:428-432`), SDF copied into every job row, cache directory never evicted (cubes can be tens of MB) | Host disk usage and DB table size trends |
| R12 | Head-of-line blocking: one worker container (`docker-compose.single.yml:42-62`, fixed `container_name`), FIFO across job types, so a millisecond conformer waits behind multi-minute Psi4 jobs | Queue-wait percentiles by `job_type` |
| R13 | Concurrent `create_all` at first boot (API and worker) races on an empty DB | Disposable compose up on an empty volume |
| R14 | Starlette echoes the request origin when `allow_origins=["*"]` and `allow_credentials=True` and cookies are present (`main.py:31-37`) | Low today, since the API has no cookie auth; dismiss if it stays cookie-free |

---

## 5. Consumer compatibility matrix (consumer @ `35dd5f4`, local source, not production proof)

| Endpoint | Consumer sends (source) | Server outcome (R6 probe) | Classification |
| --- | --- | --- | --- |
| `POST /v1/orbitals/compute` | `sdf_content`, `charge`, `multiplicity`, `method`, `basis`, `orbitals`, `isovalue`, `grid_spacing`, `density_isovalue`, `client_cache_key` or `smiles`/`inchi_key` (`orbitals/schema.ts:36-61`; `presentations/cache/route.ts:164-181`) | Accepted. Snake_case is accepted via `populate_by_name`. **Dropped:** `charge`, `multiplicity`, `density_isovalue`, `smiles`. | Confirmed drift (D8, D9) |
| Orbital response | Reads `cache_key ?? cacheKey`, `job_id ?? jobId`, … and `meta.grid_spacing_angstrom` (`orbitals/schema.ts:103-158`) | Top level is fine (camelCase emitted). `meta.gridSpacingAngstrom` is never read. | Confirmed (cosmetic) |
| Orbital key | SHA-256 of `overlay-v3` + SDF text + charge/multiplicity/method/basis/grid/iso/densityIso (`cacheKey.ts:54-85`) | Used as molecule identity. The server's own key adds method/basis/grid/iso/orbitals. | Fragmentation only (charge variants computed identically); trust issue D5 |
| `POST /v1/electrostatics/compute` | `sdfContent`, `method` (default `gasteiger`), `charge`, `multiplicity`, `no_reorient`, `no_com`, optional `surface{gridSpacing,probeRadius,padding}`, optional `potential{…}`, `client_cache_key` (`electrostatics/schema.ts:44-60`; `presentations/cache/route.ts:214-226`) | Accepted. **Dropped:** `no_reorient`, `no_com` (harmless). When `surface` is omitted the server uses 1.4/3.0 while the key assumes 0.9/1.2. 422 for consumer-valid inputs: other methods, padding < 1, \|charge\| > 10. | Confirmed drift (D8) |
| `POST /v1/conformers/compute` | `molblock2d`, `geomVersion` (= consumer `method`), `waitMs` (2000), `seed`, `numConformers`, `optimize`, `forceField` (`mmff`), `clientCacheKey`, `smiles`, `inchiKey` (`conformers/schema.ts:25-38`) | **Dropped:** everything except `molblock2d`, `geomVersion`, `waitMs`. The server runs `etkdg_v3` + `uff`. | Confirmed drift (D8): consumer intent not honored |
| Conformer response | `status`, `jobId`, `cacheKey`, `sdf3d`, `errorMessage` | Statuses are `pending\|done\|failed`; field names match | Compatible |
| `POST /v1/nomenclature/explain` | `smiles`, `molblock?`, `verifyOpsin` (default true), 10 s abort (`nomenclature.server.ts:19-39`) | Accepted with no dropped fields. Response fields match `NomenclatureExplanation` (`shared/types/index.ts:265-279`). | Compatible. Timeout and abandonment risk R8. |

---

## 6. Six-lens sweep

Duplicates are merged into the defects above. Each lens lists what is broken, then what already holds.

### 6.1 Invariants: consistency and atomicity
- **Broken:**
  - I-a, at most one active job per key: D2.
  - I-b, DONE ⇔ valid CacheEntry ⇔ complete artifact set: D1, D2, D3, P2.
  - I-c, artifact set is all-or-nothing: D3.
  - I-d, terminal states absorbing: D1.
- **Already holds:** `claim_next_job` uses `FOR UPDATE SKIP LOCKED` (`queue.py:554-561`), so two workers cannot claim the same PENDING row. The defect is everything after the claim.

### 6.2 Authority: sources of truth and cache validity
- **Broken:**
  - Identity is delegated to the client (D5).
  - The key omits engine inputs: M CHG and bonds for Gasteiger (D10).
  - The key claims parameters the engine does not honor: xTB vs Gasteiger, multiplicity, charge, densityIsovalue, conformer options (D8, D10).
  - There is no engine or identity version in the key, so fixed engines will serve old results.
  - Invalidation rewrites job history rather than recording cache validity (`queue.py:270-309`).
  - Cache metrics are per process and reset on start (`main.py:53`). That is acceptable as documented.
- **Already holds:** the orbital geometry identity is order-canonical (`sdf_parser.py:57-70`, tested). The conformer key covers the full input bytes plus parameters. Electrostatics includes charge and multiplicity.

### 6.3 Ownership: resources and background work
- **Broken:**
  - No owner proves completion: there is no lease (D4).
  - Psi4 runs in-process and cannot be timed out.
  - SIGKILL leaks `psi4_scratch/<job_id>`; the `finally` at `worker/main.py:171-179` does not run.
  - `psi4.core.clean()` is called only on success (`psi4_runner.py:304`).
  - The Psi4 options singleton persists across jobs. It is low risk today because every option is re-set.
  - Nomenclature runs JVM subprocesses with no timeout inside the API threadpool, and consumer aborts do not cancel them (R8).
  - Deploy recreates the worker mid-job.
- **Already holds:** the py2opsin temp-file race is fixed correctly. Openclatura looks up `py2opsin.py2opsin` at call time (`openclatura/opsin_verify.py:136`), so the patch applies, and `tmp_fpath` is positional argument 8, matching the `len(args) >= 8` check (py2opsin 1.2.0 signature).

### 6.4 Failures: truthful outcomes and preserved data
- **Broken:**
  - Partial results become DONE (D12).
  - Fallback chemistry is cached as if it were the requested chemistry (D10).
  - Unknown elements become carbon (D10).
  - A degenerate conformer is cached (D11).
  - A successful compute is marked failed by the IntegrityError (D2).
  - If the DB fails during `complete_job`, `fail_job` fails too, the outer loop sleeps, and the job stays RUNNING with artifacts present (`worker/main.py:181-183`).
  - DB outage is reported as 400 (`routes.py:102-104`).
  - Health is 200 when unhealthy (D15).
  - The golden test accepts missing energies and density (D12).
- **Already holds:** a missing `result.json` on a DONE job becomes an explicit `error` plus recompute rather than an empty success (`routes.py:187-198`, tested). The nomenclature route separates 422 domain errors from 500 engine errors (`nomenclature_routes.py:59-74`).

### 6.5 Semantics: meaning across boundaries
- **Broken:**
  - Basis naming (D6).
  - Charge and multiplicity (D8, D9, D10).
  - Energy index domain (D13).
  - SDF dialect coverage (D7).
  - Density isovalue (D8).
  - Conformer parameter vocabulary (D8).
  - Default surface parameters (D8).
  - Response key casing (D8).
  - Spacing assumption (D16).
- **S4, Go parity hazards for identity:** `round(v, 6)`, then negative-zero normalization, then `%.6f` formatting, then a sort of `(symbol, x, y, z)` tuples, then SHA-256 truncated to 16 hex characters (`sdf_parser.py:49-70`). In C4, Python maps −1e-7 to `0.000000`, while a naive Go `strconv.FormatFloat(-1e-7,'f',6,64)` yields `-0.000000` [KNOW]. Ties sort on rounded floats. `compute_cache_key` formats floats with `%.4f` and JSON-dumps with `sort_keys`, and Go's `encoding/json` spacing and float rendering differ from Python's `json.dumps`. Any Go key function needs a golden vector corpus.
- **Already holds:** Ångström-to-Bohr conversions are consistent (`psi4_runner.py:24`, `:150`; `cube_parser.py:13`; `charges.py:80-81`). Cube-to-world transforms are correct for orthogonal Psi4 grids (`marching_cubes.py:136-156`). Mesh encoding is float32/uint32, gzip with mtime=0, base64, and deterministic (tested).

### 6.6 Bounds: capacity and load amplification
- **Broken:**
  - Unbounded request bodies, queue and cache directory (R11).
  - One FIFO worker for all job types (R12).
  - Age-based failure of pending jobs (D4).
  - Event-loop stalls (D14).
  - `waitMs` polling at 10 Hz per waiting request.
  - Every poll of a DONE job re-reads and re-validates the full mesh JSON through Pydantic (`routes.py:235-306`).
  - Nomenclature has no concurrency cap beyond the default threadpool, plus two JVM spawns per verified call (R8).
  - Write amplification on cache hits.
- **Already holds:** surface grids are clamped to 200³ (`surface_mesh.py:96-103`). Nomenclature inputs are bounded (4096 characters of SMILES, 2 MB molblock, 256 atoms; `schemas.py:159-164`, `nomenclature.py:46`). Consumer-side timeouts exist (8–10 s).

---

## 7. State models

### 7.1 Current implicit job lifecycle (as built; accidental behavior marked)

**Owners:** API = submit and read routes. W = worker loop. Nobody owns the reaper; it is piggybacked on submit.

| Current | Event | Owner | Guard | Next | Effect / note |
| --- | --- | --- | --- | --- | --- |
| — | submit, valid CacheEntry, artifacts exist, a DONE row exists | API | none beyond read | (returns DONE job) | `hit_count++` (write on read) |
| — | submit, CacheEntry but artifacts missing | API | — | DONE→ERROR for all key rows; entry deleted | then continues to the active check |
| — | submit, no active row | API | read only, **no lock or unique constraint** | PENDING | duplicates possible (D2) |
| PENDING/RUNNING | submit, row older than 900 s | API | age only, **ignores the live owner** | ERROR | new PENDING created (zombie risk, D1) |
| PENDING | claim | W | `status='pending'` + SKIP LOCKED; **if age > 900 s, set ERROR** | RUNNING | attempt++ (never used); backlog becomes errors (D4) |
| RUNNING | engine returns | W | **none** (no status or owner check) | DONE | upsert CacheEntry; races cause IntegrityError (D2) |
| RUNNING | engine raises | W | none | ERROR | — |
| ERROR | late `complete_job` | W (zombie) | **none** | DONE | reproduced (P1) |
| DONE | late `fail_job` | W | **none** | ERROR | orphan CacheEntry (P2) |
| DONE | GET, `result.json` missing or unparseable | API | — | ERROR (all DONE rows of key) | entry deleted |
| RUNNING | worker crash, OOM or deploy | — | — | **RUNNING indefinitely** | cleared only by a later same-key submit after 900 s |
| QUEUED | — | — | — | — | dead state |

### 7.2 Minimal authoritative model (proposed; enforced in Postgres so Python and Go can share it)

Keep two records and one orthogonal dimension. Do not flatten them into one large enum.

**`jobs`:** `id`, `kind`, `cache_key`, `key_version`, `effective_spec` (the parameters the engine will actually honor), `state ∈ {queued, running, succeeded, failed}`, `attempt`, `lease_token`, `lease_expires_at`, `heartbeat_at`, `error_code`, `error_retryable`, `result_ref` (cache_key + publication id), timestamps.

**`cache_entries`:** `cache_key` (PK), `publication_id`, `manifest` (files, sha256, bytes), `engine_versions`, `effective_params`, `quality` (`exact | degraded:<reason>`), `published_at`, `invalidated_at NULL`, `producer_job_id`.

**Publication sub-state, on disk and owned by the worker or supervisor:** `staging` → `sealed` (manifest and checksums written) → `published` (atomic rename to `cache/<key>/<publication_id>/`) → `referenced` (DB row). It is orthogonal to the job state.

**Invariants (DB-enforced where possible):**
1. At most one `queued|running` job per `cache_key`: `CREATE UNIQUE INDEX … ON jobs(cache_key) WHERE state IN ('queued','running')`. Submit does `INSERT … ON CONFLICT DO NOTHING RETURNING` and falls back to selecting the existing row.
2. Only the lease holder transitions a running job: `UPDATE jobs SET … WHERE id=$1 AND state='running' AND lease_token=$2`. Zero rows updated means stale: the caller discards its staging dir and records nothing.
3. `succeeded` implies the referenced publication exists with matching checksums. CacheEntry is written in the same transaction with `ON CONFLICT (cache_key) DO NOTHING`. If an entry already exists and is valid, the job still succeeds and references it.
4. `succeeded` and `failed` are absorbing. Invalidation sets `cache_entries.invalidated_at` and never rewrites jobs.
5. Recovery is driven only by lease expiry, never by queue age. Optional client-visible expiry is a separate, explicit field (a product decision).
6. `attempt ≤ max_retries`. Exhaustion leads to `failed(retries_exhausted)`.
7. A `degraded` result is never reused for a request whose spec demands `exact`. The alternative is that the degraded result gets its own key; that choice is the product rule U2 in §10.
8. `key_version` changes whenever identity or engine semantics change. Old entries stay readable for rollback.

| Current | Event | Owner | Guard | Next | Effect |
| --- | --- | --- | --- | --- | --- |
| — | submit | Edge (API) | spec validated strictly; valid non-invalidated entry with compatible `quality` | (answer from cache) | best-effort hit counter (async or batched) |
| — | submit | Edge | no usable entry; `INSERT … ON CONFLICT` on the partial unique index | queued (new or existing row) | idempotent dedupe |
| queued | claim | Supervisor | `state='queued'` + `FOR UPDATE SKIP LOCKED` | running | new `lease_token`, `lease_expires_at=now+L`, `attempt++` |
| running | heartbeat | Supervisor | token matches and running | running | extend lease |
| running | engine ok, publication sealed and renamed | Supervisor | token matches, running, checksums verify | succeeded | same transaction: entry upsert `ON CONFLICT DO NOTHING`; `result_ref` set |
| running | engine error (retryable) | Supervisor | token matches; `attempt < max` | queued | clear lease; record error |
| running | engine error (non-retryable) or timeout policy | Supervisor | token matches | failed | `error_code` (`engine_error`, `timeout`, `invalid_input`, `unsupported`) |
| running | lease expired | Reaper | `lease_expires_at < now()` and running | queued, or failed(`retries_exhausted`) | clear lease. The old worker's later writes fail invariant 2. |
| running | complete/fail with a stale token | Zombie worker | token mismatch | unchanged | rejected; worker deletes its staging dir |
| succeeded | artifact missing or corrupt on read | Edge | entry exists | succeeded (job unchanged); entry `invalidated_at=now()` | the read returns `error: result_invalidated`; the next submit creates a new job |
| succeeded/failed | any worker event | — | — | unchanged | ignored and logged (absorbing) |
| queued | cancel | ? | **product rule needed** | — | currently undefined; do not invent it |

**Failure and recovery of the publication sub-state:**
- Crash in `staging`: an orphan staging dir, removed by an age-based GC; the lease expires and the job is requeued.
- Crash after rename, before the DB commit: an orphan `published` dir that no row references. GC it after a grace period, or let a retry adopt it after checksum verification. Publication ids make "overwrite" unnecessary.
- DB commit failure after rename: same as the previous case.
- Readers resolve only `cache_entries.publication_id`. They never see a partial directory.

**Nomenclature is stateless request/response** and does not need this model. Its required contract is a server-side deadline below the consumer's 10 s, a concurrency cap, and a distinct `opsin_status` for timeout versus error.

---

## 8. Migration options and recommendation

| Option | Description | Addresses D1–D16? | Main risks | Verdict |
| --- | --- | --- | --- | --- |
| A. Harden in Python | Stages 0–3 below, no Go | Yes, all are fixable | Keeps one language. Python edge concurrency needs care (D14). | **Mandatory first, whatever is chosen** |
| B. Go control plane + Python engines (strangler on a shared Postgres authority) | Go edge, queue control, artifact serving and supervisor; chemistry stays Python | Yes, plus strict decoding and process isolation by construction | Two languages; identity parity; shadow-diff effort | **Recommended *if* stage 0 data or operational goals justify it, after A** |
| C. Port the API only to Go, keep the worker as is | Thin rewrite | No. Ports D1–D5, and splits transition logic across languages. | Split-brain transitions; key drift | Reject |
| D. Full Go, chemistry included (cgo RDKit MinimalLib; subprocess Psi4/xTB/OPSIN) | Rewrite everything | Re-creates the chemistry layer | No Go RDKit package found; Psi4 has no Go API; large numeric parity surface | Reject |
| E. Python API + Go engine supervisor only | Go only for subprocess supervision and leases | Fixes the D4 ownership class | Small surface; still needs the §7 SQL | A viable smaller alternative to B |

### Stages: acceptance and rollback gates, ordered by dependency, with no deadline

**Stage 0. Evidence baseline (no behavior change)**
- **Do:**
  - An operator runs read-only queries for job-state histograms, queue-wait and run-time percentiles by `job_type`, top error messages, `charge_model` vs requested method, cache-directory size, running image digests, and which compose file production uses.
  - Record the five consumer request shapes as contract fixtures.
  - Capture a privacy-reviewed corpus of real SDF inputs and their current server keys.
- **Accept when:** the numbers are recorded, and the fixtures run against the current server, showing the D8 drift list as red.
- **Rollback:** none needed.

**Stage 1. Two-way contract agreement (consumer ↔ server)**
- **Do:**
  - For each D8 item, decide "implement" or "reject with a coded 422": orbital `charge`, `multiplicity`, `densityIsovalue`; conformer `optimize`, `forceField`, `seed`, `method`; electrostatics default surface; meta casing.
  - The server logs unknown fields first (observe mode), then switches to `extra="forbid"` per model.
  - The consumer stops sending loose identity fields and aligns its defaults.
  - Decide the trust model for `clientCacheKey` (U3 in §10): network-restrict, sign with HMAC, or derive identity on the server.
- **Accept when:** both repos' contract tests pass, and no unknown-field events appear over an agreed observation window.
- **Rollback:** a config flag back to `ignore`. Consumer changes are backward compatible because they send a subset.

**Stage 2. Authoritative state model in Postgres (still Python)**
- **Do:**
  1. Add an Alembic baseline that matches the production schema; verify it as a no-op.
  2. Run an additive migration for the §7.2 columns.
  3. Run a data migration that collapses existing duplicate active rows.
  4. Add the partial unique index.
  5. Implement the conditional transitions, lease, heartbeat and reaper, wire up `max_retries` and the timeout, and add staging, rename and manifest publication.
  6. Upsert CacheEntry `ON CONFLICT`.
  7. Record invalidation on the cache entry, not on jobs.
- **Accept when** a fault-injection suite on disposable Postgres passes:
  - 50 concurrent identical submits leave exactly 1 active row.
  - A zombie completion is rejected.
  - `kill -9` mid-job leads to a requeue after lease expiry and exactly one success.
  - A reader never observes a partial publication.
  - A double completion leaves both jobs succeeded and one entry.
  - A backlog older than the old 900 s threshold is not failed.
  - Then canary production metrics show 0 duplicates and 0 orphans.
- **Rollback:** migrations are reversible. The partial index must be dropped before reverting code, because old code would hit IntegrityError on duplicate inserts. Keep the old claim path behind a flag for one release.

**Stage 3. Engine semantic fixes (Python), with `key_version` bump**
- **Do:**
  - D6: canonical basis passed to Psi4.
  - D9: charge and multiplicity to Psi4, if adopted in stage 1.
  - D7: parser fixes, or parse with RDKit in the worker.
  - D10: explicit fallback policy, `uhf`, an RDKit periodic table or rejection for unknown elements, M CHG and bonds in the electrostatics identity.
  - D11: keep the 2D coordinates on failure.
  - D12: fail, or emit an explicit `partial` quality.
  - D13: `epsilon_a_subset("MO","ALL")` or `symmetry c1`.
  - D16: decide whether to fix it, since it changes outputs.
  - Pin images to lockfile and digest. Gate deploy on CI.
- **Accept when:**
  - CI runs Psi4 integration on the pinned worker image with assertions required, not optional, over: symmetric water, a C1 molecule, 6-31G*, a cation, a radical, a ≥100-atom SDF, a blank-title SDF, a metal (explicit rejection), and a forced embed failure.
  - A chemist reviews the before/after diffs.
  - Keys for currently accepted inputs are unchanged, except where `key_version` is intentionally bumped.
- **Rollback:** previous image digest plus previous `key_version`. Old entries remain because the keys differ.

**Decision point: Go or not.** Proceed to B (or E) only for a stated reason, for example:
- measured edge latency or throughput problems after D14 is fixed,
- operational need for subprocess supervision,
- the team wanting Go for the control plane,
- edge-security posture.

If none applies, stop at A. It is a complete, defensible outcome.

**Stage 4. Go edge, read-only (strangler)**
- **Do:** the Go service serves `GET /v1/*/jobs/{id}`, artifacts (with ETag and Range), readiness with correct 5xx, and metrics. A Caddy path split sends that traffic to Go, and Python handles shadow comparisons. Compare *parsed* JSON, because float text differs between `encoding/json` and `json.dumps`.
- **Accept when:** zero semantic diffs over the window, p99 is no worse, and error-rate parity holds.
- **Rollback:** flip the Caddy route back. The DB is unchanged by Go.

**Stage 5. Go submissions, one job type at a time**
- **Do:**
  - Conformers first: the key is SHA-256 over input bytes plus `json.dumps(params, sort_keys=True, separators=(",",":"))` (`conformers_routes.py:27-39`), which is trivially reproducible with test vectors.
  - Electrostatics and orbitals only after the identity function reaches **100% parity** on the stage-0 corpus plus fuzzed SDFs (negative zero, ties, formatting, CRLF, trailing data), *or* after identity moves engine-side (U5 in §10).
  - All transitions use the same SQL statements, shared as a versioned SQL file, so Python and Go cannot diverge.
- **Accept when:** key parity is 100%, dedupe and hit-rate metrics are unchanged, and the shadow diff is 0.
- **Rollback:** route flip.

**Stage 6. Engine supervisor (Go or Python)**
- **Do:** run one-shot engine subprocesses (`spec.json` in, `manifest.json` out). The supervisor enforces a hard timeout with process-group kill, heartbeats the lease, publishes per §7.2, and adds separate lanes or priorities per `job_type` to remove head-of-line blocking (R12).
- **Accept when:** kill, timeout and OOM tests pass, and throughput stays within an agreed tolerance including per-job interpreter and Psi4 import cost. That cost is unmeasured; measure it here.
- **Rollback:** restore the worker entrypoint to the in-process loop, which shares the same DB contract.

**Stage 7. Nomenclature as an engine sidecar (Python)**
- **Do:** run Openclatura, RDKit and OPSIN in their own container with a persistent JVM or bounded pool, a server deadline below 10 s, and a concurrency cap. The edge proxies to it.
- **Near-term advice for the uncommitted work:** before committing it into the API image, at least add a deadline and a cap. Prefer giving it its own container now, since it adds RDKit and a JRE to the API tier, which moves against a thin edge.
- **Accept when:** the OPSIN-verified benchmark meets its gate in an environment with Java, and p99 stays under the consumer timeout at the expected concurrency.
- **Rollback:** route back to the in-API implementation.

**Stage 8. Retire FastAPI** after an agreed window with zero fallback-route traffic.

---

## 9. Assessment of the uncommitted nomenclature work

**What is sound:**
- Explicit separation of the engine and viewer atom domains, with mapping only after a full-graph canonical-SMILES equality check (`nomenclature.py:181-209`).
- Bounded inputs and guarded domains for metals, mixtures and size.
- A correct per-call OPSIN temp file (verified against Openclatura's call site).
- Clear 422 vs 500 mapping.
- 8 unit tests passing.
- A reproducible corpus with provenance.
- Local baseline: 398 cases, 100% supported-name rate, 100% guardrail recall, about 18 ms per case without OPSIN.

**Gaps:**
- The OPSIN gate is unmeasured here and has no Java in the CI unit job, so the README's release gate is not enforced anywhere.
- There is no subprocess timeout, and each verified request launches two JVMs (R8).
- The image does not install from the lockfile, so py2opsin and RDKit versions in production can differ from what was tested (D15, R9).
- Committing it to `main` deploys it immediately, without a CI gate (D15).
- Part of the corpus comes from Openclatura's own round-trip cases, so the 100% name rate is partly self-referential. The README acknowledges this.

**Boundary placement:** chemistry engine (Python), stateless, synchronous with a deadline. Do not port it to Go.

---

## 10. Disagreements and unknowns that block a sound decision

- **U1. Production bottleneck unknown.** Nothing here shows the Python edge is a bottleneck. Psi4 dominates, by inspection. "Go for performance" is unsupported until stage 0.
- **U2. Product rules.** Should charged or open-shell orbitals be supported? Should a degraded electrostatics result (Gasteiger for an xTB request) be served at all, or keyed separately? Must conformer `optimize`/`forceField`/`seed` be honored? Is there a cancel operation? Is age-based expiry of queued work desired? These need the two-way agreement and should not be invented by either implementation.
- **U3. Trust model for `clientCacheKey`.** Is `/v1/*` a public API, or only for the Next.js backend? The answer decides between network restriction, signed keys and server-derived identity, and it blocks any identity port.
- **U4. Preserving existing cache entries.** Keys may already contain poisoned entries (D5, D10), degenerate ones (D11) or partial ones (D12). Deciding whether to keep them through `key_version` changes needs a read-only audit, which was out of scope.
- **U5. Where identity is computed.** Either the edge does it, which needs SDF parsing in Go with a parity corpus, or the engine does, which makes dedupe asynchronous and changes the "cached on submit" behavior. This is the main design fork for option B.
- **U6. Team and operations.** A second runtime adds build, observability and on-call cost. That is not measurable from the repo.
- **U7. Deployed reality.** Which compose file (`DEPLOY_COMPOSE_FILE` is a secret), image digests, whether xTB exists, and real queue behavior are all unverified. R1–R4, R7, R10 and R12 depend on them.
- **U8. Psi4 execution evidence.** D6's final failure step and D13's golden manifestation need one disposable run of the pinned worker image.

**Likely reconciliation points with Codex** (I have not seen Codex's conclusions):
1. I place DB-enforced transitions (stage 2) *before* any Go code, so Go ports a correct model rather than an implicit one.
2. I treat porting the identity function to Go as gated on a parity corpus, or replaced by engine-side identity.
3. I keep all numeric mesh code in Python.
4. I consider stopping at option A a legitimate end state.

---

## 11. Evidence index (primary sources consulted)

- Psi4 `psi4/driver/qcdb/libmintsbasisset.py` (`make_filename`): https://github.com/psi4/psi4/blob/master/psi4/driver/qcdb/libmintsbasisset.py
- Psi4 `psi4/src/export_mints.cc` (Vector bindings) and `psi4/src/psi4/libmints/vector.h` (`get(int m)` returns `v_[m]`; `dim(int h=0)` returns `dimpi_[h]`): https://github.com/psi4/psi4
- Psi4 cubeprop manual (energy-order indices, `Psi_a_5_1-B1.cube` naming): https://psicode.org/psi4manual/master/cubeprop.html
- xtb-python `xtb/interface.py` (`Calculator(..., charge=None, uhf=None, ...)`): https://github.com/grimme-lab/xtb-python/blob/main/xtb/interface.py
- RDKit CFFI (MinimalLib) for other languages: https://greglandrum.github.io/rdkit-blog/posts/2021-05-01-rdkit-cffi-part1.html. The Rust `rdkit-sys` wraps it (https://docs.rs/crate/rdkit-sys/latest). No Go package was found.
