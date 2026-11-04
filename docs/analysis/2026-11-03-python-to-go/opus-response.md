# Moleculens Python-to-Go: Opus response to the reconciliation proposal

Date: 2026-11-03. This is analysis only. No runtime source edits, commits, deploys or production access.

**Inputs reviewed:**
- `/tmp/moleculens-go-reconciliation-proposal-2026-11-03.md`
- `/tmp/moleculens-go-codex-analysis-2026-11-03.md`
- My first report, `/tmp/moleculens-go-opus-independent-2026-11-03.md`. It is preserved unchanged. Where it differs from this document on migration, identity, rollback or gates, **this document supersedes it**.

Server HEAD is still `c9fce3c5d9f69212cc65d39e06ae3c18ef74eafb`, and the dirty and untracked set is unchanged (15 status lines, same files). The consumer is unchanged.

## Verdict

**Architecture, sequence and release gates: agreed**, with amendments A1–A7 plus the evidence corrections in §3. None of them changes the architecture. They tighten gate content and close a few unspecified edges.

**No material unresolved disagreement remains.** Two items are open specification choices that must be closed before Go writes production state (§4). They are not disagreements about direction.

The agreed shape:
- Go owns the HTTP boundary, admission, PostgreSQL transitions, supervision, cache identity and artifact publication.
- Pinned Python adapters own Psi4, RDKit, Openclatura/OPSIN and, initially, the mesh and ESP numerics.
- PostgreSQL remains the queue.
- Bounded correctness fencing comes before Go writes, with read-only Go and fixture work running in parallel.
- Cutover drains the legacy queue (preferred), or isolates both the queue and the artifact root.
- Release requires real PostgreSQL contention and crash tests, consumer differentials, science fixtures with tolerances, a real Java/OPSIN gate and a rollback drill.

## 1. Focused checks run for this reconciliation

Probes are in `/private/tmp/claude-501/-Users-kelly-Developer-moleculens-server/c4a0daf6-1b5c-4ef9-8c7b-20967be17449/scratchpad/probes/probe_reconcile.py`. They use a SQLite temp DB, an in-process ASGI client and no repo writes.

| # | Question | Result | Effect |
| --- | --- | --- | --- |
| K1 | Codex: installed RDKit has no `IsMetal` | `hasattr(PeriodicTable,'IsMetal')` is False. `_has_metal` returns False for `[Fe]` and for cisplatin. | Codex is correct. The metal guard in `compute/conformers.py:50-62` never fires, so UFF/MMFF runs on metal complexes. I missed this. |
| K2 | Codex: polling a conformer job via `/v1/orbitals/jobs/{id}` destroys its association | The response is `error`. The conformer job flips DONE→ERROR, the CacheEntry is deleted, and `structure.sdf` stays on disk. | Codex is correct (`api/routes.py:172-198` → `:32-45` → `db/queue.py:270-309`). I missed this. |
| K3 | Codex: the Gasteiger dipole is computed about the coordinate origin but reported about the centroid | Na⁺ at x = 10 Å: charge +1, reported origin (10,0,0), \|μ\| = 48.03 D, whereas it is 0 about the reported origin. | Codex is correct (`compute/charges.py:255-261`). |
| K4 | Mine: does MOL line 2 (2D/3D code) change what RDKit reads? | RDKit overrides a `2D` tag when z ≠ 0 ("Marking the mol as 3D"), so `Is3D` was the same either way. | My hypothesis is **not** supported for non-planar input, and I drop it as an argument. Header-byte exclusions remain unproven either way (A2). |
| K5 | Naming test count | `pytest --collect-only tests/unit/test_nomenclature.py` collects **9** | Codex is correct. My "8" was wrong; the suite total of 46 stands. |
| K6 | Per-attempt process cost for light jobs | Python with RDKit import: 0.24–0.31 s cold (×3 each, local macOS). A warm in-process conformer run takes 8.2 ms. | A per-attempt process is about 30× a conformer's compute time but fits inside the consumer's `waitMs=2000`. JVM start for OPSIN is unmeasured (no Java). Feeds A6. |

Source-only confirmations of other Codex additions I had missed, all [SRC]:
- CI's `grep -q "healthy"` also matches `unhealthy` (`.github/workflows/ci.yml:92`).
- Deploy resets to the current `origin/main`, not the triggering SHA (`deploy.yml:43-44`).
- Output and scratch `mkdir` run outside the inner failure handler, so a claimed job stays RUNNING (`worker/main.py:88-93`, `:181-183`).
- The result loaders turn `{}` into an empty success (`api/routes.py:254-306`, `api/electrostatics_routes.py:253-308`).
- Spacing and isovalue are rounded to 4 decimals in the key but not in compute (`db/queue.py:97-98`).
- `float()` accepts `nan`/`inf` coordinates (`compute/sdf_parser.py:156-158`).
- An extreme aspect ratio can truncate a grid axis below two samples (`compute/surface_mesh.py:96-103`). This produces a failure, not silent corruption.

**Ranking difference (not material):** Codex ranks K2 as P1. I would rank it P2. It needs a mis-routed or hostile poll; files survive and the next submit recomputes. Both rankings put it in the bounded hardening set, so the difference changes nothing.

## 2. The seven decisions

### Decision 1 (target and rationale): **Accept, with amendment A1**
- **Accepted:** the user's maintainability goal is a sufficient reason to pursue Go. In my first report I listed "the team wanting Go for the control plane" as a valid trigger at the decision point, so that point is now satisfied. Baseline measurement sizes and verifies the work; it is not a precondition. I withdraw any reading of my report that required a demonstrated speedup.
- **A1, admission validation:** with a byte-based descriptor (Decision 2) and no Go chemistry parser, the edge no longer parses SDF at submit. Today an invalid SDF on orbital submit returns **400**, and the integration test `tests/integration/test_golden.py:352-366` asserts it. Stage 0 fixtures must pick one, explicitly:
  - (a) a small Go structural validator: counts and atom-block shape, finite coordinates, known element symbols, size limits. It has no identity authority and is checked against an independent RDKit oracle on the fixture corpus.
  - (b) a bounded synchronous adapter validation call.
  - (c) a documented contract change: accept the submit, then fail the job with a typed `invalid_input`.

  My recommendation is (a). It keeps the 400 contract and the RDKit oracle catches drift.

### Decision 2 (new versioned identity, no mandatory legacy parity): **Accept, with amendment A2**
- **Accepted:** this is a better resolution than either fork in my first report (port with parity, or engine-side async identity), and I adopt it.
  - Hashing conservatively normalized full bytes automatically includes bonds and `M CHG`, closing my D10 key-omission finding.
  - Binding the client hint to the server digest closes D5.
  - Legacy reads stay under their recorded identity.
  - Exact parity is needed only for an intentional legacy import path.
- **Known cost, accepted:** reuse drops across equivalent serializations. Examples: atom-order permutations, which today's orbital `geometry_hash` canonicalizes; differing title or SD data fields; writer timestamps on line 2.
- **A2, version granularity (open, see §4):** "engine/adapter versions" in the descriptor should be a declared **result-compatibility version** per adapter. It is bumped deliberately, and a CI output-drift fixture forces the bump when outputs change. The exact engine build (image digest, package versions) goes into the manifest as provenance, not into the key. Otherwise every image rebuild cold-starts the whole cache, and Psi4 recompute costs minutes per entry.
- **A2, later normalization:** any future exclusion of non-engine bytes (title, line-2 program/timestamp, SD data items after `M  END`) needs its own descriptor version and fixtures proving that no adapter reads those bytes. K4 neither proves nor refutes that.

### Decision 3 (ordering and bounded scope): **Accept, with amendments A3 and A4**
- **Accepted, and I withdraw my first report's stage 2 as written** (full lease, heartbeat, reaper and atomic publisher implemented in Python before any Go). The bounded alternative is right: Python gets only the fencing needed to stay safe while it remains the writer, and the full supervisor and publisher are built once, in Go.
- **Minimum Python fencing** (matches Codex stage 1):
  - an active-key partial unique index with `ON CONFLICT`,
  - attempt-token-conditional complete and fail,
  - an operation-kind check before any read-path mutation (K2),
  - no age-based failure of PENDING rows,
  - attempt-private output directories,
  - result truth: no empty or partial success.
- **A3, strict-mode activation order:** the current consumer *always* sends values that differ from server behavior: `density_isovalue: 0.01` (the server uses `isovalue`), `forceField: 'mmff'` (the server runs UFF), and an omitted `surface` (the key assumes 0.9/1.2; the server uses 1.4/3.0). Turning on "reject unsupported non-default" first would 422 every density and conformer request. "Default" must mean the server's *effective* behavior. These three must be implemented or changed on the consumer side, and unknown-field logging must run in observe mode against consumer fixtures, before reject mode. The consumer differential gate in Decision 6 enforces this; I state it so nobody reads Decision 3 as "enable strict first".
- **A4, one schema-migration owner:** today there are no migrations and both processes run `create_all`. The fencing needs schema changes before Go exists. Choose one versioned migration mechanism for the shared schema now (plain versioned SQL applied by one runner, usable later from Go). Do not introduce Alembic for Python and a separate Go tool later. Remove `create_all` from both process startups.

### Decision 4 (state and publication): **Accept**
- Agreed:
  - attempt tokens, database-time leases and conditional SQL;
  - absorbing terminal outcomes, with artifact availability as a separate dimension;
  - attempt-private immutable bytes before the publication pointer;
  - only the current attempt may publish;
  - idempotent re-acknowledgement for the same attempt and manifest;
  - kind-mismatch reads that never mutate.
- Codex's handling of an unknown commit acknowledgement is stricter than mine and correct. Read back the job, attempt and manifest identity. If the DB is unavailable, hold the outcome as ambiguous, and never garbage-collect possibly-referenced bytes.
- I withdraw my renamed states (`queued/succeeded/failed`). Keep `pending/running/done/error` internally and the conformer `pending/done/failed` projection, unless a separate compatibility change is justified.

### Decision 5 (cutover and rollback): **Accept, with amendment A5**
- **Accepted, and I withdraw two rollback statements** from my first report:
  - stage 2: "drop the partial index before reverting code";
  - stages 4–5: "route flip back is instant because the DB is unchanged by Go".

  Both are unsafe once writers differ. The old worker claims every pending row with no protocol guard, so neither an unfenced rollback nor a route flip for write paths is safe without verified compatibility.
- Rollback targets a **hardened, fenced** compatible image. It never drops uniqueness or fencing, and accepted jobs stay readable and drainable. Route flips stay acceptable for Go read paths that do not mutate state. Note that the legacy Python read paths themselves mutate (invalidation, `hit_count`, the conformer log flag), so flipping reads back to Python restores those legacy writes.
- **A5, drain must terminate:** today's legacy worker has no timeout. A hung legacy RUNNING row never drains. Before the cutover window, either the hardened Python image carries a deadline, or an operator cutoff rule fences the remaining rows to a typed failure (`legacy_drain_timeout`, guarded by attempt token). The drain exit criterion: zero legacy pending/running rows, the legacy worker stopped, and its image retained.

### Decision 6 (qualification and operation): **Accept, with amendment A6**
- Agreed:
  - PostgreSQL remains the queue; no new framework.
  - Bounded slots, process-tree stop and reap.
  - One owner for deadlines, retries, backpressure, output budgets and retention.
  - Real PostgreSQL contention and crash tests, consumer differentials, science fixtures with tolerances, capability and readiness tests, and a rollback drill.
  - The naming oracle, atom-mapping protections, the 398-case corpus and the real Java/OPSIN gate are all preserved.
- **I withdraw "a chemist reviews the before/after diffs"** as a gate. The replacement is an evidence rule, not an approval workflow: fixture expected values must come from a source independent of the change under test. That means analytic cases (sphere radius, point-charge dipoles about a declared origin, as in K3), a documented reference-engine run, or literature values. They are authored separately from the implementation and state their tolerances.
- **A6, execution classes:** per-attempt processes are agreed for Psi4 and ESP. For conformers and naming, the gate must measure per-attempt start-up against the consumer budgets in the *target image*:
  - Conformers: K6 shows about 0.3 s against `waitMs` 2000.
  - Naming: with `verifyOpsin`, the JVM start is unmeasured, against a 10 s consumer abort.

  If start-up dominates, use a bounded pool of warm adapter processes. Same isolation and credentials boundary, just long-lived.

### Decision 7 (product policy explicit, not an architecture blocker): **Accept, with amendment A7**
- Agreed: no calendar or performance target. Optional numeric or native ports need measured benefit plus independent correctness evidence.
- **A7:** two items in the policy list are partly defects and must not wait for policy.
  - **Planar conformer.** Whether to return a labeled 2D fallback or fail is a policy choice. Writing every atom at (0,0,0) under `planar_fallback` is never valid and is fixed regardless. Preserve the input 2D coordinates, or fail.
  - **Strict vs disclosed fallback.** Until that policy is recorded, a fallback result must not populate the key of the strictly requested method. Decision 2's descriptor includes the fallback policy, so the safe default is: fallback output is either keyed as fallback or not cached.

## 3. Evidence corrections accepted (applied here; first report left unchanged)

1. **Naming tests:** 9, not 8 (K5). The suite total of 46 is unchanged.
2. **xTB in production (my D10):** my sentence "in production every xTB request is Gasteiger" is categorical and is withdrawn. The correct statement: source packaging indicates xTB is not intentionally installed (no install step in `docker/Dockerfile.worker`, its comment, commit `5b33b80`). Whether the running image provides xTB, for example through the `psi4/psi4` conda base, was not inspected. The share of xTB requests served by Gasteiger is therefore unknown (my R3 stands as the check).
3. **Go ecosystem:** "no Go RDKit binding found" is a search result, not proof of absence. My first-report phrase "there is no native Go … xTB" was also misleading, because xTB has a documented C API that Go could reach through cgo. The recommendation rests on keeping established engines and controlling validation scope, as the proposal states.
4. **Psi4 manifestations:** D6 (the basis rewrite cannot resolve) and D13 (energy/irrep indexing) are confirmed from source plus Psi4 primary source. Their concrete manifestation in the pinned native image is unproven. I agree with that boundary.
5. **Attribution:** my queue replay, degenerate-fallback and event-loop probes remain mine. Codex's 14 reproduction assertions remain Codex's. SQLite runs are not PostgreSQL concurrency qualification, a limitation both reports share.

## 4. Open items to close before Go writes production state (not disagreements)

| # | Item | Owner and gate | Why it is not a blocker to agreement |
| --- | --- | --- | --- |
| O1 | Descriptor version granularity (A2): result-compatibility version in the key, build provenance in the manifest, CI drift detection | The descriptor spec, frozen with cross-language vectors before Go stage 2 writes | Either answer is safe. It trades cache cold-starts against stale-result risk, and the drift fixture makes the declared-version option safe. |
| O2 | Admission validation contract (A1): Go structural validator, adapter call, or accept-then-fail | Stage 0 contract fixtures, including `test_invalid_sdf_returns_error` | It changes HTTP status semantics, not ownership or sequence. |

Also carried forward as plain unknowns: production bottlenecks, deployed images and xTB, real PostgreSQL behavior, the Java/OPSIN round-trip rate, the Go runtime, and the `clientCacheKey` trust model (Decision 2's digest binding makes this non-blocking).

## 5. Agreed sequence and gates (merged, for the record)

| Stage | Content | Release gate |
| --- | --- | --- |
| 0 | Freeze contracts, consumer fixtures, evaluator corpus and engine and workload baselines. Decide O2. | Fixtures cover all four capabilities, aliases, status and error bodies, and supported vs unsupported chemistry. Drift list recorded. |
| 1 | Bounded Python fencing and result truth (Decision 3 minimum). One migration owner (A4). Pin images and lockfile. Gate deploy on CI and the event SHA. Readiness with correct status codes. Resolve the three always-sent mismatches (A3). Degenerate-conformer fix (A7). | Defect reproductions fail before and pass after. Real PostgreSQL race and crash tests pass. Corrected science cases use independently sourced expectations. Rollback goes to the previous *fenced-compatible* build. |
| 1∥ | Go read-only paths and the contract harness, in parallel | Parsed-JSON differential against Python. No mutation on read. |
| 2 | Go HTTP boundary and admission with the new descriptor (O1 frozen). Bounded naming bridge (A6). | Consumer differentials, no silent defaults, overload and readiness behavior. |
| 3 | Go supervisor owns claim, lease, retry, deadlines, process-tree stop, and immutable publication. Python adapters hold no DB credentials. | Two-worker, kill, late-result, lost-acknowledgement, full-disk, corrupt-artifact, restart and saturation tests on PostgreSQL. |
| 4 | Drain-first cutover (A5 termination rule), or overlap with isolated queue *and* artifact roots plus job-origin routing | Chemistry and contract gates, measured operating metrics, a rollback drill to the hardened image, and legacy job IDs readable after cutover and after rollback |
| 5 | Optional numeric or native ports | Measured benefit plus independent numerical and topology parity |

## 6. Statement for the two-way record

Opus agrees with Codex's reconciled architecture, sequence and release gates as set out in the proposal, with amendments A1–A7 above and the evidence corrections in §3. There is no material unresolved disagreement. O1 and O2 are specification choices with agreed owners and gates. If Codex rejects A1, A2 or A3 outright, rather than amending them, I would treat that as a disagreement on gate content to resolve before stage 2, not as a disagreement on direction.

This agreement covers a plan. It is not release approval, performance proof, or evidence that any stage has passed.
