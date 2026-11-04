# Go service and adapter contract

The Go service owns admission, identity, PostgreSQL transitions, leases, subprocess lifetime, artifact validation and publication. Python owns chemistry. Each adapter process handles one attempt, exits, and receives no database credentials. Psi4, RDKit, xTB, Openclatura and OPSIN remain native dependencies.

The frozen pre-migration OpenAPI is `tests/contracts/legacy-openapi.json`. Consumer source snapshots and revision evidence are in `tests/contracts/consumer` and the analysis handoff. `tests/contracts/requests.json` freezes real upstream request shapes and independently calculated descriptor vectors. `scripts/freeze_contract_fixtures.py` is an independent Python encoder; Go tests compare its expected keys and effective options.

## HTTP compatibility

- Existing compute, poll and artifact paths remain. Public fields remain camelCase. Existing snake_case request aliases and known consumer metadata are accepted. Conflicting aliases are rejected. Hints do not control identity or bypass validation.
- Orbital/ESP public states are pending/running/done/error (queued remains a valid legacy state). Conformer projects active states to pending and errors to failed.
- Invalid structural V2000 input is rejected immediately; V3000 is explicitly unsupported. The bounded Go check is structural, not chemical validation. RDKit/Psi4 may still reject a structurally plausible molecular graph. The historical `tests/golden/water.sdf` has a misplaced CTAB version field and is rejected by RDKit; canonical `tests/contracts/water-v2000.mol` is used for native graph qualification.
- Charge and density threshold are retained. Orbital multiplicity other than one is explicitly rejected until open-shell output has its own qualified contract. Psi4 method identifiers keep the legacy bounded identifier syntax. No xTB-to-Gasteiger fallback occurs. Gasteiger reports its condensed implicit-hydrogen approximation.
- Conformer seed/forceField/optimize/numConformers are mapped explicitly; multiple conformers are rejected. Known ETKDG and disabled-optimization aliases normalize to effective options.
- ESP omitted-surface defaults remain 0.25/1.4/3.0. Qualified consumer application `7c8254ab` explicitly sends its intended 0.25/0.9/1.2 values so they agree with its cache key; deploy the paired consumer change with the server. The original patch is a historical proposal, superseded by that integration.
- A saturated payload-read limit returns HTTP 503 `busy` with `Retry-After: 1` and leaves artifact availability unchanged. An unavailable completed artifact returns HTTP 503 `result_unavailable`; execution remains done. Wrong-kind lookup returns 404 without mutation. Readiness returns 503 when the DB, compatible native worker or storage is unavailable. `/live` is process liveness. `/metrics/cache` intentionally describes durable job counts; it does not impersonate the old process-local hit counters.

## Descriptor and adapter protocol 1

Identity hashes the exact full UTF-8 input bytes, operation kind, effective scientific options, fallback policy and per-adapter result-compatibility version. JSON whitespace around the request does not matter; molecular input bytes do. Numeric parameters use exact IEEE-754 binary64 bits, not rounded decimal geometry. The canonical encoder emits typed strings (`s` plus JSON string), floats (`f` plus 16 hex digits), booleans (`b0/b1`), null (`n`), sorted-key maps and ordered arrays, with explicit separators. Orbital sets are sorted and de-duplicated. Canonical JSON strings use Go JSON escaping, including `<`, `>`, `&` and U+2028/U+2029. `waitMs` and non-authoritative hints are excluded.

An attempt request contains protocol, kind, compatibility, input, inputSha256, params, fallback and key. Go writes it into a private `<new-root>/<key>/<random-token>` directory. The adapter is invoked by an argument vector, never a shell, with fixed request/output paths and the attempt token. It writes `response.json`, native artifacts and `manifest.json`. The manifest binds the key, token, input digest, compatibility version, elapsed time, engine/source/build provenance, and each artifact's basename, byte length and SHA-256.

Go kills/reaps the process group before validating publication. It verifies identity, regular files, byte limits, SHA-256, required response fields, compressed-buffer expansion lengths, finite float values, index bounds, vector shape and numerical dipole consistency. It fsyncs and seals files before a conditional PostgreSQL pointer update. A lost acknowledgement is read back using the same token and exact manifest; an unresolved attempt is retained and marked for manual reconciliation.

Compatibility versions initially are `psi4-2`, `esp-2`, `rdkit-conformer-2`, `openclatura-1`. Bump the relevant version by default for scientific code, dependency, compiler, BLAS, native platform or engine changes. Only an independently proven packaging-only rebuild with identical resolved scientific inputs may reuse it. Exact produced image IDs, source hashes and native versions stay in manifests; do not use image rebuild timestamps as compatibility authority. Workers claim only versions they support, and readiness checks that version map. Optional `worker --kinds orbital,electrostatics` or `--kinds conformer,nomenclature` selects a subset; readiness requires the union of live compatible workers to cover all operations. The default worker handles every kind. A missing/stale class makes readiness fail; it does not fabricate capability.

## State and authority

| State/event | Guard | Durable result |
|---|---|---|
| submit | capacity plus one active key, serialized in PostgreSQL | pending or existing compatible job |
| claim | pending row locked with SKIP LOCKED | running, new token, lease, attempt count |
| renew | same token, running, lease still live | extended lease |
| publish | immutable valid bytes, same token and live lease | done plus manifest pointer |
| fail | same token, running, live lease | error |
| lease expiry | running lease expired | pending until 3 attempts, otherwise error |
| artifact loss | completed result integrity failed | availability false; done unchanged |
| late completion/failure | token or lease mismatch | rejected |

Pending age never expires work. Done/error are absorbing. Filesystem bytes alone do not publish a result. A manifest pointer and successful integrity read establish availability. Legacy rows stay in `public.jobs`; Go rows use `moleculens_v2.jobs`, so unfenced legacy broad claims cannot consume the new queue. Historical IDs and roots remain readable.

## Bounds and operations

Default admission is 100 pending/running jobs and 4 HTTP handlers (with at most two integrity-checked payload reads at once). Input is 2 MB; HTTP JSON body is 4 MB. Default native attempt duration is 5 minutes, configurable only between 1 second and 1 hour; naming is 8 seconds. Leases are 30 seconds and renew every third of their duration. A transient renew failure stops/reaps the process conservatively. There are at most three attempts after lease loss. Scientific failures are terminal, not retried as infrastructure failures.

One native attempt runs per worker container. The provided candidate limits each worker to 2 CPUs, 3 GiB and 128 processes; the API uses a 384 MiB Go heap target inside its 512 MiB container and the worker a 256 MiB Go heap target (separate from native chemistry memory). These are qualification settings, not a Hetzner sizing claim. Per-attempt artifacts are capped at 256 MiB, payloads at 32 MiB, logs at 1 MiB; a 512 MiB free-space reserve stops admission/claiming. Use the container lifecycle (or systemd KillMode=control-group) so worker death kills its entire descendant group. The Linux parent-death signal adds protection to the direct child; it is not a substitute for container/cgroup cleanup.

`moleculens gc` lists up to 100 unreferenced retired attempt directories older than 24 hours; `--apply` removes only these after an immediate SQL authority recheck. Every done result, live attempt and ambiguous publication is retained. Completed-result eviction is deliberately not automatic: alert on storage capacity and choose an explicit retention policy rather than silently deleting scientific outputs.

The optional complete `compose.vps1.yml` profile has one long chemistry worker and one naming/conformer worker, queue capacity 12, PostgreSQL and Caddy. Its resource ceilings, measured workload and limits are documented in [vps1-capacity.md](vps1-capacity.md). The durable job/attempt state machine is unchanged; worker classes constrain claim eligibility without creating another queue authority.
