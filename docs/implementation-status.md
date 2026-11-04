# Go migration execution status

**Live follow-up:** the purchased OVH VPS runs the exact qualified Go images and passed native chemistry, recovery, HTTPS, reboot persistence and logical database restore checks. Production DNS, historical data and the paired frontend remain on their existing deployment. The old host's identity is verified through Hetzner's console. See [the launch record](ovh-launch-2026-11-03.md) and [actual-host receipts](qualification/ovh-2026-11-03/deployment.json). Sections below retain the original implementation handoff and local-only evidence.

Owner: commissioned agent thread `moleculens_go_implementation`.
Started: 2026-11-03.
Worktree: `/Users/kelly/Developer/moleculens-server-go`.
Branch: `agent/go-service-migration-20261103`.
Base HEAD: `c9fce3c5d9f69212cc65d39e06ae3c18ef74eafb`.
Preserved overlay commit: `669dc8a` (all 58 source snapshot files verified against the original checkout before indexing the pre-existing changes). Runtime databases/caches and credentials were excluded. The original checkout remains unchanged.

Decision authority: `docs/analysis/2026-11-03-python-to-go/agreement.md` and its final Opus confirmation. Applied the code-lenses and state-machine skills; review and state evidence are in `code-lens-review.md` and `go-protocol.md`.

## Implemented

- Frozen consumer/server snapshots, independently encoded descriptor vectors, explicit option/alias mapping and qualified consumer integration in an isolated companion worktree.
- Bounded Python hardening: active-key uniqueness, attempt-conditional completion/failure, no pending-age expiry, private attempt directories, operation-kind checks, truthful unavailable results and scientific/parser regressions. One explicit checksum-tracked migration owner; no startup schema mutation.
- Go HTTP/admission, full-input descriptor, separate PostgreSQL queue, version-compatible claiming, leases, bounded infrastructure retries, process-group supervision, native deadlines, output/CPU/RAM/storage bounds and immutable hash-validated publication with acknowledgement readback.
- Read-only legacy IDs/artifacts, independent artifact availability, bounded orphan collection, pinned native engine recipe, Go API/worker containers, qualification CI and finite-drain/cutover/rollback runbook. Unsafe automatic SSH deployment is retired in this branch.

## Local qualification completed

- 60 Python unit/regression tests and one real PostgreSQL legacy concurrency/fencing integration test pass.
- Go race-enabled tests pass against isolated UTF8 PostgreSQL 16.14 databases: admission/claim contention, capacity, leases/retries, stale attempts, absorbing terminal states, kind safety, dropped commit acknowledgement with real readback, legacy namespace isolation, unavailable artifacts, cutover/rollback reads, orphan retention and subprocess cancellation/descendant cleanup.
- Native Go HTTP runs pass for Psi4 1.11 HOMO/LUMO/density, RDKit 2026.3.4 Gasteiger ESP and seeded MMFF conformer, xTB 6.7.1 ESP, and Java/OPSIN-verified naming. Repeated submissions reuse the exact job.
- Independent direct Psi4 SCF invocation matches HOMO/LUMO energies exactly with four irreps. Native charged xTB translation conserves +1 charge and preserves the dipole at its declared translated origin within 1e-5 D.
- Initial qualified images at `718696e` pass native HTTP checks, charge +2 and density-only variants. Killing that worker container with SIGKILL reaps it; lease recovery completes the same public job under attempt 2 and a new token. Replacing the candidate API with a separate quiesced compatible API preserves both Go and legacy job reads and identical artifact SHA256 values while rejecting submissions with 503.
- Full 398-case naming benchmark with actual Java OPSIN passes defined gates: 392 named, six guarded, 391 matched and one explicit mismatch (99.7449% match among named). No mismatch is represented as verified.
- Native runs use linux/amd64 emulation on an ARM Mac. Their timings are correctness receipts, not production performance or hosting-cost measurements.

## Current stage

Implementation, clean-source image builds, native qualification and isolated rollback are complete locally. Initial normal-profile application `718696e22f30b658815f093b19c02d7df2567f3a` and exact images remain retained in [qualification/release.json](qualification/release.json). Final application `e2abc7bc814de5acb7ef9992029b29f3990a9235` adds optional worker classes and truthful busy-read errors; its separately qualified clean build is recorded in [qualification/vps1/final-release.json](qualification/vps1/final-release.json). Later handoff commits change documentation/evidence only.

Consumer compatibility integration is complete locally in `/Users/kelly/.openclaw/workspace/projects/moleculens-go-contract`, branch `agent/go-consumer-compatibility-20261103`, application revision `7c8254ab465032098f3ca8d55b76b895e76747fe`, documentation handoff `76c5b419e933bd0ea5d42fbf868c9704b5e59ca7`. Static checks/build, 52 focused tests and seven actual Next-handler-to-Go/PostgreSQL/native-artifact integration tests pass. The seven integration cases also pass against final server `e2abc7b` through Caddy in the full small-host stack. Original consumer source remains unchanged. See [initial consumer receipt](qualification/consumer/receipt.json) and [final paired receipt](qualification/vps1/consumer-pair.json). Paired deployment remains unperformed.

The optional two-CPU/4-GiB complete PostgreSQL/Caddy/API/native-worker trial passed. Fresh 24-atom caffeine Psi4 completed in 67.08 seconds, with Java-verified naming in 5.06 seconds and conformer work in 2.62 seconds while Psi4 remained running; dense ESP, xTB and interrupted 21-atom aspirin recovery also passed. Minimum sampled VM available memory was 2.145 GiB, with zero final-workload OOM events. An earlier unbounded test-harness PostgreSQL OOM and tuning corrections remain explicitly documented. Current plus retained rollback image storage is approximately 5.30 GB. This supports a low-traffic small-molecule deployment trial, not unrestricted chemistry or actual OVH performance. See [capacity scope and measurements](vps1-capacity.md).

See [the evidence index](qualification/README.md), [state and adapter contract](go-protocol.md), [review](code-lens-review.md) and [release runbook](go-release-runbook.md). The prepared branch is reviewable and locally runnable. No push, merge, live deployment, provider purchase, DNS/billing change or remote CI run has occurred.

## Remaining live-release gates

Trusted verification of the changed SSH host identity; actual deployed topology/roles/resources and active-job inventory; paired consumer deployment; production backup restore; finite legacy process drain; host-specific ingress/rollback rehearsal and representative capacity measurements. These cannot be inferred from local source or tests. See `go-release-runbook.md`.

No Brev resources were used. All five small-host task containers and their network were removed, the disposable host PostgreSQL cluster is stopped, and Colima profile `moleculens-qualification` is stopped. No listeners remain on qualification ports 55439, 55440, 58082, 58083 or 58443. Exact images and isolated test volumes remain in the stopped VM for review; see [final cleanup](qualification/vps1/cleanup.json).
