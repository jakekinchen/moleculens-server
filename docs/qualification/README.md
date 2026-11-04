# Qualification evidence

Initial normal-profile application revision: `718696e22f30b658815f093b19c02d7df2567f3a`, built from a clean worktree. The root-level receipts below refer to that revision; `release.json` records its exact binary and image identities. Final application `e2abc7bc814de5acb7ef9992029b29f3990a9235` adds optional worker classes and truthful busy-read error classification; it has its own clean-build and full small-host receipts under `vps1/`. The native engine image and scientific compatibility are unchanged. Commits after `e2abc7b` record documentation/evidence only. Timings use amd64 emulation on an ARM Mac and do not measure Hetzner or OVH performance/cost.

| Receipt | Proven scope |
|---|---|
| `release.json`, `final-container-policy.json` | Clean build identity; actual final API/worker image IDs, non-root user, read-only filesystem, memory limits, native CPU/PID limit and mount permissions. |
| `go-tests.txt` | Uncached race-enabled Go suite with real isolated UTF8 PostgreSQL 16.14: contention, capacity, fencing, lease/retry, dropped commit acknowledgement/readback, subprocess cleanup, artifact integrity/availability, legacy routing and GC retention. |
| `python-tests.txt` | 60 unit/regression tests and one real PostgreSQL legacy concurrency/fencing test pass. |
| `final-native-http.json` | Exact final Go API/worker images complete native Psi4 orbitals/density, RDKit Gasteiger ESP and MMFF conformer, xTB ESP, actual Java/OPSIN naming, and exact job-ID cache reuse. |
| `final-charge-density.json` | Final images honor density-only selection/threshold and charge +2, producing expected distinct meshes/energies. |
| `psi4-independent-reference.json`, `psi4-charged-reference.json` | Independent direct native SCF values match served neutral and charged HOMO/LUMO/total energies; four-irreps fixture exercises orbital energy ordering. |
| `xtb-charged-translation.json` | Native +1 ammonium charge conservation and dipole translation at the declared origin within 1e-5 D. |
| `final-crash-recovery.json` | Final worker container SIGKILL/restart: stopped PID 0/exit 137, same public job ID, attempt 1 to 2, new token and matching publication authority, successful native result. |
| `final-rollback.json` | Actual local candidate shutdown/replacement with a separate quiesced compatible API; both Go and seeded legacy job IDs and identical artifact hashes retained; admission rejected with 503. |
| `nomenclature-opsin.json` | Full 398-case corpus using local RDKit/Openclatura and Java OPSIN: 392 named, six guarded, 391 matched, one disclosed mismatch; defined gates pass. Linux native smoke proof is separate in `final-native-http.json`. |
| `consumer-compatibility.patch` | Historical prepared/apply-checked patch, superseded by the actual qualified consumer branch below. It is not the final implementation. |
| `original-preservation.json` | All 58 original source snapshot hashes and original HEAD remain unchanged after implementation. |

Earlier `native-*`, `packaged-*` and intermediate receipts are retained as development history. The `final-*` receipts are the final image qualification. The original analysis agreement/validation files preserve the pre-implementation state and are not rewritten as current implementation claims.

## Actual consumer and optional full small-host qualification

Consumer application `7c8254ab465032098f3ca8d55b76b895e76747fe` is implemented and qualified in isolated branch `agent/go-consumer-compatibility-20261103`, with a docs-only handoff at `76c5b419e933bd0ea5d42fbf868c9704b5e59ca7`. Its original checkout remains unchanged. The final server pairing supersedes the limitation of an unapplied patch.

| Receipt | Proven scope |
|---|---|
| `consumer/receipt.json`, `consumer/check.txt`, `consumer/build.txt`, `consumer/targeted-tests.txt`, `consumer/live-go.txt` | Consumer static checks, build, 52 focused tests and seven actual Next-handler/HTTP integration tests against the initial qualified server. |
| `vps1/consumer-pair.json`, `vps1/consumer-integration.txt` | The same seven actual integration tests pass through Caddy against final server `e2abc7b`; only premium/access dependencies are stubbed. |
| `vps1/final-release.json`, `vps1/final-containers.json` | Clean final application build and exact API/worker/PG/Caddy identities, actual limits and container policy. |
| `vps1/go-tests.txt` | Uncached Go race tests on actual PostgreSQL, serial packages and bounded direct helper pool; worker-class/readiness and result-read saturation regressions included. |
| `vps1/final-workload.json` | Fresh 24-atom caffeine Psi4 plus concurrent Java naming and conformer, dense ESP, xTB, status responsiveness, and 21-atom aspirin worker-SIGKILL/restart/lease recovery. |
| `vps1/final-native-http.json` | Final API reads the canonical native cases through complete stack; the same-engine basic artifacts were retained from earlier trial population. |
| `vps1/final-memory.json`, `vps1/disk.json` | Actual two-CPU/no-swap VM cgroup peaks and sampled VM headroom/attempt storage; Docker installed/shared-layer accounting for current and retained rollback images. |
| `vps1/initial-test-harness-oom.txt`, `vps1/initial-unbounded-test-pools.txt`, `vps1/initial-postgres.log` | Preserved early PostgreSQL cgroup OOM from uncontrolled test-helper connections; not omitted from the final result. |
| `vps1/initial-short-cpu-timeout.txt`, `vps1/initial-result-busy-mislabel.*` | Preserved tuning failures corrected before the final workload: short-worker CPU budget and transient read error classification. |
| `vps1/original-preservation.json` | Rehashed all 58 original server and 24 original consumer source files plus unchanged original HEADs after all work. |
| `vps1/final-job-inventory.json`, `vps1/cleanup.json` | No unfinished task jobs; five containers/network removed, host PostgreSQL stopped, VM stopped and qualification ports closed. |

See [the measured capacity scope](../vps1-capacity.md). `vps1/initial-*`, `images.json`, `workload.json` and other non-final files preserve intermediate runs; use `final-*` plus the final disk/consumer receipts for the conclusion. The original root `cleanup.json` records initial cleanup; `vps1/cleanup.json` records the final cleanup after resuming the VM.

The implementation owner ran these checks. This is not an independent final advisor sign-off. Remote CI, complete production topology, trusted replacement SSH identity, live consumer deployment, production backup restore, ingress/TLS cutover and representative measurements on the actual host remain release gates. Local evidence supports the measured small-molecule envelope, not every structurally accepted input or unknown production data size. No production or paid hosting changes were made.
