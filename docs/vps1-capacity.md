# Optional 2-CPU/4-GiB backend capacity result

**Actual-host follow-up:** the user purchased the Virginia VPS-1. Native amd64 qualification on Ubuntu 26.04 passed, with approximately 2.11 GiB minimum sampled available RAM, responsive mixed work and no recorded OOM in that run. Caffeine Psi4 completed in 87.39 seconds. HTTPS and reboot persistence passed. See [the live launch record](ovh-launch-2026-11-03.md) for exact scope; production data fit and sustained traffic remain unproven. Measurements below describe the earlier local ARM/Rosetta trial.

The complete local backend passes the measured small-molecule workload within 4 GiB and two CPUs. This supports a **conditional go for a low-traffic deployment trial** of the proposed OVH VPS-1 using the optional `compose.vps1.yml` profile. It does not establish that every accepted molecule or an unknown production backlog fits. No hosting was purchased, moved or deployed.

The candidate supplied by the user is OVH VPS-1 2027, US-WEST-OR, two vCores, 4 GB RAM, 40 GB NVMe and Ubuntu 26.04. The actual test was a disposable Colima VM running Ubuntu 24.04.4 on ARM, with linux/amd64 containers under Rosetta, two CPUs, 4,095,635,456 bytes of kernel-visible RAM and no swap. The actual OVH CPU, Ubuntu 26.04, network, HTTPS certificate issuance and production data remain unqualified. In particular, verify the purchased VM's kernel-visible RAM rather than treating the label “4 GB” as an exact byte guarantee.

## Exact revisions and isolation

- Server application: `e2abc7bc814de5acb7ef9992029b29f3990a9235`; clean build and exact image/binary identities in [final-release.json](qualification/vps1/final-release.json).
- Consumer application: `7c8254ab465032098f3ca8d55b76b895e76747fe`; documentation handoff: `76c5b419e933bd0ea5d42fbf868c9704b5e59ca7`. All seven actual Next-handler integration tests passed through Caddy against the final server; [paired receipt](qualification/vps1/consumer-pair.json).
- Native engine: unchanged pinned Psi4/RDKit/xTB/Openclatura/Java image `sha256:57b7901dbaee58a29fa98c48b54f9c62267d68bd0bd60a89b680a40bc4293500`.
- Initial normal-profile application `718696e22f30b658815f093b19c02d7df2567f3a` and its qualified images remain retained. No original server or consumer source was edited; preservation receipts cover 58 and 24 files respectively.

The optional profile runs PostgreSQL, Caddy, Go API, one long chemistry worker and one short-operation worker. Worker kind selection separates orbitals/ESP from naming/conformers using the same durable queue and attempt protocol. Readiness requires the union of live compatible workers to cover every operation. The default worker still handles all kinds. Normal-profile limits remain unchanged; its tmpfs YAML was corrected to keep the comma-containing mount option in one string.

## Measured workload

Inputs have fresh titles, so the principal mixed-load cases compute new results rather than reuse an existing cache entry. See [final-workload.json](qualification/vps1/final-workload.json) and `scripts/qualify_small_host.py`.

| Case | Observed result |
|---|---|
| Caffeine, 24 atoms, B3LYP/6-31g(d), HOMO/LUMO and density at 0.25 Å | Completed in 67.08 s; total energy -680.3712676247128 hartree. |
| Fresh benzoic-acid naming while caffeine Psi4 remained running | Java/OPSIN matched in 5.06 s using a distinct worker owner. |
| Caffeine MMFF conformer during the same Psi4 calculation | Optimized in 2.62 s. |
| Concurrent status polling during mixed work | 637 background polls; maximum 0.208 s, mean 0.00357 s; no busy responses in this final run. |
| Dense caffeine Gasteiger ESP, 0.1 Å grid | Completed in 2.95 s with 44,622 surface vertices. |
| Caffeine xTB 6.7.1 ESP | Completed in 2.36 s with 7,320 surface vertices. |
| Aspirin, 21 atoms, B3LYP/6-31g(d), interrupted during execution | SIGKILL reaped the worker container; after manual container restart, automatic lease recovery completed the same public job under attempt 2 with a new token/owner. 81.55 s including lease wait; final computation 50.68 s. |

The basic native HTTP qualifier and seven consumer integration cases also passed against the final API. Basic water artifacts were populated earlier with the identical native engine and reused; the larger cases above are fresh final-worker executions. Earlier independent neutral/charged Psi4 comparisons, charged xTB translation and the 398-case naming benchmark remain separate scientific receipts, not new small-host throughput claims.

The actual PostgreSQL Go race suite passed with serial test packages and a bounded direct test-helper pool. It includes worker-class claims/readiness and a deterministic result-read saturation test. At capacity, result reads now return HTTP 503 `busy` with `Retry-After: 1`; saturation does not falsely mark a completed artifact unavailable. Scientific/result compatibility versions did not change.

## Measured memory and storage

[Final memory telemetry](qualification/vps1/final-memory.json) records kernel cgroup peaks, plus VM and disk samples at 0.5-second intervals. Per-container peaks need not occur simultaneously. The PostgreSQL peak includes the controlled regression suite; VM minima and scratch maxima are sampled, not exact instantaneous high-water marks.

| Component | Configured memory ceiling | Observed kernel cgroup peak |
|---|---:|---:|
| Long chemistry worker | 1,792 MiB | 1,156.95 MiB |
| Naming/conformer worker | 768 MiB | 280.94 MiB |
| Go API | 512 MiB | 59.42 MiB |
| PostgreSQL | 256 MiB | 181.52 MiB |
| Caddy | 96 MiB | 27.71 MiB |

Maximum sampled simultaneous container memory was **1.360 GiB**. Minimum VM `MemAvailable` was **2.145 GiB**. No final-container OOM events occurred, and no swap was configured or used. Psi4's 1,024 MB setting is an engine allocation setting, not total worker RSS; the cgroup peak above is the relevant process-tree evidence. Configured ceilings total 3,424 MiB and do not substitute for these measurements or guarantee that simultaneous worst-case allocations are safe.

Logical artifact-volume high-water was 56,692,416 bytes; PostgreSQL-volume high-water was 73,059,588 bytes; largest sampled attempt directory was 42,301,489 bytes. These figures cover this fixture corpus, not indefinite retention or historical production data.

Docker's stored-layer accounting gives a conservative **5,252,992,377-byte** upper bound for the current four images, plus **50,817,506 bytes** of unique layers for the retained initial API/worker rollback images: **5,303,809,883 bytes total**, about 5.30 GB or 4.94 GiB. This includes shared native layers once. `docker image inspect` content/blob sizes are smaller and are not installed disk usage. [disk.json](qualification/vps1/disk.json) preserves the accounting. The test VM's Docker partition used 7.12 GB including additional development images; its root used 1.25 GB. Its provisioned virtual disk sizes are not workload footprint.

Forty GB is sufficient for the software plus these measured fixtures with substantial space remaining. Production database/artifact size, backups, build caches, logs, old releases and retained outputs must be inventoried before claiming the whole historical deployment fits. Build images elsewhere and transfer immutable releases. The 512 MiB free-space admission guard prevents unsafe new work but is not a storage-retention plan.

## Tuning failures and limits

The first test-harness run exceeded the tiny database's connection/memory budget: parallel Go packages and direct pgx pools sized from the Mac bypassed HTTP admission bounds. PostgreSQL hit its **256 MiB cgroup OOM limit**; [kernel evidence](qualification/vps1/initial-test-harness-oom.txt) and initial test logs are retained. Rerunning with `-p 1` and `pool_max_conns=4` bounded the test helper's directly constructed pool and passed. `Store.Open` itself explicitly caps service pools at 12; the URL parameter does not override that service cap. Final application traffic used those service pools, four HTTP handlers and one active attempt per worker. Do not run unrestricted administrative/test concurrency against this small live database or add worker replicas without requalification.

A 0.5-CPU short worker missed the eight-second native naming deadline. Its final ceiling is one CPU, alongside one CPU for long chemistry; two historical naming failures remain explicit in the database inventory. Caddy's initial 64 MiB ceiling was too close to its observed peak and was raised to 96 MiB. Concurrent result reads exposed the misleading `result_unavailable` classification, corrected and deterministically regression-tested before the final build. Final application cases passed after these changes; “no final-workload OOM” must not be read as “no OOM occurred during tuning.”

The eight-second naming deadline, five-minute long-job deadline and queue capacity 12 remain finite. A short-job class prevents a long Psi4 job from monopolizing naming execution, but does not guarantee zero waiting behind other short jobs. The structural validator accepts up to 999 atoms; that syntactic limit is **not a promise that 999-atom chemistry completes on this host**. Large basis sets, dense meshes on larger molecules, sustained traffic, multiple users and production-sized retention have not been qualified. Maintain one long and one short worker until representative workload evidence supports changes.

## Optional profile operation and reproduction

Use the immutable API and worker identities in `final-release.json`, the pinned PostgreSQL/Caddy images in `compose.vps1.yml`, and a private environment file with a new database password and encoded `DATABASE_URL`. Set `API_HOST`, `ALLOWED_ORIGINS` and ingress bindings for the intended environment. The profile defaults to loopback bindings. The local trial used `API_HOST=:80`, ports 58083/58443 and HTTP only. A public trial would require an explicitly authorized hostname, ports 80/443, public binding, firewall and TLS checks; none was performed here. `/metrics/*` is blocked at Caddy.

Initialize only new isolated storage, migrate explicitly once, then start services:

```sh
docker compose --env-file <private-env> -f compose.vps1.yml run --rm init-storage
docker compose --env-file <private-env> -f compose.vps1.yml up -d postgres
docker compose --env-file <private-env> -f compose.vps1.yml run --rm migrate
docker compose --env-file <private-env> -f compose.vps1.yml up -d api chemistry short caddy
```

The `legacy` volume is read-only to the API. Existing production data must be restored/imported deliberately after backup and verified finite drain; the empty test volume is not a migration of historical results. Follow [the release runbook](go-release-runbook.md) for host trust, backup/restore, cutover and quiesced compatible rollback.

For repeated local proof, run `scripts/qualify_go_service.py` to populate canonical native cases, `scripts/qualify_small_host.py` for fresh mixed-load/interruption cases, and `scripts/measure_vps1_vm.py` as root inside the disposable VM for actual cgroup/VM evidence. The latter takes an output path and fresh stop-marker path; its sampling must span the workload. The small-host script takes API URL, output and Docker context. Run consumer integration with the explicit local URL as described in its verification document. Use `TEST_DATABASE_URL='<isolated DSN with pool_max_conns=4>' go test -race -p 1 -count=1 ./...` for bounded local regression checks. Test fixtures may create/drop only isolated test databases.

Stop containers and the disposable VM when proof is complete. Exact images and isolated volumes remain retained in the stopped local VM for reproducibility; [cleanup.json](qualification/vps1/cleanup.json) records the final observed state. No Brev resource was used. Approval to purchase or deploy this host is outside this qualification.
