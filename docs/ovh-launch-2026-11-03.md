# OVH production launch — November 3, 2026

## Live state

**https://api.moleculens.com** now serves the Go backend from OVH **135.148.205.154**. **https://www.moleculens.com** serves the compatible frontend. Cloudflare's existing `api` A record remains DNS-only with a 300-second TTL; both authoritative nameservers and two public resolvers return OVH. No API AAAA record exists. Trusted HTTPS, readiness `go-2`, CORS preflight (204), public metrics blocking (404), actual frontend Go-job reads and historical conformer/ESP reads pass. The browser renders the production 3D viewer. See [cutover receipt](qualification/ovh-2026-11-03/production-cutover.json).

The old Hetzner API and worker are stopped with restart disabled. Its Caddy forwards remaining old-IP clients to the fixed OVH address over verified TLS, with the provider hostname as TLS server name. Its source database is stable and every historical file matches the original snapshot across a 30-second observation interval. Twelve `conformer_log_emitted` flags added by old readback were synchronized; all scientific parameters and other job fields stayed unchanged. Permanent resource retirement is still pending action-time confirmation. Private account/billing details remain outside source control.

## Exact runtime

- OVH order 9375564, host `vps-12a3a763.vps.ovh.us`, Virginia `os-us-east-va-2`; IPv6 `2604:2dc0:101:200::47f0`.
- Ubuntu 26.04.1, kernel `7.0.0-14-generic`, amd64, two vCores, approximately 4.00 billion bytes kernel-visible RAM, no swap, 40 GB disk. Approximately 20 GB is free after restored data and qualification cleanup.
- Docker 29.9.3 and Compose 5.7.0 from Docker's official repository. SSH is key-only as `ubuntu`; root/password SSH is disabled. Inbound TCP 22/80/443 is allowed for both address families. Only Caddy publishes application ports.
- Runtime `/opt/moleculens`, project `moleculens-vps1`, root-private `.env` (0600). Legacy storage is separate and read-only to the API; new storage is owned by UID/GID 65532.

The API and both workers use the clean public build **deca10b457436e428c17ccb31bc9269aaddf98dc**:

| Component | Image |
|---|---|
| API | `sha256:cabeb6d6e30617fde35c45be5d7deec9207bec1195994d447783bb1991fada7f` |
| Worker | `sha256:18bfc72609c3edf4232843d30066bf5c03432b764ec8301f393fb2bbfbbcc3f8` |
| Unchanged native engine | `sha256:57b7901dbaee58a29fa98c48b54f9c62267d68bd0bd60a89b680a40bc4293500` |
| PostgreSQL 16.14 | `postgres:16.14-alpine@sha256:57c72fd2a128e416c7fcc499958864df5301e940bca0a56f58fddf30ffc07777` |
| Caddy 2.10.2 | `caddy:2.10.2-alpine@sha256:4c6e91c6ed0e2fa03efd5b44747b625fec79bc9cd06ac5235a779726618e530d` |

After the bounded large-conformer replay ended, the conformer/naming worker was replaced with the final worker image. A fresh optimized conformer completed in 3.27 seconds, cache reuse passed, and OPSIN-verified naming completed in 4.55 seconds. All five services run, readiness is healthy, no requests remain active, and both operations timers are enabled. See [final acceptance](qualification/ovh-2026-11-03/final-acceptance.json).

The final image archive is retained off-host: 1,212,602,880 bytes, SHA256 `b0e4eb3c2103f06777dcebe4d9ec1e2464187aaf7d38c5c88606ef78a9d5b514`. It was verified after transfer and loaded by immutable IDs; redundant server transfer archives were removed. The initial qualified Go API `sha256:a54c9f890d208e0a0a4b06f7b1c1096847d695ecc32250bac38930ec96f62574` remains available for quiesced rollback.

## Data and recovery evidence

All **448 historical jobs** (439 completed, nine failed), **439 cache entries**, and **5,414 files** (11,548,152,870 bytes) are restored. Every historical HTTP response preserves previously exposed scientific values, IDs, timestamps and timings; extra stored fields and omitted optional nulls are permitted. All restored historical file SHA256 values match the source. Cube/JSON downloads preserve bytes, MIME types and attachment filenames. See [HTTP comparison](qualification/ovh-2026-11-03/legacy-http-comparison.json), [artifact restore](qualification/ovh-2026-11-03/legacy-artifact-restore.json) and [downloads](qualification/ovh-2026-11-03/artifact-downloads.json).

A final logical dump was restored into a separate disposable database. Historical tables match, all terminal Go snapshot rows and immutable persisted inputs match, and referenced completed Go files match restored manifests. The snapshot retained one then-active qualification request, which subsequently reached its deadline. Its private working directory was excluded from the artifact archive because a persisted active request resumes with fresh attempt ownership after recovery. All 5,414 historical files were rechecked. The new Go artifact archive was decompressed/read back on the server and Mac, verifying 201 files. Configuration/TLS archives and the database dump also match off-host hashes. Private backups remain under the original workspace's ignored `build/legacy-backup/` and `build/ovh-live/final-candidate/`. See [production restore receipt](qualification/ovh-2026-11-03/production-backup-restore.json).

Daily root-only database backups retain seven dumps under `/var/backups/moleculens`. Provider Standard automated backup is activated; a full OVH snapshot restore is untested. Conservative artifact GC retains completed results, active attempts and uncertain publications, collecting only retired unreferenced attempts older than 24 hours, at most 100 per run. Its operations service mounts only new artifact storage. Units and the verifier are in [ops/ovh](../ops/ovh/README.md).

An actual HTTPS ingress replacement with the retained quiesced Go API preserved Go and legacy responses byte for byte and returned 503 for new admission, then restored normal routing. See [rollback rehearsal](qualification/ovh-2026-11-03/rollback-ingress.json). Rollback retains the Go queue and both roots; never restart the old unfenced Python worker against the migrated database.

## Capacity and limits

One chemistry worker has 1 CPU/1,792 MiB; one conformer/naming worker has 1 CPU/768 MiB. API is 512 MiB, PostgreSQL 256 MiB and Caddy 96 MiB. Chemistry/private scratch is capped at 8 GiB; conformer/naming at 1 GiB. A worker requires its entire budget plus 2 GiB free before claiming, while API readiness/admission requires 2 GiB. The queue permits 12 pending/running requests. Published artifacts stay capped at 256 MiB, response payloads at 32 MiB, outer attempts at two hours and naming execution at eight seconds.

- Native small Psi4, RDKit ESP, xTB, conformer, Java/OPSIN and cache reuse checks passed. Fresh chemistry under the final 8-GiB setting also passed, concurrently with the large conformer.
- The largest historical density input (74 atoms, B3LYP/6-31g(d)) passed in 2,030.91 seconds. Maximum sampled attempt storage was 4,495,376,302 bytes. The run used a 16-GiB allowance; 8 GiB was selected from the measured peak, with regression coverage for that measured size. Minimum sampled available RAM was 1,775,759,360 bytes. Chemistry reached its cgroup ceiling, with zero recorded OOM events. [Density replay](qualification/ovh-2026-11-03/historical-density.json), [telemetry](qualification/ovh-2026-11-03/historical-density-memory.json).
- The initial 256-MiB, 1-GiB and 4-GiB scratch allowances were inadequate. Only regenerable scratch from three failed migration fixtures was reclaimed after checking terminal rows, no manifests/uncertain publication/live reference/open writer; original inputs/logs/hashes remain private. Historical user artifacts were untouched.
- The largest 785-atom conformer replay failed at the two-hour deadline with `adapter process: context deadline exceeded`; its earlier one-hour replay also failed. No OOM events were recorded during the two-hour observation. Very large fresh conformers remain unqualified on this setup. The historical random seed was unrecorded, so its previous runtime is context rather than an identical stochastic trajectory. All saved historical conformers remain readable without recomputation. The final acceptance receipt records this failure separately from passing operational checks.
- Earlier actual-OVH mixed work completed caffeine B3LYP in 87.39 seconds, concurrent naming in 5.26 seconds and a conformer in 3.33 seconds. Interrupted-worker recovery preserved the same public ID under attempt 2. An earlier real host reboot retained 11 Go rows and 29 artifact hashes and restarted all five services. The host has not been rebooted again after the later data import and image replacements; logical restore and rolling container replacements provide the later persistence evidence.

This is evidence for the measured low-traffic workload, not every syntactically accepted molecule or sustained concurrent traffic. Giant conformers can queue other work in their class and may reach the finite deadline. The 999-atom structural limit is not a performance guarantee. Do not add worker replicas without revisiting RAM, CPU and storage reserves.

## Source and operation

Backend [PR 30](https://github.com/jakekinchen/moleculens-server/pull/30) merged as `27931a359543dea2a409a2641619df7846ba7722`. The paired frontend [PR 127](https://github.com/jakekinchen/moleculens/pull/127) merged as `a56e1a23fc36391dfc94b74703c9e2ad443b9c8b`; it contains only compatibility work on the actual production baseline, preserving unrelated unpublished UI work. The frontend deployment from merged main is Ready, and actual new Go and historical conformer reads passed again after that deployment. Production `MOLECULENS_API_URL` is explicitly `https://api.moleculens.com`.

The actual PostgreSQL/native CI qualification and secret scan pass for the released application source; local Go race tests/vet and whole-tree Python lint/format checks pass. Authentication/premium guards in the earlier Next-handler integration tests were mocked; public production reads and the browser check are separate evidence, not a signed-in premium-flow test.

Use the original workspace's private `build/ovh-live/ssh-config` with host `moleculens-ovh` for administration. Keep the old host's pinned trust separate. See [the release runbook](go-release-runbook.md) and [operations instructions](../ops/ovh/README.md). The disposable local Colima builder is stopped. No Brev resource was used.
