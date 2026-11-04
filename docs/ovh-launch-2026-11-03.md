# OVH VPS launch — November 3, 2026

## Current state

The purchased OVH VPS is running a Go candidate at **https://vps-12a3a763.vps.ovh.us**. Production `api.moleculens.com` still resolves to Hetzner `5.161.53.230`; no production traffic switch or Hetzner cancellation has occurred.

Historical migration is staged: all 448 public job rows and 439 cache rows match the source; all 5,414 historical artifact file hashes match. Every historical API response has been compared, preserving previously exposed scientific values (additive fields and omitted optional nulls are permitted). See [HTTP comparison](qualification/ovh-2026-11-03/legacy-http-comparison.json). An off-host archive and isolated database/artifact restore are verified. The frontend compatibility change is based on the current production revision, pushed as [draft PR 127](https://github.com/jakekinchen/moleculens/pull/127), and built as an unpromoted production deployment.

Capacity remains a gate. The largest 74-atom density replay exceeded both 1-GiB and 4-GiB private scratch limits; the 4-GiB attempt allocated 4,377,811,800 bytes before termination. Inspection of the same pinned Psi4 1.11 engine found 694 primary and 3,676 auxiliary basis functions (three-index tensor estimates: 14,163,951,488 bytes full; 7,092,180,320 bytes triangular). These are sizing estimates, not measured successful peak usage. The largest density replay subsequently passed in 2,030.91 seconds under a 16-GiB allowance. Maximum sampled attempt storage was 4,495,376,302 bytes; minimum sampled available VM memory was 1,775,759,360 bytes. The chemistry cgroup reached its 1,792-MiB ceiling, with zero recorded OOM events. See [density replay](qualification/ovh-2026-11-03/historical-density.json) and [memory/storage](qualification/ovh-2026-11-03/historical-density-memory.json). The production profile now gives chemistry 8 GiB and conformer/naming 1 GiB, and requires each worker's budget plus a 2-GiB reserve before claiming. The large replay used the earlier 16-GiB cap; the final 8-GiB cap was selected from the measured peak, with a sparse-file regression covering that measured size. Fresh small native cases are checked under the final cap. API readiness/admission uses the 2-GiB reserve, so an owned large workspace does not itself disable small operations. Workers wait for space; referenced results are never evicted to admit work. Published artifacts remain capped at 256 MiB and responses at 32 MiB. Native scientific parameters, dependencies, CPU and RAM limits remain unchanged across these scratch/runtime tuning builds. Outer deadlines become two hours to accommodate historical conformers recorded above one hour; naming retains its eight-second deadline. The 24-hour orphan retention remains longer than the maximum process deadline.

Only scratch from three terminal failed migration fixtures was reclaimed (6,588,078,936 bytes), after verifying no manifest, uncertain publication, live row reference or open writer; request/log/error metadata and scratch hashes are retained privately. Original historical artifacts remain untouched. The preserved one-hour 785-atom conformer reached its execution deadline. A fresh replay under the two-hour profile is running. This uses the new contract's deterministic seed 0; the historical deployment did not record its random seed, so its previous runtime is context rather than an identical stochastic trajectory. Historical geometries remain preserved byte for byte. Large-conformer qualification remains pending.

Hetzner Cloud contains one production VM, two primary IPs and one firewall; the other Cloud project and Robot server inventory are empty. Account closure is separate from resource retirement and requires final billing settlement. Private billing records and the account-specific closure status remain outside source control. Paid resource retirement can stop future usage after migration acceptance; permanent deletion requires action-time confirmation. The actual OVH HTTPS rollback/ingress rehearsal passed with byte-identical Go and legacy job responses and HTTP 503 for quiesced admission; normal ingress was restored. See [rollback receipt](qualification/ovh-2026-11-03/rollback-ingress.json). Daily backup units and a full logical restore verifier are retained under `ops/ovh/`. A live cube download comparison found identical bytes but missing legacy MIME/attachment headers; the API now preserves these download headers, with focused HTTP regression coverage.

- OVH order `9375564`; host `vps-12a3a763.vps.ovh.us`; IPv4 `135.148.205.154`; IPv6 `2604:2dc0:101:200::47f0`.
- Virginia `os-us-east-va-2`, Ubuntu 26.04, Linux amd64, two vCores, 4,000,071,680 bytes kernel-visible RAM, no swap, 40 GB disk. Approximately 30.5 GB remained after qualification, including the transfer archive.
- Standard provider automated backup is activated; complete snapshot restore is untested.
- The user approved reinstalling this new VPS with this Mac's existing Moleculens SSH key. Key-only administration as `ubuntu` and passwordless sudo were verified before disabling password authentication and root SSH. Access survived reboot.
- Docker Engine 29.9.3 and Compose v5.7.0 came from Docker's official Ubuntu repository. OS updates and reboot completed; kernel `7.0.0-14-generic`.
- Host incoming policy allows TCP 22/80/443 for both address families. Only Caddy publishes application ports; PostgreSQL and the API remain on the private Compose network. Caddy blocks `/metrics/*`.
- TLS 1.3 and a trusted Let's Encrypt certificate were verified externally. `/ready` reports healthy `go-2`; `/live` returns 200; `/metrics/cache` returns 404.

## Initial native qualification deployment

Clean application: `e2abc7bc814de5acb7ef9992029b29f3990a9235`. The deployed `compose.vps1.yml` profile runs PostgreSQL, Go API, one long chemistry worker, one short-operation worker and Caddy. A separate `compose.public.yml` override publishes ports. Runtime files are under `/opt/moleculens`; credentials remain private and excluded from source.

| Component | Exact image |
|---|---|
| API | `sha256:a54c9f890d208e0a0a4b06f7b1c1096847d695ecc32250bac38930ec96f62574` |
| Worker | `sha256:8a5608c6b4d068e1970e621c644eec069c1cab4726dfbec69b6cb0e458dc7e14` |
| PostgreSQL 16.14 Alpine | `sha256:57c72fd2a128e416c7fcc499958864df5301e940bca0a56f58fddf30ffc07777` |
| Caddy 2.10.2 Alpine | `sha256:4c6e91c6ed0e2fa03efd5b44747b625fec79bc9cd06ac5235a779726618e530d` |

The linux/amd64 image archive is 1,352,139,776 bytes, SHA256 `6b5476fe76745fb4a7bfd9a43edda3ca031aa712d30e51a4d7f53b93fca88044`. All 66 blob hashes and both qualified application manifests/configs/layers were checked; the transferred archive matched. Images were built elsewhere. See [deployment.json](qualification/ovh-2026-11-03/deployment.json).

## Actual OVH qualification

These are native Linux amd64 measurements from the purchased machine, separate from the earlier Mac/Rosetta rehearsal.

- Native HTTP fixtures passed for Psi4 orbitals/density, RDKit Gasteiger electrostatics, xTB electrostatics, MMFF conformers, Java/OPSIN naming and exact cache reuse. Water orbital energies agreed with the earlier reference within approximately 2.4e-12 eV.
- Fresh 24-atom caffeine B3LYP/6-31g(d) completed in 87.39 seconds. Concurrent naming took 5.26 seconds and conformer work 3.33 seconds while long chemistry remained busy. Dense ESP and xTB took 5.19 and 3.09 seconds.
- 851 background polls had maximum latency 0.0723 seconds and mean 0.00503 seconds. This is fixture evidence, not a sustained traffic benchmark.
- SIGKILL terminated the long worker; after manual container restart, automatic fenced lease recovery completed the same aspirin job under attempt 2, in 99.87 seconds including lease wait.
- Minimum sampled available RAM was 2.11 GiB. Original cgroup memory peaks were approximately 1,428 MiB long worker, 411 MiB short worker, 17 MiB API, 79 MiB PostgreSQL and 20 MiB Caddy. The sampler followed original cgroups and may miss peaks after worker replacement. Final container OOM flags were also checked; no recorded OOM occurred in this run.
- Reboot changed the boot ID and all five services restarted. All 11 job rows were byte-identical and 29 artifact file hashes were unchanged. Public native checks reused the same IDs after reboot; no failed systemd unit remained.

Receipts: [native checks](qualification/ovh-2026-11-03/native-http.json), [mixed workload](qualification/ovh-2026-11-03/mixed-workload.json), [memory](qualification/ovh-2026-11-03/native-memory.json), [HTTPS](qualification/ovh-2026-11-03/https.json), [reboot persistence](qualification/ovh-2026-11-03/after-reboot.json).

## Historical workload limits found during migration

The production inventory contains a 74-atom B3LYP/6-31g(d) density job that took approximately 20 minutes on Hetzner, and 785-atom conformer jobs, including one recorded at approximately 61 minutes. Replaying the largest density job on the initial OVH profile failed after 105.67 seconds because its private Psi4 scratch file grew to approximately 698 MB, exceeding the 256 MiB attempt budget. The largest conformer replay reached the initial 45-second deadline. Neither failure recorded an OOM. Receipts: [initial historical replay](qualification/ovh-2026-11-03/historical-initial.json), [memory](qualification/ovh-2026-11-03/historical-initial-memory.json).

The executor separates temporary attempts from the unchanged 256 MiB published-artifact and 32 MiB response limits. A first 1 GiB scratch allowance still stopped the largest orbital job when its scratch file reached approximately 1.51 GB. The 4-GiB attempt also failed; the current-state section above supersedes this intermediate allowance. The optional VPS profile permits one-hour outer attempt deadlines for the historical long tail; naming retains its separate eight-second execution deadline. Memory/CPU/concurrency and scientific parameters remain unchanged. A sparse-file regression failed before separating the budgets and passed afterward; engine, contract and HTTP race tests plus focused vet passed. Large conformers can occupy their worker and delay other short-class work; the initial caffeine/naming responsiveness result does not establish zero queueing behind those larger conformers.

## Initial backup and recovery scope

Root-only daily `moleculens-database-backup.timer` is enabled. The script takes a custom-format PostgreSQL dump under `/var/backups/moleculens`, validates its restore listing, prevents overlapping runs, requires more than 5 GiB free and retains seven completed dumps. The first dump was restored into a separate disposable database: all 11 job rows matched, and 29 retained live artifact files matched the restored manifests. The generated restore database was then removed.

This proves logical database restore against retained artifacts. It does **not** prove a complete off-host database/artifact archive restore or OVH snapshot restore. Historical production backup remains outstanding. See [backup-restore.json](qualification/ovh-2026-11-03/backup-restore.json).

## Hetzner access and remaining production work

The authenticated Hetzner console identifies `moleculens-prod`, server `125033108`, project `13948104`, CCX23, four vCPUs, 16 GB RAM, 160 GB disk, Ashburn. The cloud project contains this server, two Primary IPs and one firewall; no attached volume was seen. Other account products and final billing have not been fully inventoried.

The initial read-only console inventory found four containers (Caddy, Python API, worker, PostgreSQL 16), with bind-mounted data under `/home/agent/moleculens-server/var`: cache approximately 11.55 GB, PostgreSQL 67.93 MB and logs 4 KiB. The earlier console reading of 11 MiB/55 MiB was incorrect; these sizes come from the authenticated SSH inventory. They fit the new disk with reserve, subject to backup and retention accounting. See [initial inventory](qualification/ovh-2026-11-03/legacy-initial-inventory.json).

The user personally reset the old root password. Root SSH remains disabled. The existing admin is `agent`. The user approved adding this Mac's existing key; SSH and Docker administration now work, preserving the original key. Sudo still requires a password, so application work uses the account's existing Docker access. The trusted provider console independently established ED25519 fingerprint `SHA256:LwGCFWlEeT8upiEzt3GOg4xOCa/306Zufm2jDNtOGIE`, matching the network key. A task-specific pinned trust file preserves historical global known-host entries unchanged. See [host identity](qualification/ovh-2026-11-03/legacy-host-identity.json). The installed key fingerprint matches the Mac's existing key. The old deployment is revision `c9fce3c5d9f69212cc65d39e06ae3c18ef74eafb`, with 439 completed and nine failed jobs, no active work, and 439 cache entries. The database has been backed up off-host and restored into an isolated OVH database; all public job/cache rows match byte for byte. Historical artifact backup is in progress. See [SSH inventory](qualification/ovh-2026-11-03/legacy-inventory.json) and [database restore](qualification/ovh-2026-11-03/legacy-database-restore.json).

Before retiring Hetzner:

1. Establish migration access; inventory deployed services, active jobs, DB schema/roles, artifact paths/sizes, every supervisor and live API route. Consumer references to `/v1/figure` require checking because the new API qualifies four other operations.
2. Back up production DB, artifacts and deployment configuration; verify an isolated restore and capacity. Preserve original bytes and published job IDs.
3. Compare paired frontend revision `7c8254ab465032098f3ca8d55b76b895e76747fe` with the actual deployment and deploy required compatibility changes. Local checks are not production frontend proof.
4. Follow finite admission drain, unfenced-process shutdown, explicit schema ownership and compatible rollback in [the runbook](go-release-runbook.md). Keep new and legacy artifact roots separate, with legacy storage read-only.
5. Rehearse ingress and rollback, configure `api.moleculens.com` and change the intended Cloudflare record. Verify HTTPS, actual frontend requests, historical results and restart behavior.
6. After production acceptance, recovery verification and account inventory, retire superseded Hetzner resources and close the account as requested. Preserve invoices and any unrelated resources.

The conservative `gc` operations profile and daily artifact-GC units were added; the live dry run exited successfully with no eligible attempts and left no maintenance container. Timer activation follows final persistence checks.

Private launch scripts/logs currently reside in ignored `build/ovh-vps1/` and `build/ovh-live/`. The initial password helper is obsolete. No push or merge was performed. Colima `moleculens-qualification` is stopped; no Brev resource was used. OVH and old production Hetzner remain running for their current tasks.
