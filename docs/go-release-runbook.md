# Go release and cutover runbook

Live follow-up: the purchased OVH candidate is deployed and qualified on native Linux amd64; production remains on Hetzner. The old host's changed SSH identity was independently verified through the authenticated provider console, with a separate pinned trust file. See [the current launch record](ovh-launch-2026-11-03.md) for live evidence and remaining gates. Historical local qualification does not establish production data or frontend migration. Do not use StrictHostKeyChecking=no, delete historical known_hosts entries, or infer trust from a network-presented key alone.

## Build and qualification

Use Go 1.27.3, linux/amd64, the pinned Psi4 base, hashed Python lock, exact Java/Debian package versions and checksum-verified xTB 6.7.1 release. Run:

```sh
uv sync --frozen --extra dev
uv run ruff check .
uv run ruff format --check .
uv run pytest tests/unit -q
TEST_DATABASE_URL='<isolated UTF8 PostgreSQL admin DSN>' go test -race ./...
go vet ./...
scripts/build_go_release.sh
```

Tests create and drop their own uniquely named databases. Never use production credentials for qualification. `build/release.json` records the source revision, dirty flag, Go binary hash and exact produced image IDs. A dirty build is not a release. Package those exact images with `docker save` or publish them only through an explicitly authorized registry route; deploy their immutable IDs/digests. A missing pinned apt package must fail the build, not silently upgrade the engine. Any updated scientific dependency/native platform triggers compatibility-version review and qualification.

The CI workflow runs PostgreSQL race tests and native HTTP/xTB/naming checks. Its definition is delivered here; no remote CI run is claimed. The old automatic SSH deployment workflow and legacy deploy script are retired in this branch. The replacement manual workflow packages images without contacting production.

## Before touching the existing host

1. Verify the SSH fingerprint through a trusted existing control plane or console and record the evidence. Inventory every API, worker, timer, supervisor, replica, DB role and artifact mount; source defaults do not prove deployed topology.
2. Record deployed image IDs, PostgreSQL version/encoding/schema, queue counts by state, artifact root and ownership, disk reserve, CPU/RAM quotas, active process trees and current client revision. Preserve a database backup plus artifact snapshot and verify restore on an isolated target.
3. Review and deploy paired consumer application `7c8254ab465032098f3ca8d55b76b895e76747fe` from branch `agent/go-consumer-compatibility-20261103`. It integrates ESP defaults, orbital metadata and shared conformer options/cache identity. Static checks, build, 52 focused tests and seven actual HTTP integration tests pass; the latter also pass against final server `e2abc7b` through Caddy. The earlier patch is superseded by this real integration. See `docs/qualification/vps1/consumer-pair.json`; no consumer deployment has occurred.
4. Establish the candidate's actual CPU/RAM budget from representative workloads on the intended architecture. The optional full two-CPU/4-GiB local trial passes the measured small-molecule envelope in `docs/vps1-capacity.md`; its Ubuntu 24.04 ARM/Rosetta timings are not an OVH benchmark. The subsequently purchased OVH Ubuntu 26.04 candidate passed native fixture capacity and HTTPS checks, documented in the launch record. Production data fit and actual traffic remain separate gates.

## Finite legacy drain and explicit schema owner

1. Stop admission at the existing ingress while preserving GET/poll/artifact reads. Record the timestamp and pending/running job IDs. Stop queue consumers from claiming more work. Choose and record a finite completion deadline for already running jobs; do not extend it indefinitely.
2. Let those jobs finish until the deadline. Stop/reap all remaining **unfenced legacy processes**, including children and every supervisor that could restart them. Verify process/container inventory is empty and artifact writes have ceased across an observation interval. Where distinct DB roles exist, revoke legacy worker write privileges as an additional fence; never assume roles are distinct until inventoried.
3. Only after that proof, finalize interrupted active legacy rows with an explicit reason. Keep original IDs, request metadata and completed artifacts. Do not mark a still-writing legacy job failed and start a replacement against the same output directory.
4. Audit duplicate active cache keys. The migration intentionally fails if these remain; reconcile them after verified drain, preserving the historical records. No migration silently deletes jobs.
5. Run the **single** explicit schema owner: `moleculens migrate`. API and worker startup only check schema. It transactionally installs checksum-tracked migrations, adds the legacy attempt fence/active-key constraint and creates the separate `moleculens_v2` queue. Require UTF8. Take a migration receipt.
6. Before any legacy work is resumed, use the hardened Python revision from this branch with its fenced completion/failure and private attempt directories. Never restart the old unfenced worker image against the migrated deployment.

## Isolated cutover, then live candidate

Create a **new** artifact root. Do not point Go at the legacy cache as its write root. Set owner UID/GID 65532 for the new root only; mount legacy storage read-only. Retain both roots through cutover and rollback. Start one qualified worker container and the Go API on a private port, using the same image revision and compatibility map. `compose.go.yml` supplies bounded candidate settings and requires explicit image/DB/root values; it does not provision PostgreSQL or change DNS.

For a complete small-host stack, use the separately qualified optional `compose.vps1.yml` and its [capacity/operation notes](vps1-capacity.md). It retains one long and one short worker, preserving naming/conformer responsiveness during Psi4. Do not reduce the normal profile limits and assume equivalent evidence.

Verify `/ready`, all native operation fixtures, cache identity distinctions, charged xTB translation, naming/atom mapping, duplicate admission, stale completion rejection, missing artifact behavior, process-group termination, DB outage and result readback. Rehearse swapping the ingress target and returning to a quiesced compatible API while reading the same IDs and artifact bytes. Protocol tests and an actual local container replacement drill pass: the candidate API/worker stopped, a separate quiesced compatible API took the same local port, both Go and legacy IDs remained readable, and their artifact SHA256 values were unchanged. The final-image SIGKILL/restart drill also recovered the same job under a new attempt token. Receipts are in `docs/qualification`. Production ingress, service manager, DB role revocation and backup restore remain host-specific release gates.

After those gates, route a bounded canary to the Go API. Monitor admission count, running attempts, leases, native exit errors, HTTP 503/429 rates, storage reserve, CPU/RAM, latency and scientific result comparisons. Expand only with a recorded acceptance result. The user subsequently authorized the OVH hosting migration and retirement of superseded Hetzner hosting; domain and billing changes still follow the acceptance and recovery gates in the launch record.

## Rollback

Stop new admission first (`api --quiesced` on the retained compatible Go image or equivalent ingress routing). Stop/reap candidate workers and verify the processes ended. Retain the Go queue schema, unique indexes, attempt fences, both artifact roots and every published manifest. Serve existing Go IDs with the last qualified compatible Go API; its read-only legacy support also serves old IDs.

A rollback that restores an unfenced old Python worker, drops the uniqueness/fencing constraints, reuses legacy output directories, deletes the new schema, or abandons new Go IDs is invalid. If legacy admission must resume, use the hardened legacy revision and retain the compatible Go GET/artifact route for already-issued new IDs. That mixed routing needs a host-specific rehearsal; the default rollback is quiesced Go until forward repair.

## Retention and cleanup

Each worker container owns one attempt process group. Keep the container/cgroup lifecycle so a worker crash terminates descendants; bare detached processes are not the release topology. `moleculens gc` defaults to a dry run and never removes referenced completed outputs or ambiguous publications. Investigate an `publication-uncertain` marker by reading the SQL row/token/manifest; preserve it until the outcome is established. Alerts, explicit retention and storage sizing must prevent the bounded free-space guard becoming a recurring outage.

No Brev resources were used. The local PostgreSQL/Colima qualification environment is disposable and must be stopped after verification; its native results and compact receipts are retained in this branch, not live databases or secrets.
