# Psi4 hosting assessment

Date checked: November 3, 2026. This is a shortlist and qualification plan, not a provisioned or benchmarked migration.

## Recommendation

Evaluate lower-cost CPU virtual machines before changing the application architecture. Keep the current Psi4 computation and qualify the same pinned container on each candidate. Start with OVHcloud's 8–12 GB VPS offerings for low or intermittent load, and compare a netcup dedicated-resource VM if calculations are sustained or shared-CPU variability is unacceptable. The actual workload and peak memory remain unknown, so 8–12 GB describes candidate sizes, not a certified requirement.

A hosting move can happen separately from the proposed Go migration. The worker already runs from a `psi4/psi4` Docker base. Go service ownership does not remove the need to package and operate the chemistry engine.

## Current evidence

- Updated by the user after this assessment: the provider is **Hetzner**, not Hostinger; the estimated bill is approximately $40/month. Plan, invoice, machine location, resource allocation and utilization remain unconfirmed.
- The inspected single-box Compose file runs the API, PostgreSQL, worker and Caddy together, with persistent cache/database directories and separate Psi4 scratch.
- The source/Compose defaults allocate Psi4 one thread and 512 MB, while the README describes three threads and 4096 MB. Live environment overrides were not inspected. This is a configuration question worth resolving before buying more CPU.
- Live public DNS returned `5.161.53.230` for `api.moleculens.com`. [RIPE RDAP](https://rdap.db.ripe.net/ip/5.161.53.230) identifies the address range as `CLOUD-ASH`, country `US`, with registrant Hetzner Online GmbH. This establishes the public endpoint's network allocation, not the bill payer, backend worker location or complete hosting topology.
- Public `/health` returned `healthy`, version `1.0.0`, database/cache flags true. As established in the earlier audit, this does not qualify Psi4 computation or worker health.
- SSH inspection was blocked by a changed host key. No trust records were changed and no server commands executed. Machine specs, disk use, job volume and current engine identity therefore remain unverified.

## Current advertised options

Prices are starting/list prices as displayed on the linked pages. Region, tax, contract length, availability, backups and other options must be included in the actual quote; advertised CPU counts are not a Psi4 benchmark.

| Option | Advertised resources | Displayed starting price | Assessment |
| --- | --- | --- | --- |
| [OVHcloud VPS-2](https://us.ovhcloud.com/vps/) | 4 vCores, 8 GB RAM, 75 GB NVMe, daily backup | $8.50/month | First low-cost candidate; measure sustained CPU and largest supported jobs |
| [OVHcloud VPS-3](https://us.ovhcloud.com/vps/) | 6 vCores, 12 GB RAM, 100 GB NVMe, daily backup | $12.32/month | More RAM headroom for the same single-box layout; still needs workload qualification |
| [netcup RS 1000 G12.5](https://www.netcup.com/en/server/root-server) | 4 dedicated cores as advertised, 8 GB RAM, 128 GB NVMe | From €21.73/month on the displayed 12-month offer | Candidate for predictable CPU allocation. Quote the actual location, term and tax treatment before comparing with dollars |
| [DigitalOcean Basic](https://www.digitalocean.com/pricing/droplets) | 4 shared vCPUs, 8 GiB RAM, 160 GiB SSD | $48/month | Does not reduce the reported $40 baseline at this size; dedicated CPU offerings cost more |

[Hetzner's cost-optimized CX33](https://www.hetzner.com/cloud/cost-optimized/) advertises 4 vCPUs and 8 GB, but the retrieved official page marked it unavailable and did not render a price. Its stated fit is variable CPU performance and low-to-medium use. Do not promise a move to an unavailable cheap plan or use older pricing as a quote. The current endpoint's Hetzner allocation makes identifying the existing account/plan especially useful: a right-size may be simpler than a provider move.

At a confirmed $40/month baseline, an $8.50–$12.32 base VM would save $27.68–$31.50/month before tax/options and migration effort, **if** it meets the same service requirements. This is conditional arithmetic, not established equivalent performance.

## Psi4-specific considerations

[Psi4's runtime documentation](https://psicode.org/psi4manual/master/external.html) describes CPU threading through BLAS/LAPACK, OpenMP and C++ threads, and recommends fast local scratch storage. Prioritize real CPU performance, sufficient RAM and scratch capacity/throughput. There is no evidence that the inspected SCF/cubeprop path benefits from renting a GPU.

Use a Linux x86-64 candidate first to minimize a simultaneous architecture/packaging change; an ARM option needs separate image and numerical qualification. Leave memory for PostgreSQL, API, OS and mesh post-processing rather than assigning all host RAM to Psi4. Measure peak process RSS and scratch bytes, since a Psi4 memory setting alone is not a total-container memory guarantee.

Pay-per-job hosting is a later option if demand is very sparse. For illustration, [Cloud Run Jobs pricing](https://cloud.google.com/run/pricing) at the displayed Iowa rates is $0.000018 per vCPU-second plus $0.000002 per GiB-second. A 4-vCPU/8-GiB allocation for 300 billed seconds is $0.0264 before free allowances and other services. That is allocation arithmetic, not a measured job runtime. The API, database, persistent artifacts, startup and transfer costs still count. The [default writable filesystem](https://docs.cloud.google.com/run/docs/container-contract#file_system) uses memory and is ephemeral, so scratch and persistence require an explicit design. The current shared-disk worker is not a direct unchanged Cloud Run deployment.

## Qualification before moving

1. Verify the real provider, invoice, machine specs and trusted SSH host identity; read effective Psi4 threads/memory without exposing secrets.
2. Capture representative small, typical and largest supported jobs, plus observed job frequency, cache-hit rate, peak RAM/scratch and queue time. No production stress test is needed.
3. Pin the same engine image and run those fixtures on candidates. Compare numerical outputs, meshes, failure behavior, wall time, peak memory, scratch use and sustained shared-CPU variability. Include full host costs and backups in the comparison.
4. If a candidate meets the service requirements at lower total cost, perform a backup/restore check, stop new submissions during a bounded drain, preserve old job IDs/artifacts, and switch with a tested rollback. Retire the old paid server only after the new service and data are verified.

No infrastructure was purchased, created, stopped or migrated during this assessment. No Brev resources were used.
