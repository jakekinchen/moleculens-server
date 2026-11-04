# Codex–Opus two-way agreement

Date: November 3, 2026. Status: **agreed plan, implementation not started**.

Codex and **Opus 5.5 Extra** agree on the architecture, sequence and release gates below. Opus reviewed the local source independently, then challenged the Codex proposal, reproduced additional findings and explicitly confirmed the final decisions. Both initially open specification choices are closed. No material disagreement remains; a severity-ranking difference for wrong-endpoint polling does not change its inclusion in the prerequisite fixes.

The scope is `moleculens-server` at `c9fce3c5d9f69212cc65d39e06ae3c18ef74eafb` plus the original local changes recorded in [the source snapshot](source-snapshot.json). The selected local consumer files at `35dd5f4b4bb174b2005074b8e1e0f16a2702dd5b` are separately identified in [the consumer snapshot](consumer-snapshot.json). These identities do not establish production state.

## Decisions

| Topic | Agreed decision |
| --- | --- |
| Migration target | Go owns HTTP, admission, PostgreSQL job transitions, execution supervision, cache identity and validated artifact publication. Keep established Python/native chemistry adapters initially. |
| Rationale | The user's Go maintainability preference is sufficient to pursue this target. Measure operating behavior; make no unmeasured chemistry-speedup or Python-free claim. |
| State authority | Conditional SQL, attempt tokens and database-time leases; terminal outcomes are absorbing. Result availability is separate from job outcome. Keep existing public state projections. PostgreSQL remains the queue. |
| Publication | Private immutable attempt output is validated before committing a pointer. A late attempt cannot publish. An unknown commit acknowledgement requires identity readback; retain ambiguous output while the database is unavailable. |
| Identity | A new server-owned descriptor hashes conservatively normalized full input bytes, effective options, operation kind, fallback policy and a per-adapter result-compatibility version. Bind client hints to that digest. Preserve legacy job/cache reads under their existing identity; do not require new keys to reproduce a defective old recipe. |
| Build compatibility, formerly O1 | Record exact build provenance in manifests. Scientific code, engine, dependency or numerical-platform changes get a new compatibility version by default. Reuse across such changes requires independent qualification. Packaging-only reuse requires identical resolved scientific artifacts, not just unchanged source. Current unpinned rebuilds do not qualify. |
| Admission, formerly O2 | A small Go structural validator preserves immediate invalid-SDF errors without claiming graph-identity authority. Use fixed-width V2000 counts, valid blank headers, finite coordinates, known elements and bounded sizes. Make V3000 support/rejection explicit. Validate requests carrying a client cache key too. Check fixtures against RDKit independently. |
| Compatibility activation | Resolve ignored density, force-field and ESP-default differences in the consumer/engine contract before strict rejection. Observe unknown fields against fixtures first. Preserve known aliases and response shapes. |
| Bounded prerequisite work | Harden Python's active-job uniqueness, attempt-conditional completion/failure, kind-safe reads, pending-job handling, attempt-private output and result truth. Build the full supervisor/lease/reaper once in Go. Read-only Go and fixture work can proceed in parallel. Use one versioned SQL migration owner. |
| Execution classes | Psi4 uses per-attempt processes. Measure conformer/naming startup in the target image; a warm pool is allowed if justified, with bounded single-job slots, private output, reset/state-leak tests and kill/replacement on timeout. Chemistry adapters have no database credentials or shared publication authority. |
| Cutover | Prefer a finite drain. If overlap is required, isolate both queue schema/tables and artifact roots, and route old job IDs by origin. A cache prefix alone is insufficient. |
| Drain and rollback | At a legacy cutoff, quiesce submissions/claims, stop and reap unfenced old processes, verify writes cease, then finalize their rows. Rollback uses a pinned, fenced, compatible image and preserves accepted job reads. Never drop uniqueness/fencing to restore an older worker. |
| Scientific truth | Never publish partial/corrupt output or all-zero fallback geometry as valid results. A fallback cannot populate a strict requested-method key. Independent analytic/reference fixtures define expected science and tolerances; matching two implementations alone is insufficient. |

## Implementation order and gates

1. Freeze route and consumer contracts, corrected scientific fixtures, adapter protocol and workload/capability baselines.
2. Apply bounded correctness fixes and reproducible builds. Turn defect reproductions into desired-behavior regressions. Verify the relevant PostgreSQL contention and interrupted-writer cases. Keep unrelated Go read-only work moving in parallel.
3. Introduce the Go boundary and new computation descriptor, then transfer worker supervision and publication ownership. Qualify consumer compatibility, admission/backpressure, two-worker races, crashes, late results, uncertain commits, full disk, corruption, restart and process cleanup.
4. Perform a bounded shadow comparison and drain-first cutover. Verify chemistry, operating metrics, legacy job-ID access and an actual rollback drill.
5. Consider individual numerical/native ports only after profiling and independent numerical/topology qualification.

The full contracts, transition tables, findings and stage-specific acceptance criteria are in [the consolidated analysis](analysis.md).

## What remains to specify during implementation

Record supported charged/open-shell chemistry, valid planar-fallback policy, strict versus disclosed engine fallback, resource/retention limits and privacy boundaries before changing the affected behavior. These are implementation policy decisions with explicit gates, not unresolved architectural disagreements. Fixing degenerate coordinates and preventing fallback cache mislabeling does not wait for those policy decisions.

No two-week deadline was imposed. The user clarified that the request was for a **two-way agreement**.

## Evidence and limits

Codex ran the 46-test unit baseline, fourteen isolated defect-reproduction tests, lint, and a 398-case naming benchmark without OPSIN (392 names, six guarded rejections). Passing reproduction tests means the defects exist. Opus independently ran source-focused probes, including degenerate conformer fallback and event-loop interference, then reproduced additional Codex findings during reconciliation.

Actual PostgreSQL concurrency, native Psi4/xTB computation, Java/OPSIN round trips, Go runtime behavior, production capacity and deployed capabilities remain unqualified. No implementation, deployment or performance result is implied by this agreement. All 58 recorded source/configuration/fixture hashes remained unchanged at final verification. No Brev resources were used.

## Preserved review record

- [Independent Opus assessment](opus-independent.md), written before reading the Codex analysis.
- [Seven-decision reconciliation proposal](reconciliation-proposal.md).
- [Opus response with amendments and evidence corrections](opus-response.md).
- [Opus's final explicit confirmation](opus-final-confirmation.md), closing O1/O2 and clarifying pinned builds and legacy process shutdown.
- [Validation receipt](validation.md).

Codex accepts the final confirmation, including stopping/reaping unfenced legacy processes before finalizing their rows. This decision record and the consolidated analysis supersede conflicting initial recommendations in the preserved independent report. That report's categorical production xTB claim, mandatory legacy-key parity, unconditional unsafe rollback suggestions and negative Go-ecosystem inference were corrected during the exchange, rather than silently edited out of its history.
