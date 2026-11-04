# Implementation agent commission

Commissioned November 3, 2026 at the user's explicit request to start a new agent thread and work through the agreed implementation plan.

## Objective and authority

Implement the Codex–Opus agreed Python-to-Go migration for MolecuLens, preserving Psi4 as a critical capability. Start with contract fixtures and bounded correctness hardening, then implement the Go service/worker boundary and qualify it. This is an execution commission, not another request for a plan-only report.

The decision authority is [agreement.md](agreement.md), read together with [analysis.md](analysis.md) and [opus-final-confirmation.md](opus-final-confirmation.md). The original independent Opus report contains explicitly superseded recommendations; do not adopt those over the final agreement.

Read and apply the code-lenses and state-machine skills during implementation and validation. If consulting an advisor, read the shared advisor prompting guide and distinguish independently verified results from reported/unproven claims. No deadline or two-week schedule was imposed.

## Starting checkout

Repository: `/Users/kelly/Developer/moleculens-server`.

Analysis baseline HEAD: `c9fce3c5d9f69212cc65d39e06ae3c18ef74eafb`. Verify the current identity and relevant instructions before editing. [source-snapshot.json](source-snapshot.json) describes the pre-existing dirty implementation; it is part of the input, not disposable work. In particular, preserve the uncommitted nomenclature capability, benchmark corpus, nine naming tests, modified dependency/configuration files, and existing runtime data.

Use an isolated branch/worktree or checkout with the required dirty source overlay preserved. Do not reset, clean or overwrite the user's changes. Do not copy live databases/caches or credentials into fixtures or commits. Record the actual worktree and branch in a new execution status artifact so the user can locate the work.

## Execution sequence

1. Freeze actual consumer/server request and response fixtures and the adapter protocol. Inspect the narrow consumer paths identified in `consumer-snapshot.json`; coordinate changes needed for density threshold, conformer force field and ESP defaults before strict rejection.
2. Implement bounded Python correctness fixes: active-job uniqueness, guarded attempt completion/failure, operation-kind checks, correct pending-job handling, private attempt output and truthful validated results. Correct parser and scientific semantics with independent expected results. Introduce one versioned schema migration owner and reproducible engine builds.
3. Implement Go ownership of HTTP/admission, the new computation descriptor, PostgreSQL state transitions, supervision, deadlines/retries, process cleanup and immutable artifact publication. Keep pinned Python/native chemistry adapters without database credentials. Build the full supervisor once in Go.
4. Run isolated PostgreSQL contention/crash/publication tests, consumer differential checks and native Psi4/chemistry qualification. Exercise the staged cutover and rollback in a disposable environment. Prepare a concrete deployable release and migration/rollback runbook, with precise remaining live-release requirements.
5. Only pursue optional numerical/native ports if profiling and independent correctness evidence justify them. They are not prerequisites to delivering Go service ownership.

Preserve immediate invalid-SDF errors using the agreed limited structural validator. Use a new full-input descriptor and per-adapter result-compatibility version; retain exact build provenance in manifests and preserve legacy identities. Implement the agreed conservative compatibility-version bump rules.

## Critical evidence rules

- The existing unit baseline is 46 passing tests. Fourteen saved audit probes reproduce defects; convert relevant probes into regressions asserting corrected behavior. Do not count the original probes passing as fixes.
- The naming benchmark has 398 cases: 392 named and six guarded. Its existing baseline did not run OPSIN. Preserve the corpus and mapping guarantees and qualify real Java/OPSIN when available.
- At analysis time Go, native Psi4/xTB, Java, Docker daemon and PostgreSQL qualification were unavailable/unexecuted. Recheck capabilities. Install or configure standard local tools where useful and authorized, and use isolated data. Do not silently substitute mocks for native engine or real PostgreSQL proof.
- Runtime performance, complete production topology and current engine versions are unknown. Never present source inspection or unit tests as production or scientific integration proof.
- Keep progress and validation receipts by stage. Persist through the authorized implementation; report a concrete external blocker only after completing independent work that remains possible.

## Hosting and live boundaries

The user confirmed **Hetzner**, with an estimated $40/month bill. [hosting-assessment.md](hosting-assessment.md) is a separate, unbenchmarked cost shortlist; the migration must not depend on purchasing another host. Current machine specs and utilization are unverified. SSH to `api.moleculens.com` failed host-key verification during read-only inspection. Verify that identity through a trusted route; do not bypass host-key checking or replace trust records based only on the newly presented key.

This commission authorizes implementation and qualification. It does not commission paid infrastructure purchases, account/billing changes, cancellation of the existing host, or a public DNS/provider migration. Prepare those decisions concretely if needed. Preserve active jobs and artifacts; exercise cutover/rollback in an isolated environment before any live action.

No Brev resources are needed for the plan. If an authorized task later uses Brev, apply the supplied cost-control protocol: inspect the authenticated inventory before finishing and stop/delete resources no longer needed for the verified task.

## Completion criteria

Deliver source changes in the recorded worktree/branch, focused regression and integration evidence, a runnable Go-operated service retaining Psi4 and other required chemistry, compatibility fixtures, versioned migrations, documented budgets and release/cutover/rollback instructions. Clearly separate implemented, locally qualified, native-engine qualified and live-deployed states. The user should receive a concrete reviewable implementation, not a claim that the entire plan passed when external qualification is still missing.
