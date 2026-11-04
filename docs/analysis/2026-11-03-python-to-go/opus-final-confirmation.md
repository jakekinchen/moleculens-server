# Moleculens Python-to-Go: Opus final confirmation

Date: 2026-11-03. This is analysis only, with no new broad checks. It closes O1 and O2 from `/tmp/moleculens-go-opus-agreement-2026-11-03.md`. Both earlier Opus documents are preserved unchanged:
- `/tmp/moleculens-go-opus-independent-2026-11-03.md`
- `/tmp/moleculens-go-opus-agreement-2026-11-03.md`

## Confirmation

**Confirmed.** Opus and Codex agree on the architecture, the sequence and the release gates, with A1–A7 and the evidence corrections accepted by both sides. O1 and O2 are closed as stated below. No material disagreement remains. This is agreement on a plan, not release approval or performance proof.

## O2 (closed): admission validation, A1(a)

A small Go structural validator returns immediate invalid-SDF errors. It has no authority over graph identity.

**Requirements:**
- explicit supported formats;
- finite coordinates;
- known elements;
- size and count bounds;
- fixtures checked against an independent RDKit oracle.

**Precision so the fixtures catch the known parser defects:**
- Parse V2000 counts as fixed-width columns (atoms in cols 1–3, bonds in 4–6). Preserve the header line positions, so blank titles are valid.
- Treat V3000 as an explicit decision: accept it or reject it with a coded error. It must never be rejected by accident.
- Apply the validator uniformly. Today, orbital and electrostatics submits that carry a `clientCacheKey` skip parsing entirely. That is a deliberate contract change and needs its own fixture.

## O1 (closed): identity versioning

- The key holds a per-adapter **result-compatibility version**.
- The manifest records exact build provenance.
- Changes to scientific code, engines, dependencies or the numerical platform get a new version by default.
- Keeping a version across such a change requires explicit independent qualification. A small passing fixture is not enough.

This default is stricter than my "drift fixture forces a bump" proposal. I adopt it.

**One clarification.** A rebuild counts as "packaging-only" only if its *resolved* environment is identical: the lockfile-resolved package set, the engine base-image digest and native libraries. An unchanged source diff is not enough. Today's images ignore `uv.lock` and use `psi4/psi4:latest` with unpinned pip, so no current rebuild can qualify. The retention path becomes available only after stage 1 pins builds. Until then, every rebuild implies a new compatibility version, or existing entries stay legacy.

## Warm pools and drain (accepted as stated)

**Warm pools**, if measurements justify them:
- bounded single-job slots;
- attempt-private output;
- reset and state-leak tests;
- kill and replacement after timeout;
- no DB credentials.

Psi4 stays per-attempt, because its in-process option and global state is a known leak vector.

**Drain cutoff:** stop and reap legacy processes *before* fencing their rows. Then verify that the legacy artifact roots receive no further writes. The legacy worker writes unfenced into final directories, so empty database states alone prove nothing.

Also accepted: A3 (resolve the three always-sent consumer mismatches before strict activation) and a single schema-migration owner.
