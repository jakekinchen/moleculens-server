# Python-to-Go analysis and Opus agreement

Analysis date: November 3, 2026. This directory contains analysis and reproducible audit evidence; it does not implement the migration.

Start with [the two-way agreement](agreement.md) and [the consolidated analysis](analysis.md). They cover the migration boundary, all six code lenses, current and proposed state transitions, ranked findings, stages, qualification gates and rollback.

## Review record

- [Opus's independent report](opus-independent.md) is preserved byte for byte. It was written before Opus read the Codex analysis. Its initial recommendations are historical proposals; the final reconciliation supersedes any differences.
- [The concrete reconciliation proposal](reconciliation-proposal.md) records the seven decisions sent back for challenge, along with evidence limits and corrections.
- [Opus's reconciliation response](opus-response.md) records the accepted amendments and evidence corrections.
- [Opus's final confirmation](opus-final-confirmation.md) closes the two specification choices. Codex accepted its final clarifications in [the agreement record](agreement.md).

The review used the [code-lenses skill](/Users/kelly/.codex/skills/code-lenses/SKILL.md), [state-machine skill](/Users/kelly/.codex/skills/state-machine/SKILL.md), and [advisors skill](/Users/kelly/.codex/skills/advisors-skill/SKILL.md), with the shared [advisor prompting guide](/Users/kelly/.codex/skills/advisors-skill/references/advisor-prompting.md). The advisor was Opus 5.5 Extra in Claude Desktop's local Code session `local_182bddf0-7b25-413f-801e-2340256df144`.

## Evidence

- [Validation receipt](validation.md) distinguishes executed checks from unexecuted release qualification.
- [Server snapshot](source-snapshot.json) identifies the exact dirty checkout inspected; [consumer snapshot](consumer-snapshot.json) identifies the narrow local frontend sample.
- [Audit probes](audit_probes.py) reproduce fourteen current defects. Passing means the undesirable behavior was observed. They are outside the ordinary test suite and are not desired-behavior regression tests.
- [Initial probe output](audit-probes.log) and [additional probe output](audit-probes-additional.log) preserve the executed results.
- [Naming baseline](nomenclature-baseline.json) and [its run output](nomenclature-baseline.log) cover 398 cases without OPSIN. They do not establish structural round-trip correctness.

No production runtime, native Psi4/xTB computation, PostgreSQL contention test, Java/OPSIN gate or Go build was executed. Existing source edits and runtime data were preserved.
