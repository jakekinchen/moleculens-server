# Validation receipt

Date: 2026-11-03. Server HEAD: `c9fce3c5d9f69212cc65d39e06ae3c18ef74eafb`, plus existing user changes recorded in `source-snapshot.json`.

Codex ran:

1. `.venv/bin/python -m pytest tests/unit -q -p no:cacheprovider`
   - 46 passed in 4.47 seconds.
   - One warning: Java unavailable for py2opsin.
   - Existing database tests use SQLite. No production database was used.
2. `python -m pytest .../audit_probes.py -q -s -p no:cacheprovider`
   - Initial 11 probes passed in 4.50 seconds, from `/tmp` using the repository virtualenv.
   - Output: `audit-probes.log`.
   - These assertions reproduce defects; passing means the defect was observed.
3. `python -m pytest .../audit_probes.py -q -p no:cacheprovider -k 'client_options or fixed_width_counts or element_and_metal'`
   - Three added probes passed in 1.64 seconds, 11 previous probes deselected.
   - Output: `audit-probes-additional.log`.
4. `.venv/bin/python -m ruff check src tests scripts --output-format concise`
   - All checks passed.
5. Naming benchmark without OPSIN, launched from `/tmp`:
   - `python -m moleculens.compute.nomenclature_benchmark .../benchmarks/nomenclature_corpus.json --output .../nomenclature-baseline.json`
   - 398 cases in 6.444 seconds: 392 names, six guarded rejections.
   - Supported-name rate and guardrail recall: 1.0 each.
   - OPSIN not requested; its release gate is **not passed or measured** by this run.
   - Outputs: `nomenclature-baseline.json`, `nomenclature-baseline.log`.
6. Compared current server source/config/fixture bytes to the recorded source snapshot after the probes.
   - No mismatches at that check.

Local packages observed: Python 3.11.14, FastAPI 0.128.0, SQLAlchemy 2.0.45, NumPy 2.4.0, scikit-image 0.26.0, RDKit 2026.3.4, Openclatura 0.1.2, pytest 9.0.2.

Unexecuted: real PostgreSQL locking/races, full Psi4/xTB integration, Java/OPSIN round trips, Go build/runtime, production load/canary/browser behavior. Psi4 and xTB Python modules are absent; Java is unavailable; the Docker daemon is not running; Go was not found on PATH. No installations or service restarts were needed for this analysis.

The local consumer sample is separately identified in `consumer-snapshot.json`. Its three API schema files and overlay cache-key helper are clean at the recorded HEAD; the nomenclature caller is untracked. No consumer implementation was changed or executed.

All temporary database and artifact writes from audit probes use disposable pytest directories. The checked-in-style probe artifact is intentionally under this analysis directory, outside `tests/`; it must not be added to production CI as a desired-behavior test suite. When implementing fixes, convert each relevant probe to a regression asserting the corrected behavior.

Opus's independently executed probes and their limits are attributed in `opus-independent.md` and `opus-response.md`. Its final agreement is in `opus-final-confirmation.md`. Those original reports are preserved byte for byte; `agreement.md` records the reconciled decisions and corrections. Opus's local startup and ASGI interference measurements are not production benchmarks.

The final source snapshot comparison checked all 58 recorded files with no mismatches. Final lint also covers the disposable audit probe file. No runtime code, deployment, commit or push was performed.
