"""Benchmark the guarded nomenclature adapter against a versioned corpus."""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

from moleculens.compute.nomenclature import (
    NomenclatureError,
    UnsupportedNomenclatureError,
    explain_nomenclature,
)

_MIN_CORPUS_SIZE = 300
_MAX_CORPUS_SIZE = 500


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("corpus", type=Path, help="Versioned nomenclature corpus JSON")
    parser.add_argument(
        "--verify-opsin",
        action="store_true",
        help="Run OPSIN structural round-trip checks (requires Java)",
    )
    parser.add_argument("--output", type=Path, help="Optional JSON report path")
    parser.add_argument(
        "--enforce",
        action="store_true",
        help="Exit non-zero when the documented rollout gates are not met",
    )
    parser.add_argument("--min-name-rate", type=float, default=0.90)
    parser.add_argument("--min-opsin-match-rate", type=float, default=0.95)
    return parser


def _load_cases(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    cases = document.get("cases")
    if not isinstance(cases, list) or not all(isinstance(case, dict) for case in cases):
        raise ValueError("Corpus must contain a cases array")
    return document, cases


def run_benchmark(
    document: dict[str, Any],
    cases: list[dict[str, Any]],
    *,
    verify_opsin: bool,
    min_name_rate: float,
    min_opsin_match_rate: float,
) -> dict[str, Any]:
    """Run cases through the same adapter used by the API and return a compact report."""

    started = time.perf_counter()
    outcome_counts: Counter[str] = Counter()
    opsin_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()
    failures: list[dict[str, Any]] = []
    supported_total = 0
    supported_named = 0
    supported_opsin_matched = 0
    guarded_total = 0
    guarded_rejected = 0

    for case in cases:
        case_id = str(case.get("id", "unknown"))
        smiles = str(case.get("smiles", ""))
        category = str(case.get("category", "uncategorized"))
        expected = str(case.get("expected", "supported"))
        category_counts[category] += 1
        if expected == "guarded":
            guarded_total += 1
        else:
            supported_total += 1

        try:
            result = explain_nomenclature(smiles, verify_opsin=verify_opsin)
            outcome_counts["named"] += 1
            opsin_counts[result.opsin_status] += 1
            if verify_opsin and result.opsin_status != "matched":
                failures.append(
                    {
                        "id": case_id,
                        "category": category,
                        "kind": "opsin_non_match",
                        "status": result.opsin_status,
                        "displayName": result.display_name,
                        "smiles": smiles,
                    }
                )
            if expected == "guarded":
                failures.append(
                    {
                        "id": case_id,
                        "category": category,
                        "kind": "guardrail_miss",
                        "displayName": result.display_name,
                    }
                )
            else:
                supported_named += 1
                if result.opsin_status == "matched":
                    supported_opsin_matched += 1
        except UnsupportedNomenclatureError as exc:
            outcome_counts["guarded"] += 1
            if expected == "guarded":
                guarded_rejected += 1
            else:
                failures.append(
                    {
                        "id": case_id,
                        "category": category,
                        "kind": "unexpected_guardrail",
                        "code": exc.code,
                        "message": str(exc),
                    }
                )
        except NomenclatureError as exc:
            outcome_counts["naming_failed"] += 1
            failures.append(
                {
                    "id": case_id,
                    "category": category,
                    "kind": "naming_failed",
                    "message": str(exc),
                }
            )
        except Exception as exc:  # pragma: no cover - defensive benchmark boundary
            outcome_counts["engine_error"] += 1
            failures.append(
                {
                    "id": case_id,
                    "category": category,
                    "kind": "engine_error",
                    "message": str(exc),
                }
            )

    name_rate = supported_named / supported_total if supported_total else 0.0
    guarded_rate = guarded_rejected / guarded_total if guarded_total else 1.0
    opsin_match_rate = (
        supported_opsin_matched / supported_named if verify_opsin and supported_named else None
    )

    gates = {
        "corpusSize": {
            "passed": _MIN_CORPUS_SIZE <= len(cases) <= _MAX_CORPUS_SIZE,
            "actual": len(cases),
            "required": f"{_MIN_CORPUS_SIZE}-{_MAX_CORPUS_SIZE}",
        },
        "supportedNameRate": {
            "passed": name_rate >= min_name_rate,
            "actual": round(name_rate, 6),
            "required": min_name_rate,
        },
        "guardrailRecall": {
            "passed": guarded_rate == 1.0,
            "actual": round(guarded_rate, 6),
            "required": 1.0,
        },
        "opsinMatchRate": {
            "passed": bool(
                verify_opsin
                and opsin_match_rate is not None
                and opsin_match_rate >= min_opsin_match_rate
            ),
            "actual": round(opsin_match_rate, 6) if opsin_match_rate is not None else None,
            "required": min_opsin_match_rate,
            "measured": opsin_match_rate is not None,
        },
    }
    return {
        "corpusVersion": document.get("version"),
        "caseCount": len(cases),
        "verifyOpsin": verify_opsin,
        "durationSeconds": round(time.perf_counter() - started, 3),
        "outcomes": dict(sorted(outcome_counts.items())),
        "opsinStatuses": dict(sorted(opsin_counts.items())),
        "categories": dict(sorted(category_counts.items())),
        "supportedNameRate": round(name_rate, 6),
        "guardrailRecall": round(guarded_rate, 6),
        "opsinMatchRate": round(opsin_match_rate, 6) if opsin_match_rate is not None else None,
        "gates": gates,
        "failures": failures,
    }


def main() -> None:
    """CLI entry point."""

    args = _parser().parse_args()
    document, cases = _load_cases(args.corpus)
    report = run_benchmark(
        document,
        cases,
        verify_opsin=args.verify_opsin,
        min_name_rate=args.min_name_rate,
        min_opsin_match_rate=args.min_opsin_match_rate,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if args.enforce and not all(gate["passed"] for gate in report["gates"].values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
