"""Build the versioned nomenclature benchmark corpus from locked upstream cases."""

from __future__ import annotations

import ast
import json
from importlib.metadata import version
from pathlib import Path
from typing import TypedDict

import openclatura
from rdkit import Chem


class CorpusCase(TypedDict):
    id: str
    label: str
    smiles: str
    category: str
    expected: str
    source: str


CURATED_CASES: tuple[tuple[str, str, str, str], ...] = (
    ("water", "O", "moleculens_catalog", "supported"),
    ("ethanol", "CCO", "moleculens_catalog", "supported"),
    ("glucose", "OCC1OC(O)C(O)C(O)C1O", "moleculens_catalog", "supported"),
    (
        "sucrose",
        "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O",
        "moleculens_catalog",
        "supported",
    ),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "moleculens_catalog", "supported"),
    ("acetaminophen", "CC(=O)Nc1ccc(O)cc1", "moleculens_catalog", "supported"),
    ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "moleculens_catalog", "supported"),
    ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "moleculens_catalog", "supported"),
    ("benzene", "c1ccccc1", "moleculens_catalog", "supported"),
    ("acetic-acid", "CC(=O)O", "moleculens_catalog", "supported"),
    ("citric-acid", "O=C(O)CC(O)(CC(=O)O)C(=O)O", "moleculens_catalog", "supported"),
    ("lactic-acid", "CC(O)C(=O)O", "moleculens_catalog", "supported"),
    ("urea", "NC(N)=O", "moleculens_catalog", "supported"),
    ("glycine", "NCC(=O)O", "moleculens_catalog", "supported"),
    ("dopamine", "NCCc1ccc(O)c(O)c1", "moleculens_catalog", "supported"),
    ("serotonin", "NCCc1c[nH]c2ccc(O)cc12", "moleculens_catalog", "supported"),
    (
        "cholesterol",
        "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
        "moleculens_catalog",
        "supported",
    ),
    ("salicylic-acid", "O=C(O)c1ccccc1O", "moleculens_catalog", "supported"),
    (
        "penicillin-g",
        "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
        "moleculens_catalog",
        "supported",
    ),
    ("n-phenylacetamide", "CC(=O)Nc1ccccc1", "name_builder", "supported"),
    ("l-alanine", "C[C@H](N)C(=O)O", "stereochemistry", "supported"),
    ("benzyl-alcohol", "c1ccccc1CO", "moleculens_catalog", "supported"),
    ("2-3-dimethylbutane", "CC(C)C(C)C", "fallback_smiles", "supported"),
    ("acetate", "CC(=O)[O-]", "charged", "supported"),
    ("methylammonium", "C[NH3+]", "charged", "supported"),
    ("nitrobenzene", "[O-][N+](=O)c1ccccc1", "charged", "supported"),
    ("chiral-halide", "F[C@H](Cl)Br", "stereochemistry", "supported"),
    ("trans-2-butene", "C/C=C/C", "stereochemistry", "supported"),
    ("pyridine", "n1ccccc1", "heterocycle", "supported"),
    ("imidazole", "c1ncc[nH]1", "heterocycle", "supported"),
    ("sodium-chloride", "[Na+].[Cl-]", "guardrail_salt", "guarded"),
    ("ammonium-chloride", "[NH4+].[Cl-]", "guardrail_salt", "guarded"),
    ("ethanol-methylamine-mixture", "CCO.CN", "guardrail_mixture", "guarded"),
    ("cisplatin", "N[Pt](Cl)(Cl)N", "guardrail_metal", "guarded"),
    ("iron-ion", "[Fe+2]", "guardrail_metal", "guarded"),
    ("bismuthane", "[BiH3]", "guardrail_metal", "guarded"),
)


def _canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid curated SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _upstream_cases() -> list[CorpusCase]:
    root = Path(openclatura.__file__).resolve().parent / "tests_roundtrip"
    if not root.is_dir():
        raise RuntimeError("The locked Openclatura wheel does not include its round-trip corpus")

    cases: list[CorpusCase] = []
    for path in sorted(root.glob("test_*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        values: list[str] = []
        for node in tree.body:
            if not isinstance(node, ast.Assign):
                continue
            if not any(
                isinstance(target, ast.Name) and target.id == "SMILES" for target in node.targets
            ):
                continue
            if isinstance(node.value, (ast.List, ast.Tuple)):
                values.extend(
                    item.value
                    for item in node.value.elts
                    if isinstance(item, ast.Constant) and isinstance(item.value, str)
                )
        category = path.stem.removeprefix("test_")
        for index, smiles in enumerate(values, start=1):
            _canonical(smiles)
            cases.append(
                {
                    "id": f"openclatura-{category}-{index:03d}",
                    "label": f"Openclatura {category} case {index}",
                    "smiles": smiles,
                    "category": category,
                    "expected": "supported",
                    "source": f"openclatura=={version('openclatura')}:{path.name}",
                }
            )
    return cases


def main() -> None:
    destination = Path(__file__).resolve().parents[1] / "benchmarks/nomenclature_corpus.json"
    seen: set[str] = set()
    cases: list[CorpusCase] = []
    for label, smiles, category, expected in CURATED_CASES:
        canonical = _canonical(smiles)
        if canonical in seen:
            continue
        seen.add(canonical)
        cases.append(
            {
                "id": f"moleculens-{label}",
                "label": label.replace("-", " ").title(),
                "smiles": smiles,
                "category": category,
                "expected": expected,
                "source": "moleculens curated product and guardrail cases",
            }
        )
    for case in _upstream_cases():
        canonical = _canonical(case["smiles"])
        if canonical in seen:
            continue
        seen.add(canonical)
        cases.append(case)

    if not 300 <= len(cases) <= 500:
        raise RuntimeError(f"Expected 300-500 unique cases, found {len(cases)}")
    document = {
        "version": 1,
        "description": "Initial MolecuLens nomenclature rollout corpus",
        "openclaturaVersion": version("openclatura"),
        "methodology": (
            "Curated MolecuLens catalog, fallback, stereochemistry, charged, and guardrail cases "
            "plus the locked Openclatura round-trip corpus, deduplicated by canonical isomeric SMILES."
        ),
        "cases": cases,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(cases)} cases to {destination}")


if __name__ == "__main__":
    main()
