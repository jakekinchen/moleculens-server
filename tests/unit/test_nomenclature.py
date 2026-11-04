"""Tests for the guarded Openclatura nomenclature adapter."""

from pathlib import Path

import pytest
from rdkit import Chem

from moleculens.api.nomenclature_routes import explain
from moleculens.api.schemas import NomenclatureRequest
from moleculens.compute.nomenclature import (
    UnsupportedNomenclatureError,
    _safe_py2opsin,
    explain_nomenclature,
)


def _molblock(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return Chem.MolToMolBlock(mol)


def test_explains_name_and_maps_atoms_to_viewer_order() -> None:
    smiles = "CC(=O)Nc1ccccc1"
    response = explain_nomenclature(smiles, molblock=_molblock(smiles), verify_opsin=False)

    assert response.display_name == "N-phenylacetamide"
    assert response.source == "openclatura"
    assert response.mapping_status == "matched"
    assert response.opsin_status == "not_requested"
    assert response.rules_hit
    assert response.components
    assert any(component.phase == "prefix" for component in response.components)
    assert any(component.phase == "suffix" for component in response.components)
    assert all("atom " not in component.text.lower() for component in response.components)
    assert response.atom_locants
    assert any(step.atoms for step in response.trace)
    assert all(atom >= 0 for step in response.trace for atom in step.atoms)


def test_maps_an_arbitrary_viewer_atom_permutation() -> None:
    response = explain_nomenclature(
        "CCO",
        molblock=_molblock("OCC"),
        verify_opsin=False,
    )

    parent_step = next(step for step in response.trace if step.phase == "parent_selection")
    assert response.mapping_status == "matched"
    assert parent_step.engine_atoms == [0, 1]
    assert parent_step.atoms == [1, 2]
    assert response.atom_locants == {"1": 1, "2": 2}


def test_never_exposes_engine_atoms_as_viewer_atoms_when_mapping_fails() -> None:
    response = explain_nomenclature(
        "CCO",
        molblock=_molblock("c1ccccc1"),
        verify_opsin=False,
    )

    assert response.mapping_status == "unmatched"
    assert response.atom_locants == {}
    assert any(step.engine_atoms for step in response.trace)
    assert all(step.atoms == [] for step in response.trace)


def test_rejects_partial_substructure_matches_as_viewer_mappings() -> None:
    response = explain_nomenclature(
        "CCO",
        molblock=_molblock("CCO.C"),
        verify_opsin=False,
    )

    assert response.mapping_status == "unmatched"
    assert response.atom_locants == {}
    assert all(step.atoms == [] for step in response.trace)


@pytest.mark.parametrize(
    ("smiles", "code"),
    [
        ("CCO.CN", "disconnected_structure"),
        ("N[Pt](Cl)(Cl)N", "metal_structure"),
        ("[BiH3]", "metal_structure"),
    ],
)
def test_rejects_guarded_domains(smiles: str, code: str) -> None:
    with pytest.raises(UnsupportedNomenclatureError) as exc_info:
        explain_nomenclature(smiles, verify_opsin=False)

    assert exc_info.value.code == code


def test_opsin_wrapper_uses_unique_writable_tempfile() -> None:
    captured: list[str] = []

    def delegate(names: list[str], *, tmp_fpath: str) -> list[str]:
        path = Path(tmp_fpath)
        assert path.exists()
        path.write_text("ethanol", encoding="utf-8")
        captured.append(tmp_fpath)
        return names

    wrapped = _safe_py2opsin(delegate)
    assert wrapped(["ethanol"]) == ["ethanol"]
    assert wrapped(["propan-1-ol"]) == ["propan-1-ol"]

    assert len(set(captured)) == 2
    assert all(not Path(path).exists() for path in captured)


def test_route_serializes_camel_case_contract() -> None:
    response = explain(
        NomenclatureRequest(
            smiles="CC(=O)Nc1ccccc1",
            molblock=_molblock("CC(=O)Nc1ccccc1"),
            verifyOpsin=False,
        )
    )
    payload = response.model_dump(by_alias=True)

    assert payload["displayName"] == "N-phenylacetamide"
    assert payload["source"] == "openclatura"
    assert payload["mappingStatus"] == "matched"
    assert payload["opsinStatus"] == "not_requested"
    assert payload["engineVersion"] == "0.1.2"
    assert payload["trace"][0]["engineAtoms"]
