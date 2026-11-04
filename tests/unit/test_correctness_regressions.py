"""Corrected assertions for the saved fourteen audit findings and native boundaries."""

import asyncio
import hashlib

import numpy as np
import pytest
from fastapi import HTTPException
from rdkit import Chem

from moleculens.api import electrostatics_routes, routes
from moleculens.api.schemas import ComputeRequest, ConformerRequest
from moleculens.compute.charges import _element_to_number, calculate_charges_gasteiger
from moleculens.compute.conformers import _has_metal, run_conformer_computation
from moleculens.compute.sdf_parser import SDFParseError, parse_sdf
from moleculens.compute.surface_mesh import generate_molecular_surface
from moleculens.db import queue as q
from moleculens.db.models import Job, JobStatus
from tests.unit import test_job_queue as queue_fixtures

isolated_job_queue = queue_fixtures.isolated_job_queue
reset_cache_metrics_registry = queue_fixtures.reset_cache_metrics_registry


def submit(queue, kind="orbital"):
    return queue.submit_job_with_cache_key("key", {"job_type": kind})[0]


def set_running(job_id, token):
    with q.db_session() as session:
        job = session.get(Job, job_id)
        job.status = JobStatus.RUNNING
        job.attempt_token = token


def test_stale_completion_and_failure_are_fenced(isolated_job_queue):
    job_id = submit(isolated_job_queue)
    set_running(job_id, "current")
    assert isolated_job_queue.complete_job(job_id, {}, [], 1, attempt_token="old") is False
    assert isolated_job_queue.fail_job(job_id, "late", attempt_token="old") is False
    assert isolated_job_queue.get_job(job_id).status == JobStatus.RUNNING
    assert isolated_job_queue.fail_job(job_id, "real", attempt_token="current") is True
    assert isolated_job_queue.complete_job(job_id, {}, [], 1, attempt_token="current") is False
    assert isolated_job_queue.get_job(job_id).status == JobStatus.ERROR


def test_missing_artifacts_cannot_complete(isolated_job_queue):
    job_id = submit(isolated_job_queue)
    set_running(job_id, "attempt")
    with pytest.raises(ValueError, match="Incomplete"):
        isolated_job_queue.complete_job(job_id, {}, ["result.json"], 1, attempt_token="attempt")
    assert isolated_job_queue.get_job(job_id).status == JobStatus.RUNNING


def test_empty_json_is_not_success(isolated_job_queue):
    root = isolated_job_queue.cache_dir / "key"
    root.mkdir()
    (root / "result.json").write_text("{}")
    (root / "electrostatics_result.json").write_text("{}")
    assert routes._load_result("key", {}) is None
    assert electrostatics_routes._load_electrostatics_result("key", {}) is None


def test_wrong_kind_is_read_only_404(isolated_job_queue):
    job_id = submit(isolated_job_queue, "conformer")
    with q.db_session() as session:
        session.get(Job, job_id).status = JobStatus.DONE
    with pytest.raises(HTTPException) as caught:
        asyncio.run(routes.get_job_status(job_id))
    assert caught.value.status_code == 404
    assert isolated_job_queue.get_job(job_id).status == JobStatus.DONE


def test_client_identity_cannot_alias_distinct_inputs(isolated_job_queue, water_sdf):
    args = {
        "method": "scf",
        "basis": "sto-3g",
        "grid_spacing": 0.25,
        "isovalue": 0.05,
        "orbitals": ["homo"],
        "cache_identity": "same",
    }
    first = isolated_job_queue.submit_job(sdf_content=water_sdf, **args)[0]
    second = isolated_job_queue.submit_job(
        sdf_content=water_sdf.replace("0.7586", "0.7587", 1), **args
    )[0]
    assert first != second
    with pytest.raises(SDFParseError):
        isolated_job_queue.submit_job(sdf_content="invalid", **args)


def test_full_input_hash_retains_bond_graph():
    blocks = []
    charges = []
    for smiles in ("CO", "C=O"):
        mol = Chem.MolFromSmiles(smiles)
        conf = Chem.Conformer(2)
        conf.SetAtomPosition(0, (0, 0, 0))
        conf.SetAtomPosition(1, (1.2, 0, 0))
        mol.AddConformer(conf)
        block = Chem.MolToMolBlock(mol)
        blocks.append(hashlib.sha256(block.encode()).hexdigest())
        charges.append(calculate_charges_gasteiger(parse_sdf(block), block).charges)
    assert blocks[0] != blocks[1]
    assert not np.allclose(charges[0], charges[1])


def test_surface_coordinates_match_independent_sphere():
    center = np.array([2.3, -1.7, 0.4])
    _, verts = generate_molecular_surface(np.array([center]), ["C"])
    midpoint = (verts.min(axis=0) + verts.max(axis=0)) / 2
    np.testing.assert_allclose(midpoint, center, atol=1e-5)
    radius = np.linalg.norm(verts - center, axis=1)
    # Carbon VDW 1.7 + probe 1.4. Linear interpolation error is bounded by grid.
    assert np.max(np.abs(radius - 3.1)) < 0.03


def test_ionic_dipole_matches_declared_origin():
    mol = Chem.MolFromSmiles("[NH4+]")
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, (10, 2, -3))
    mol.AddConformer(conf)
    block = Chem.MolToMolBlock(mol)
    result = calculate_charges_gasteiger(parse_sdf(block), block)
    assert result.charges.sum() == pytest.approx(1)
    np.testing.assert_equal(result.dipole_origin, [10, 2, -3])
    assert result.dipole_magnitude == 0


def test_blank_header_and_fixed_width_counts_are_valid():
    mol = Chem.MolFromSmiles("C" * 101)
    block = Chem.MolToMolBlock(mol)
    assert block.startswith("\n") and block.splitlines()[3].startswith("101100")
    assert parse_sdf(block).num_atoms == 101


def test_options_reach_engine_or_fail_explicitly(water_sdf):
    r = ComputeRequest(sdf_content=water_sdf, charge=2, density_isovalue=0.013)
    assert r.charge == 2 and r.density_isovalue == 0.013
    with pytest.raises(ValueError):
        ComputeRequest(sdf_content=water_sdf, multiplicity=2)
    c = ConformerRequest(
        molblock2d=water_sdf, geomVersion="etkdgv3", seed=42, optimize=False, forceField="mmff"
    )
    assert c.params.opt == "none" and c.params.seed == 42
    with pytest.raises(ValueError):
        ConformerRequest(molblock2d=water_sdf, geomVersion="etkdgv3", numConformers=3)


def test_element_and_metal_classification():
    assert _element_to_number("Au") == 79
    assert _has_metal(Chem, Chem.MolFromSmiles("[Fe]")) is True
    assert _has_metal(Chem, Chem.MolFromSmiles("CCO")) is False
    with pytest.raises(ValueError):
        _element_to_number("Xx")


def test_embedding_failure_preserves_original_coordinates(monkeypatch, tmp_path):
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(mol)
    block = Chem.MolToMolBlock(mol)
    expected = mol.GetConformer().GetPositions()

    def fail_embed(mol, params):
        mol.RemoveAllConformers()
        return -1

    monkeypatch.setattr(AllChem, "EmbedMolecule", fail_embed)
    result = run_conformer_computation(block, opt="none", output_dir=tmp_path)
    actual = (
        Chem.MolFromMolBlock((tmp_path / "structure.sdf").read_text()).GetConformer().GetPositions()
    )
    assert result.meta["quality"] == "planar_fallback"
    np.testing.assert_allclose(actual[:, :2], expected[:, :2], atol=5e-5)
    assert np.ptp(actual[:, 0]) > 1


def test_legacy_health_outage_returns_503(monkeypatch):
    from fastapi.testclient import TestClient

    from moleculens.api.main import app

    def unavailable():
        raise RuntimeError("isolated DB outage")

    monkeypatch.setattr(routes, "get_db_engine", unavailable)
    assert TestClient(app).get("/health").status_code == 503


def test_unsupported_surface_radius_is_not_synthesized():
    from moleculens.compute.surface_mesh import get_vdw_radius

    with pytest.raises(ValueError, match="qualified"):
        get_vdw_radius("Og")
