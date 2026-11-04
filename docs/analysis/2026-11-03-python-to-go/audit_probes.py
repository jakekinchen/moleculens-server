"""Disposable audit probes: assertions confirm current defects, not desired behavior.

Run from any directory with the repository virtualenv:
  python -m pytest /absolute/path/to/audit_probes.py -q -p no:cacheprovider
All database and result-file writes use pytest temporary directories.
"""

import asyncio
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime, timedelta
from threading import Barrier

import numpy as np
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from moleculens.api import electrostatics_routes, routes
from moleculens.compute.charges import calculate_charges_gasteiger
from moleculens.compute.sdf_parser import Atom, SDFMolecule, SDFParseError, parse_sdf
from moleculens.compute.surface_mesh import generate_molecular_surface
from moleculens.db import queue as q
from moleculens.db.models import Base, Job, JobStatus


@pytest.fixture
def queue(tmp_path, monkeypatch):
    engine = create_engine(f"sqlite:///{tmp_path / 'audit.sqlite'}")
    Base.metadata.create_all(engine)
    monkeypatch.setattr(q, "_engine", engine)
    monkeypatch.setattr(q, "_session_factory", sessionmaker(bind=engine))
    monkeypatch.setattr(q.settings, "molecule_cache_dir", tmp_path)
    queue = q.JobQueue(tmp_path)
    monkeypatch.setattr(routes, "job_queue", queue)
    yield queue
    engine.dispose()


def submit(queue, key="audit", kind="orbital"):
    return queue.submit_job_with_cache_key(key, {"job_type": kind})[0]


def test_concurrent_submission_creates_two_active_jobs(queue, monkeypatch):
    barrier = Barrier(2)
    original = queue._get_reusable_active_job_id

    def synchronized_lookup(*args, **kwargs):
        result = original(*args, **kwargs)
        barrier.wait(timeout=5)
        return result

    monkeypatch.setattr(queue, "_get_reusable_active_job_id", synchronized_lookup)
    with ThreadPoolExecutor(max_workers=2) as workers:
        ids = list(workers.map(lambda _: submit(queue), range(2)))
    assert len(set(ids)) == 2
    assert all(queue.get_job(i).status == JobStatus.PENDING for i in ids)


def test_stale_attempt_can_complete_after_replacement(queue):
    old_id = submit(queue)
    with q.db_session() as session:
        old = session.get(Job, old_id)
        old.status = JobStatus.RUNNING
        old.started_at = datetime.now(UTC) - timedelta(seconds=q.settings.job_stale_seconds + 1)
    new_id = submit(queue)
    assert old_id != new_id
    assert queue.get_job(old_id).status == JobStatus.ERROR
    queue.complete_job(old_id, {"old": True}, ["result.json"], 1)
    assert queue.get_job(old_id).status == JobStatus.DONE
    assert queue.get_job(old_id).error_message is not None
    assert queue.get_job(new_id).status == JobStatus.PENDING
    assert queue.get_cache_entry("audit").result_meta == {"old": True}


def test_complete_accepts_missing_artifact(queue):
    job_id = submit(queue)
    queue.complete_job(job_id, {"any": True}, ["result.json"], 1)
    assert queue.get_job(job_id).status == JobStatus.DONE
    assert not (queue.cache_dir / "audit" / "result.json").exists()


def test_empty_json_is_loaded_as_successful_result(queue):
    root = queue.cache_dir / "audit"
    root.mkdir()
    (root / "result.json").write_text("{}")
    (root / "electrostatics_result.json").write_text("{}")
    orbital = routes._load_result("audit", {})
    esp = electrostatics_routes._load_electrostatics_result("audit", {})
    assert orbital is not None and orbital.orbitals == {}
    assert esp is not None and esp.surface.mesh.vertexCount == 0
    assert esp.dipole.magnitude_debye == 0


def test_wrong_job_route_invalidates_completed_conformer(queue):
    job_id = submit(queue, kind="conformer")
    root = queue.cache_dir / "audit"
    root.mkdir()
    (root / "structure.sdf").write_text("present conformer artifact")
    queue.complete_job(job_id, {"quality": "embedded"}, ["structure.sdf"], 1)
    response = asyncio.run(routes.get_job_status(job_id))
    assert response.status == "error"
    assert queue.get_job(job_id).status == JobStatus.ERROR
    assert queue.get_cache_entry("audit") is None


def test_client_identity_reuses_distinct_inputs(queue):
    kwargs = {
        "method": "scf",
        "basis": "sto-3g",
        "grid_spacing": 0.25,
        "isovalue": 0.05,
        "orbitals": ["homo"],
        "cache_identity": "same-client-key",
    }
    first = queue.submit_job(sdf_content="geometry A", **kwargs)[0]
    second = queue.submit_job(sdf_content="geometry B", **kwargs)[0]
    assert first == second
    assert queue.get_job(second).request_params["sdf_content"] == "geometry A"


def test_electrostatic_geometry_hash_ignores_bond_graph():
    from rdkit import Chem

    blocks = []
    partial_charges = []
    for smiles in ["CO", "C=O"]:
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", "audit")
        conf = Chem.Conformer(2)
        conf.SetAtomPosition(0, (0, 0, 0))
        conf.SetAtomPosition(1, (1.2, 0, 0))
        mol.AddConformer(conf)
        sdf = Chem.MolToMolBlock(mol)
        parsed = parse_sdf(sdf)
        blocks.append(parsed.geometry_hash())
        partial_charges.append(calculate_charges_gasteiger(parsed, sdf).charges)
    assert blocks[0] == blocks[1]
    assert not np.allclose(partial_charges[0], partial_charges[1])


def test_surface_coordinates_do_not_match_sampling_grid():
    _, verts = generate_molecular_surface(np.array([[0.0, 0.0, 0.0]]), ["C"])
    midpoint = (verts.min(axis=0) + verts.max(axis=0)) / 2
    print("sphere midpoint (expected zero):", midpoint.tolist())
    assert np.linalg.norm(midpoint) > 0.01


def test_ionic_dipole_does_not_match_reported_origin():
    from rdkit import Chem

    mol = Chem.MolFromSmiles("[Na+]")
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, (10, 0, 0))
    mol.AddConformer(conf)
    molecule = SDFMolecule("ion", [Atom("Na", 10, 0, 0)])
    result = calculate_charges_gasteiger(molecule, Chem.MolToMolBlock(mol))
    assert result.charges[0] == 1
    assert np.array_equal(result.dipole_origin, [10, 0, 0])
    assert result.dipole_magnitude > 48  # About declared origin it would be zero.


def test_valid_unnamed_rdkit_molblock_fails_sdf_parser():
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CO")
    block = Chem.MolToMolBlock(mol)
    assert block.startswith("\n")
    assert Chem.MolFromMolBlock(block) is not None
    with pytest.raises(SDFParseError):
        parse_sdf(block)


def test_health_failure_still_has_success_http_status(monkeypatch):
    from fastapi.testclient import TestClient

    from moleculens.api.main import app

    def unavailable():
        raise RuntimeError("isolated injected DB outage")

    monkeypatch.setattr(routes, "get_db_engine", unavailable)
    # No context manager: do not run application's database startup handler.
    response = TestClient(app).get("/health")
    assert response.status_code == 200
    assert response.json()["status"] == "unhealthy"


def test_client_options_are_silently_discarded():
    from moleculens.api.schemas import ComputeRequest, ConformerRequest

    orbital = ComputeRequest.model_validate(
        {
            "sdf_content": "placeholder sdf",
            "charge": 1,
            "multiplicity": 2,
            "density_isovalue": 0.01,
        }
    )
    assert not {"charge", "multiplicity", "density_isovalue"} & orbital.model_dump().keys()
    conformer = ConformerRequest.model_validate(
        {
            "molblock2d": "placeholder molblock",
            "geomVersion": "etkdgv3",
            "seed": 42,
            "numConformers": 3,
            "optimize": False,
            "forceField": "mmff",
        }
    )
    assert conformer.params.opt == "uff"
    assert not {"seed", "numConformers", "optimize", "forceField"} & conformer.model_dump().keys()


def test_fixed_width_counts_are_misread():
    from rdkit import Chem

    mol = Chem.MolFromSmiles("C" * 101)
    mol.SetProp("_Name", "audit")
    block = Chem.MolToMolBlock(mol)
    assert block.splitlines()[3].startswith("101100")
    assert Chem.MolFromMolBlock(block).GetNumAtoms() == 101
    with pytest.raises(SDFParseError):
        parse_sdf(block)


def test_element_and_metal_detection_misclassifies():
    from rdkit import Chem

    from moleculens.compute.charges import _element_to_number
    from moleculens.compute.conformers import _has_metal

    assert _element_to_number("Au") == 6
    assert _has_metal(Chem, Chem.MolFromSmiles("[Fe]")) is False
