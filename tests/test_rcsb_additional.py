import os
import sys
import types

from fastapi.testclient import TestClient

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

# Stub heavy rdkit import to keep tests lightweight
chem_mod = types.ModuleType("rdkit.Chem")
chem_mod.Fragments = object  # type: ignore[attr-defined]
chem_mod.Descriptors = object  # type: ignore[attr-defined]
chem_mod.AllChem = object  # type: ignore[attr-defined]
rdkit_mod = types.ModuleType("rdkit")
rdkit_mod.Chem = chem_mod  # type: ignore[attr-defined]
sys.modules["rdkit"] = rdkit_mod
sys.modules["rdkit.Chem"] = chem_mod

from agent_management.agents.rcsb_agent import RCSBAgent

from api.main import app

client = TestClient(app)


def test_sequence_annotations_endpoint():
    original = RCSBAgent.fetch_sequence_annotations

    def fake(self, identifier: str) -> dict:  # pylint: disable=unused-argument
        return {"annotations": []}

    RCSBAgent.fetch_sequence_annotations = fake
    try:
        resp = client.get("/rcsb/annotations/1XYZ")
        assert resp.status_code == 200
        assert resp.json()["annotations"] == {"annotations": []}
    finally:
        RCSBAgent.fetch_sequence_annotations = original


def test_computed_model_endpoint():
    original = RCSBAgent.fetch_graphql_model

    def fake(
        self, identifier: str, model_id: str
    ) -> dict:  # pylint: disable=unused-argument
        return {"entry": {"rcsb_id": identifier, "model_id": model_id}}

    RCSBAgent.fetch_graphql_model = fake
    try:
        resp = client.post(
            "/rcsb/computed-model/", json={"identifier": "1DEF", "model_id": "1"}
        )
        assert resp.status_code == 200
        assert resp.json()["metadata"]["entry"]["rcsb_id"] == "1DEF"
    finally:
        RCSBAgent.fetch_graphql_model = original


def test_fetch_esm_model_endpoint():
    original = RCSBAgent.fetch_esmf_model

    def fake(
        self, uniprot_id: str, file_format: str = "pdb"
    ) -> str:  # pylint: disable=unused-argument
        return "ESM_MODEL_DATA"

    RCSBAgent.fetch_esmf_model = fake
    try:
        resp = client.post(
            "/rcsb/fetch-esm-model/", json={"uniprot_id": "Q9XYZ1", "format": "pdb"}
        )
        assert resp.status_code == 200
        assert resp.json()["data"] == "ESM_MODEL_DATA"
    finally:
        RCSBAgent.fetch_esmf_model = original


def test_upload_structure_endpoint():
    original = RCSBAgent.upload_structure

    def fake(
        self, file_bytes: bytes, filename: str = "upload.pdb"
    ) -> str:  # pylint: disable=unused-argument
        return "UPLOAD_ID_123"

    RCSBAgent.upload_structure = fake
    try:
        resp = client.post(
            "/rcsb/upload-structure/",
            json={"data": "ATOM data", "filename": "mock.pdb"},
        )
        assert resp.status_code == 200
        assert resp.json()["upload_id"] == "UPLOAD_ID_123"
    finally:
        RCSBAgent.upload_structure = original
