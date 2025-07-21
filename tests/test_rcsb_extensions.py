import os
import sys
import types

from fastapi.testclient import TestClient

# stub rdkit to avoid heavy imports
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
chem_mod = types.ModuleType("rdkit.Chem")
chem_mod.Fragments = object
chem_mod.Descriptors = object
chem_mod.AllChem = object
rdkit_mod = types.ModuleType("rdkit")
rdkit_mod.Chem = chem_mod
sys.modules["rdkit"] = rdkit_mod
sys.modules["rdkit.Chem"] = chem_mod

from agent_management.agents.rcsb_agent import RCSBAgent

from api.main import app

client = TestClient(app)


def test_fetch_model_endpoint():
    # stub the agent method to avoid network
    original = RCSBAgent.fetch_alphafold_model

    def fake(self, uniprot_id: str, file_format: str = "pdb") -> str:
        assert uniprot_id == "P12345"
        assert file_format == "pdb"
        return "MODEL_DATA"

    RCSBAgent.fetch_alphafold_model = fake
    try:
        resp = client.post(
            "/rcsb/fetch-model/", json={"uniprot_id": "P12345", "format": "pdb"}
        )
        assert resp.status_code == 200
        assert resp.json()["data"] == "MODEL_DATA"
    finally:
        RCSBAgent.fetch_alphafold_model = original


def test_entry_metadata_endpoint():
    # stub the agent method to avoid network
    original = RCSBAgent.fetch_entry_metadata

    def fake_meta(self, identifier: str) -> dict:
        assert identifier == "1ABC"
        return {"rcsb_id": identifier}

    RCSBAgent.fetch_entry_metadata = fake_meta
    try:
        resp = client.get("/rcsb/entry/1ABC")
        assert resp.status_code == 200
        assert resp.json()["metadata"]["rcsb_id"] == "1ABC"
    finally:
        RCSBAgent.fetch_entry_metadata = original
