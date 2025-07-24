import os
import sys

from fastapi.testclient import TestClient

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import types

chem_mod = types.ModuleType("rdkit.Chem")
chem_mod.Fragments = object  # type: ignore[attr-defined]
chem_mod.Descriptors = object  # type: ignore[attr-defined]
chem_mod.AllChem = object  # type: ignore[attr-defined]
rdkit_mod = types.ModuleType("rdkit")
rdkit_mod.Chem = chem_mod  # type: ignore[attr-defined]
sys.modules["rdkit"] = rdkit_mod
sys.modules["rdkit.Chem"] = chem_mod

# Stub out the heavy `openai` dependency used in some modules
openai_mod = types.ModuleType("openai")
openai_mod.OpenAI = object  # type: ignore[attr-defined]
types_mod = types.ModuleType("openai.types")
chat_mod = types.ModuleType("openai.types.chat")
chat_completion_mod = types.ModuleType("openai.types.chat.chat_completion")
setattr(types_mod, "Completion", object)
setattr(chat_mod, "ChatCompletion", object)
setattr(chat_mod, "ChatCompletionMessage", object)
setattr(chat_mod, "ChatCompletionMessageParam", dict)
setattr(chat_completion_mod, "Choice", object)
sys.modules["openai"] = openai_mod
sys.modules["openai.types"] = types_mod
sys.modules["openai.types.chat"] = chat_mod
sys.modules["openai.types.chat.chat_completion"] = chat_completion_mod
from agent_management.agents.rcsb_agent import RCSBAgent

from api.main import app

client = TestClient(app)


def test_fetch_structure_pdb():
    # stub the agent method to avoid network
    original = RCSBAgent.fetch_structure

    def fake(self, identifier: str, file_format: str = "pdb") -> str:
        assert identifier == "1STP"
        assert file_format == "pdb"
        return "ATOM  GENERATED MOCK DATA"

    RCSBAgent.fetch_structure = fake
    try:
        resp = client.post(
            "/rcsb/fetch-structure/", json={"identifier": "1STP", "format": "pdb"}
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "ATOM" in data["data"]
    finally:
        RCSBAgent.fetch_structure = original
