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
from api.main import app

client = TestClient(app)


def test_cube_endpoint_returns_jsx():
    resp = client.get("/geometry/cube")
    assert resp.status_code == 200
    data = resp.json()
    assert "jsx" in data and "<Canvas" in data["jsx"]


def test_html_page_returns_html():
    resp = client.get("/geometry/html-test-page")
    assert resp.status_code == 200
    data = resp.json()
    assert "html" in data and "<html" in data["html"].lower()
