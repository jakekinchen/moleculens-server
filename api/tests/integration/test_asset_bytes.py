import importlib.util
import sys
import types
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient


@pytest.mark.integration
def test_render_model_bytes(tmp_path, monkeypatch):
    """Render a model using stubbed PyMOL and ensure bytes are returned."""

    # Stub minimal pymol.cmd used by the router
    pymol_stub = types.ModuleType("pymol")

    class _CmdStub:
        def reinitialize(self):
            pass

        def save(self, path: str, *args, **kwargs):
            Path(path).write_text("ATOM\n" * 300)

        def png(self, path: str, *args, **kwargs):
            Path(path).write_text("PNG")

        def do(self, *args, **kwargs):
            pass

    pymol_stub.cmd = _CmdStub()
    monkeypatch.setitem(sys.modules, "pymol", pymol_stub)

    # Stub translator and utilities
    translator_stub = types.ModuleType("api.agent_management.pymol_translator")
    translator_stub.translate = lambda desc: ["cmd.save('out.pdb')"]
    agent_root = types.ModuleType("api.agent_management")
    agent_root.pymol_translator = translator_stub
    monkeypatch.setitem(sys.modules, "agent_management", agent_root)
    monkeypatch.setitem(
        sys.modules, "agent_management.pymol_translator", translator_stub
    )
    monkeypatch.setitem(sys.modules, "api.agent_management", agent_root)
    monkeypatch.setitem(
        sys.modules, "api.agent_management.pymol_translator", translator_stub
    )

    utils_stub = types.ModuleType("api.utils")
    utils_stub.cache = types.SimpleNamespace(
        CACHE_DIR=tmp_path, get=lambda key: None, set=lambda key, meta: None
    )
    utils_stub.security = types.SimpleNamespace(validate_commands=lambda c: None)
    monkeypatch.setitem(sys.modules, "api.utils", utils_stub)

    # Load router after stubs are in place
    root = Path(__file__).resolve().parents[3]
    spec = importlib.util.spec_from_file_location(
        "render_routes", root / "api" / "routers" / "render" / "routes.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)

    app = FastAPI()
    app.include_router(module.router)

    with TestClient(app) as client:
        resp = client.post("/render", json={"description": "demo", "format": "model"})
        assert resp.status_code == 200
        assert len(resp.content) > 1024
