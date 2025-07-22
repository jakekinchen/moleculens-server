from pathlib import Path

import pytest

import importlib.util
import types
import sys

CELERY_PATH = Path(__file__).resolve().parents[2] / "celery_app.py"
spec = importlib.util.spec_from_file_location("celery_app", CELERY_PATH)
celery_module = importlib.util.module_from_spec(spec)
assert spec and spec.loader

pymol_stub = types.ModuleType("pymol")
class _CmdStub:
    def reinitialize(self):
        pass
    def save(self, *args, **kwargs):
        Path(args[0]).write_text("placeholder")
    def do(self, *args, **kwargs):
        pass
pymol_stub.cmd = _CmdStub()
sys.modules["pymol"] = pymol_stub

agent_mgmt_stub = types.ModuleType("api.agent_management")
translator_stub = types.ModuleType("api.agent_management.pymol_translator")
translator_stub.translate = lambda desc: ["cmd.save('demo.obj', format='obj')"]
agent_mgmt_stub.pymol_translator = translator_stub
sys.modules["api.agent_management"] = agent_mgmt_stub
sys.modules["api.agent_management.pymol_translator"] = translator_stub

utils_stub = types.ModuleType("api.utils")
utils_stub.cache = types.SimpleNamespace(CACHE_DIR=Path("/tmp"), get=lambda k: None, set=lambda k, m: None)
utils_stub.security = types.SimpleNamespace(validate_commands=lambda c: None)
sys.modules["api.utils"] = utils_stub

spec.loader.exec_module(celery_module)
render_scene = celery_module.render_scene
from api.utils import cache


@pytest.mark.unit
def test_render_scene_creates_output(tmp_path, monkeypatch):
    monkeypatch.setattr(cache, "CACHE_DIR", tmp_path)
    path = render_scene("foo", "gltf")
    assert Path(path).exists()
    assert path.endswith(".gltf")

