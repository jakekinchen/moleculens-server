import importlib.util
import os
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.append(str(ROOT))


def _load_scene_spec() -> type:
    spec_s = importlib.util.spec_from_file_location(
        "api.agent_management.scene_spec",
        ROOT / "api" / "agent_management" / "scene_spec.py",
    )
    module = importlib.util.module_from_spec(spec_s)
    assert spec_s and spec_s.loader
    spec_s.loader.exec_module(module)
    return module.SceneSpec


def _load_translator(monkeypatch: any) -> any:
    pkg = types.ModuleType("api.agent_management")
    pkg.__path__ = [str(ROOT / "api" / "agent_management")]
    api_pkg = types.ModuleType("api")
    api_pkg.agent_management = pkg
    monkeypatch.setitem(sys.modules, "api", api_pkg)
    monkeypatch.setitem(sys.modules, "api.agent_management", pkg)

    spec_t = importlib.util.spec_from_file_location(
        "api.agent_management.pymol_translator",
        ROOT / "api" / "agent_management" / "pymol_translator.py",
    )
    module_t = importlib.util.module_from_spec(spec_t)
    assert spec_t and spec_t.loader
    spec_t.loader.exec_module(module_t)
    monkeypatch.setitem(sys.modules, "api.agent_management.pymol_translator", module_t)
    return module_t


def test_translate_with_stub(monkeypatch):
    SceneSpec = _load_scene_spec()
    module_t = _load_translator(monkeypatch)
    translate = module_t.translate

    stub_spec = SceneSpec(
        op="mutation",
        structure_id="1ubq",
        selection="resi 50 and chain A",
    )
    monkeypatch.setattr(module_t, "_spec_from_prompt", lambda _: stub_spec)

    cmds = translate("any text here")
    assert cmds[0].startswith("fetch 1ubq")
    assert any("magenta" in c for c in cmds)


def test_translate_mutation_focus(monkeypatch):
    SceneSpec = _load_scene_spec()
    module_t = _load_translator(monkeypatch)
    translate = module_t.translate

    stub_spec = SceneSpec(
        op="mutation_focus",
        structure_id="1abc",
        selection="resi 10",
    )
    monkeypatch.setattr(module_t, "_spec_from_prompt", lambda _: stub_spec)

    cmds = translate("focus")
    assert any("zoom" in c for c in cmds)
