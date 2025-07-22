import importlib.util
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.append(str(ROOT))

spec = importlib.util.spec_from_file_location(
    "api.agent_management.pymol_templates",
    ROOT / "api" / "agent_management" / "pymol_templates.py",
)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)
sys.modules["api.agent_management.pymol_templates"] = module
mutation_scene = module.mutation_scene
mutation_focus_scene = module.mutation_focus_scene


def test_mutation_scene_commands():
    cmds = mutation_scene("1ubq", "resi 50 and chain A")
    # Basic sanity checks
    assert isinstance(cmds, list)
    assert cmds[0].startswith("fetch 1ubq")
    # Ensure mutation selection is present
    assert any("mutation_site" in cmd for cmd in cmds)


def test_mutation_focus_scene_commands():
    cmds = mutation_focus_scene("1ubq", "resi 25")
    assert isinstance(cmds, list)
    assert cmds[0].startswith("fetch 1ubq")
    assert any("zoom" in cmd for cmd in cmds)
