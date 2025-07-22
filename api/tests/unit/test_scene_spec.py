import importlib.util
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.append(str(ROOT))

spec = importlib.util.spec_from_file_location(
    "api.agent_management.scene_spec",
    ROOT / "api" / "agent_management" / "scene_spec.py",
)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)
SceneSpec = module.SceneSpec


def test_scene_spec_validation():
    spec = SceneSpec(
        op="mutation",
        structure_id="1ubq",
        selection="resi 50 and chain A",
        opts={"original_residue": "F"},
    )
    assert spec.op == "mutation"
    assert spec.selection.startswith("resi")


def test_scene_spec_mutation_focus():
    spec = SceneSpec(op="mutation_focus", structure_id="1abc", selection="resi 12")
    assert spec.op == "mutation_focus"
