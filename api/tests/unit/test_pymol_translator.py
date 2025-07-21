import os
import sys
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.append(str(ROOT))

spec_t = importlib.util.spec_from_file_location(
    "api.agent_management.pymol_translator", ROOT / "api" / "agent_management" / "pymol_translator.py"
)
module_t = importlib.util.module_from_spec(spec_t)
assert spec_t and spec_t.loader
spec_t.loader.exec_module(module_t)
sys.modules["api.agent_management.pymol_translator"] = module_t
translate = module_t.translate


def test_translate_mutation_prompt():
    cmd_list = translate("mutation 1ubq resi 50 and chain A")
    assert isinstance(cmd_list, list)
    assert cmd_list[0].startswith("fetch 1ubq")
    # The translator should have used the mutation template (looks for colour magenta)
    assert any("magenta" in c for c in cmd_list)
