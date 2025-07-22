import json
from typing import Any, List

try:
    from ..utils.openai_client import get_client
except ImportError:  # pragma: no cover - fallback when package context missing
    import importlib.util
    import os

    utils_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../utils/openai_client.py")
    )
    spec = importlib.util.spec_from_file_location("openai_client", utils_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    get_client = module.get_client
from . import pymol_templates
from .scene_spec import SceneSpec

_DISPATCH = {
    "overview": pymol_templates.overview_scene,
    "binding_site": pymol_templates.binding_site_scene,
    "mutation": pymol_templates.mutation_scene,
    "mutation_focus": pymol_templates.mutation_focus_scene,
}


def _spec_from_prompt(prompt: str) -> SceneSpec:
    """Convert free-form prompt to validated SceneSpec via OpenAI function calling."""
    print("\n" + "=" * 80)
    print("TRANSLATING PROMPT:")
    print(prompt)
    print("-" * 80)

    client = get_client()
    functions = [
        {
            "name": "build_scene_request",
            "parameters": SceneSpec.model_json_schema(),
            "description": "Return a SceneSpec describing the requested PyMOL operation.",
        }
    ]

    try:
        print("Calling OpenAI API...")
        completion = client.chat.completions.create(
            model="o3-mini",  # Using a version known to work well with function calling
            messages=[{"role": "user", "content": prompt}],
            functions=functions,
            function_call={"name": "build_scene_request"},  # Force function call
            # Note: o3-mini doesn't support temperature parameter
        )

        print("LLM response received")
        data = json.loads(completion.choices[0].message.function_call.arguments)
        print("\nGenerated spec:")
        print(json.dumps(data, indent=2))
        return SceneSpec(**data)

    except Exception as e:
        print(f"\nError in LLM translation: {str(e)}")
        raise


def translate(prompt: str) -> List[str]:
    """Return PyMOL command list for *prompt*."""
    spec = _spec_from_prompt(prompt)
    print(f"\nOperation type: {spec.op}")

    if spec.op == "raw":
        commands = spec.raw_cmds or []
    else:
        builder = _DISPATCH[spec.op]
        kwargs: dict[str, Any] = {"structure_id": spec.structure_id}

        if spec.selection:
            key = "selection" if spec.op == "binding_site" else "mutation_selection"
            kwargs[key] = spec.selection

        kwargs.update(spec.opts)
        commands = builder(**kwargs)

    print("\nGenerated PyMOL commands:")
    for cmd in commands:
        print(f"  {cmd}")
    print("=" * 80)
    return commands
