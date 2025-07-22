from __future__ import annotations

import importlib.util
import json
import os
import sys
from typing import Any, List

try:
    from ..utils.openai_client import get_client
except ImportError:
    # Fallback for direct module loading in tests
    import os
    import sys

    sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
    from utils.openai_client import get_client


# Additional fallback: create client directly if get_client fails
def _get_openai_client():
    """Get OpenAI client with fallback."""
    try:
        return get_client()
    except Exception:
        # Direct import and creation as fallback
        import os

        from openai import OpenAI

        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            raise ValueError("OpenAI API key is required")
        return OpenAI(api_key=api_key)


# Import pymol_templates directly
templates_path = os.path.join(os.path.dirname(__file__), "pymol_templates.py")
spec = importlib.util.spec_from_file_location("pymol_templates", templates_path)
pymol_templates = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pymol_templates)

# Import scene_spec directly
scene_spec_path = os.path.join(os.path.dirname(__file__), "scene_spec.py")
spec = importlib.util.spec_from_file_location("scene_spec", scene_spec_path)
scene_spec = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scene_spec)

_DISPATCH = {
    "overview": pymol_templates.overview_scene,
    "binding_site": pymol_templates.binding_site_scene,
    "mutation": pymol_templates.mutation_scene,
    "mutation_focus": pymol_templates.mutation_focus_scene,
}


def _spec_from_prompt(prompt: str) -> scene_spec.SceneSpec:
    """Convert free-form prompt to validated SceneSpec via OpenAI function calling."""
    print("\n" + "=" * 80)
    print("TRANSLATING PROMPT:")
    print(prompt)
    print("-" * 80)

    client = _get_openai_client()
    functions = [
        {
            "name": "build_scene_request",
            "parameters": scene_spec.SceneSpec.model_json_schema(),
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
        return scene_spec.SceneSpec(**data)

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
