from __future__ import annotations

import json
from typing import List, Any

import openai

from .pymol_templates import overview_scene, binding_site_scene, mutation_scene
from .scene_spec import SceneSpec

_DISPATCH = {
    "overview": overview_scene,
    "binding_site": binding_site_scene,
    "mutation": mutation_scene,
}


def _spec_from_prompt(prompt: str) -> SceneSpec:
    """Convert free-form prompt to validated SceneSpec via OpenAI function calling."""
    client = openai.OpenAI()
    functions = [
        {
            "name": "build_scene_request",
            "parameters": SceneSpec.model_json_schema(),
            "description": "Return a SceneSpec describing the requested PyMOL operation.",
        }
    ]
    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "user", "content": prompt}],
        functions=functions,
        function_call="auto",
        temperature=0,
    )
    data = json.loads(resp.choices[0].message.function_call.arguments)
    return SceneSpec(**data)


def translate(prompt: str) -> List[str]:
    """Return PyMOL command list for *prompt*."""
    spec = _spec_from_prompt(prompt)

    if spec.op == "raw":
        return spec.raw_cmds or []

    builder = _DISPATCH[spec.op]
    kwargs: dict[str, Any] = {"structure_id": spec.structure_id}

    if spec.selection:
        key = "selection" if spec.op == "binding_site" else "mutation_selection"
        kwargs[key] = spec.selection

    kwargs.update(spec.opts)
    return builder(**kwargs)
