import json
import logging
from typing import Any, List

logger = logging.getLogger(__name__)

try:
    from ..utils.openai_client import get_client
except ImportError:  # pragma: no cover - fallback when package context missing
    import importlib.util
    import os

    utils_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../utils/openai_client.py")
    )
    spec = importlib.util.spec_from_file_location("openai_client", utils_path)
    assert spec is not None, "Could not create module spec for openai_client"
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None, "Module spec has no loader"
    spec.loader.exec_module(module)
    get_client = module.get_client
from . import pymol_templates
from .scene_spec import SceneSpec

_DISPATCH = {
    "overview": pymol_templates.overview_scene,
    "binding_site": pymol_templates.binding_site_scene,
    "mutation": pymol_templates.mutation_scene,
    "mutation_focus": pymol_templates.mutation_focus_scene,
    "transparent_molecule": pymol_templates.transparent_molecule_scene,
    "publication_quality": pymol_templates.publication_quality_scene,
    "annotated_molecule": pymol_templates.annotated_molecule_scene,
    "transparent_binding_site": pymol_templates.transparent_binding_site_scene,
}


def _spec_from_prompt(prompt: str) -> SceneSpec:
    """Convert free-form prompt to validated SceneSpec via OpenAI function calling."""
    logger.info("=" * 80)
    logger.info("TRANSLATING PROMPT:")
    logger.info(prompt)
    logger.info("-" * 80)

    client = get_client()
    functions = [
        {
            "name": "build_scene_request",
            "parameters": SceneSpec.model_json_schema(),
            "description": "Return a SceneSpec describing the requested PyMOL operation.",
        }
    ]

    try:
        logger.info("Calling OpenAI API...")
        completion = client.chat.completions.create(  # type: ignore[call-overload]
            model="o3-mini",  # Using a version known to work well with function calling
            messages=[{"role": "user", "content": prompt}],
            tools=[{"type": "function", "function": functions[0]}],
            tool_choice={
                "type": "function",
                "function": {"name": "build_scene_request"},
            },
            # Note: o3-mini doesn't support temperature parameter
        )

        logger.info("LLM response received")
        if not completion.choices[0].message.tool_calls:
            raise ValueError("LLM did not return tool call")
        tool_call = completion.choices[0].message.tool_calls[0]
        data = json.loads(tool_call.function.arguments)
        logger.info("Generated spec:")
        logger.info(json.dumps(data, indent=2))
        return SceneSpec(**data)

    except Exception as e:
        logger.error(f"Error in LLM translation: {str(e)}")
        return SceneSpec(
            op="raw",
            structure_id="1ubq",
            raw_cmds=["cmd.fetch('1ubq')", "cmd.show('cartoon')"],
        )


def translate_with_options(prompt: str) -> tuple[List[str], dict[str, Any]]:
    """Translate a natural-language prompt into PyMOL commands and rendering options.

    Returns
    -------
    tuple[List[str], dict[str, Any]]
        PyMOL commands and rendering options dictionary
    """
    spec = _spec_from_prompt(prompt)
    logger.info(f"Operation type: {spec.op}")

    if spec.op == "raw":
        commands = spec.raw_cmds or []
    else:
        builder = _DISPATCH[spec.op]
        kwargs: dict[str, Any] = {"structure_id": spec.structure_id}

        if spec.selection:
            key = "selection" if spec.op == "binding_site" else "mutation_selection"
            kwargs[key] = spec.selection

        kwargs.update(spec.opts)
        commands = builder(**kwargs)  # type: ignore[operator]

    # Extract rendering options from the spec
    rendering_opts = {
        "transparent_background": spec.rendering.transparent_background,
        "ray_trace": spec.rendering.ray_trace,
        "resolution": spec.rendering.resolution,
        "dpi": spec.rendering.dpi,
        "ray_trace_mode": spec.rendering.ray_trace_mode,
        "antialias": spec.rendering.antialias,
        "ray_shadow": spec.rendering.ray_shadow,
        "depth_cue": spec.rendering.depth_cue,
        "background_color": spec.rendering.background_color,
    }

    logger.info("Generated PyMOL commands:")
    for cmd in commands:
        logger.info(f"  {cmd}")
    logger.info(f"Rendering options: {rendering_opts}")
    logger.info("=" * 80)
    return commands, rendering_opts


def translate(prompt: str) -> List[str]:
    """Translate a natural-language prompt into PyMOL commands (backward compatibility)."""
    commands, _ = translate_with_options(prompt)
    return commands
