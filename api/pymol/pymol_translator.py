import logging
from typing import Any

from api.llm.openai_provider import generate_structured

from . import pymol_templates
from .scene_spec import SceneSpec

logger = logging.getLogger(__name__)

try:
    from ..llm.openai_client import get_client
except ImportError:  # pragma: no cover - fallback when package context missing
    import importlib.util
    import os

    llm_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../llm/openai_client.py"))
    spec = importlib.util.spec_from_file_location("openai_client", llm_path)
    assert spec is not None, "Could not create module spec for openai_client"
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None, "Module spec has no loader"
    spec.loader.exec_module(module)
    get_client = module.get_client

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

    # Note: client and functions are defined but not used in current implementation
    # They may be needed for future function calling approach
    # client = get_client()
    # functions = [
    #     {
    #         "name": "build_scene_request",
    #         "parameters": SceneSpec.model_json_schema(),
    #         "description": "Return a SceneSpec describing the requested PyMOL operation.",
    #     }
    # ]

    try:
        logger.info("Calling OpenAI API...")
        # Generate structured SceneSpec via LLM
        scene_spec: SceneSpec = generate_structured(
            user_prompt=prompt,
            response_model=SceneSpec,
            system_prompt="You are a helpful assistant that translates a natural-language prompt into a PyMOL command.",
        )
        logger.info("LLM response received")
        return scene_spec

    except Exception as e:
        logger.error(f"Error in LLM translation: {str(e)}")
        raise e


def translate_with_options(prompt: str) -> tuple[list[str], dict[str, Any]]:
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


def translate(prompt: str) -> list[str]:
    """Translate a natural-language prompt into PyMOL commands (backward compatibility)."""
    commands, _ = translate_with_options(prompt)
    return commands
