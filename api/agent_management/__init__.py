from . import (
    agent_factory,
    agent_model_config,
    llm_service,
    llm_service_extension,
    model_config,
    models,
    molecule_visualizer,
    pymol_templates,
    pymol_translator,
    scene_packager,
    scene_spec,
)

try:
    from . import code_checker
except ImportError:
    # code_checker requires playwright which may not be available in all environments
    code_checker = None  # type: ignore[assignment]
from . import debug_utils, model_registry

try:
    from . import html_to_caption_json
except ImportError:
    # html_to_caption_json requires beautifulsoup4 which may not be available
    html_to_caption_json = None  # type: ignore[assignment]
from . import diagram_renderer

__all__ = [
    "pymol_translator",
    "pymol_templates",
    "llm_service",
    "scene_packager",
    "models",
    "model_config",
    "agent_model_config",
    "agent_factory",
    "scene_spec",
    "molecule_visualizer",
    "llm_service_extension",
    "model_registry",
    "debug_utils",
    "diagram_renderer",
]

# Add optional modules to __all__ only if they were successfully imported
if code_checker is not None:
    __all__.append("code_checker")
if html_to_caption_json is not None:
    __all__.append("html_to_caption_json")
