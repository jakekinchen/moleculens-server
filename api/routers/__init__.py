from . import geometry, graphic, prompt, rcsb, render

# Expose all router subpackages for convenient access
__all__ = ["prompt", "geometry", "render", "rcsb", "graphic"]

from .geometry.routes import router as geometry_router
from .graphic.routes import router as graphic_router

# Ensure submodules are properly exposed
from .prompt.routes import router as prompt_router
from .rcsb.routes import router as rcsb_router
from .render.routes import router as render_router
