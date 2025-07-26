from . import graphic, prompt, rcsb, render

# Note: geometry module doesn't exist, commented out
# from .geometry.routes import router as geometry_router
from .graphic.routes import router as graphic_router
from .prompt.routes import router as prompt_router
from .rcsb.routes import router as rcsb_router
from .render.routes import router as render_router

# Expose all router subpackages for convenient access
__all__ = ["prompt", "render", "rcsb", "graphic", "graphic_router", "prompt_router", "rcsb_router", "render_router"]
