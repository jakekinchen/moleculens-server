from . import (  # noqa: F401 – ensures /sequence-coordinates route is added
    sequence_coordinates,
)
from .routes import router

__all__ = ["router", "sequence_coordinates"]
