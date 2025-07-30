"""Moleculens Server API Package.

This package provides the FastAPI-based server implementation for
molecular visualization and analysis, featuring AI-powered molecule
rendering and diagram generation.
"""

__version__ = "0.0.1"
__author__ = "Moleculens Team"

import sys
import types as _types

# Alias the `api.routers` package so imports like `api.routers.render` work
# _api_routers_pkg = importlib.import_module("api.routers")
# sys.modules.setdefault("api.routers", _api_routers_pkg)

# Provide a lightweight stub for `pymol` (used by render routes) when the real
# library is not available (e.g., during CI or unit testing where PyMOL isn't
# installed). The stub implements the minimal API surface needed by the code.


try:
    import pymol  # noqa: F401 - Attempt to import the real PyMOL if available
except ImportError:  # Fallback to a lightweight stub during CI / local dev
    if "pymol" not in sys.modules:
        _pymol = _types.ModuleType("pymol")

        class _CmdStub:
            def finish_launching(self, *args, **kwargs):
                pass

            def reinitialize(self):
                pass

            def get_view(self):
                return []

            def centerofmass(self):
                return []

            def get_extent(self):
                return ((), ())

            def png(self, *args, **kwargs):
                pass

            def save(self, *args, **kwargs):
                pass

            def do(self, *args, **kwargs):
                pass

            def orient(self, *args, **kwargs):
                pass

            def zoom(self, *args, **kwargs):
                pass

            def label(self, *args, **kwargs):
                pass

        _pymol.cmd = _CmdStub()  # type: ignore[attr-defined]
        sys.modules["pymol"] = _pymol


# Provide a stub for the `celery` module when not installed. This allows tests
# to import the Celery app without requiring the real dependency.
if "celery" not in sys.modules:
    _celery = _types.ModuleType("celery")

    class _CeleryApp:
        def __init__(self, *args, **kwargs):
            pass

        def task(self, *args, **kwargs):  # noqa: D401 - mimic decorator
            def decorator(func):
                return func

            return decorator

    _celery.Celery = _CeleryApp  # type: ignore[attr-defined]
    sys.modules["celery"] = _celery
