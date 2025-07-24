import datetime
import os
import sys
import threading
import types
from pathlib import Path
from typing import Any, Dict

import pymol
import redis

from api import routers

# Import and initialize the model registry at startup
from api.agent_management.model_config import register_models
from api.utils.rate_limit import RateLimitMiddleware
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

# Register all models
register_models()

app = FastAPI(
    title="AI Backend",
    version="0.0.1",
)

# launch PyMOL once at startup
pymol_lock = threading.Lock()


@app.on_event("startup")
def startup_pymol() -> None:
    """Initialize PyMOL in headless mode."""
    pymol.finish_launching(["pymol", "-cq"])


# Construct an absolute path to the static directory
# Assuming main.py is in the api/ directory, and static is api/static/
STATIC_DIR = Path(__file__).parent / "static"
# Ensure the directory exists (though it should from test setup or deployment)
STATIC_DIR.mkdir(parents=True, exist_ok=True)

app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

# Determine if we're in development or production mode
is_development = os.environ.get("ENVIRONMENT", "development").lower() == "development"

# Configure CORS based on environment
if is_development:
    # Development mode - use specific origins
    origins = [
        "http://localhost:3000",  # React development server
        "http://localhost:8000",  # Backend server (for serving frontend in production)
    ]
    allow_credentials = True
else:
    # Production mode - include specific origins
    origins = [
        "https://moleculens.com",  # Primary domain
        "https://www.moleculens.com",  # www subdomain
        "https://api.moleculens.com",  # api environment
    ]
    allow_credentials = True  # Enable credentials for specific origins

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=allow_credentials,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["Content-Type", "Authorization", "X-Requested-With"],
    expose_headers=["Content-Type"],
    max_age=86400,  # Cache preflight response for 24 hours
)

# rate limiting - temporarily disabled for testing
# app.add_middleware(RateLimitMiddleware)


@app.get("/health")
async def health_check() -> Dict[str, Any]:
    """Health check endpoint."""
    try:
        # Check Redis connection - temporarily disabled
        # redis_client = redis.Redis(
        #     host=os.environ.get("REDIS_HOST", "localhost"),
        #     port=int(os.environ.get("REDIS_PORT", 6379)),
        # )
        # redis_client.ping()

        # Check PyMOL
        with pymol_lock:
            pymol.cmd.reinitialize()

        return {
            "status": "healthy",
            "timestamp": datetime.datetime.now().isoformat(),
            "services": {"redis": "disabled", "pymol": "running"},
        }
    except Exception as e:
        return {
            "status": "unhealthy",
            "timestamp": datetime.datetime.now().isoformat(),
            "error": str(e),
        }


# user management related endpoints
app.include_router(routers.prompt.router)
app.include_router(routers.geometry.router)
app.include_router(routers.render.router)
app.include_router(routers.rcsb.router)
app.include_router(routers.graphic.router)

# ---------------------------------------------------------------------------
# Optional heavy dependencies
# ---------------------------------------------------------------------------
# These libraries are large binary wheels that may be absent in sandboxed or CI
# environments (e.g., Codex, Cloud runners).  Import failures would normally
# crash the process and trigger an infinite restart loop.  To prevent that we
# create minimal stub modules exposing just the attributes referenced by the
# rest of the codebase.  When the real packages are available locally, the
# `try` imports succeed and the stubs are not used.

# ----- PyMOL ---------------------------------------------------------------
try:
    import pymol  # type: ignore
except ModuleNotFoundError:
    pymol = types.ModuleType("pymol")

    def _noop(*args, **kwargs):
        """No-op replacement for PyMOL functions when library is missing."""

    pymol.finish_launching = _noop  # type: ignore

    class _CmdStub:  # noqa: D401
        def reinitialize(self, *args, **kwargs):
            pass

    pymol.cmd = _CmdStub()  # type: ignore
    sys.modules["pymol"] = pymol

# ----- Redis --------------------------------------------------------------
try:
    import redis  # type: ignore
except ModuleNotFoundError:
    redis = types.ModuleType("redis")

    class _DummyRedis:  # noqa: D401
        def __init__(self, *args, **kwargs):
            pass

        def ping(self):
            return True

    redis.Redis = _DummyRedis  # type: ignore
    sys.modules["redis"] = redis

# ----- PubChemPy ----------------------------------------------------------
try:
    import pubchempy  # type: ignore
except ModuleNotFoundError:
    pubchempy = types.ModuleType("pubchempy")

    class _DummyCompound:  # noqa: D401
        def __getattr__(self, name):
            return None

    def _unavailable(*args, **kwargs):
        raise ImportError("pubchempy is not installed in this environment.")

    pubchempy.Compound = _DummyCompound  # type: ignore
    pubchempy.get_compounds = _unavailable  # type: ignore
    sys.modules["pubchempy"] = pubchempy
