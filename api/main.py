import datetime
import os
import threading
from pathlib import Path
from typing import Any, Dict

import pymol
import redis
from fastapi import Depends, FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from api import routers

# Import and initialize the model registry at startup
from api.agent_management.model_config import register_models
from api.utils.rate_limit import RateLimitMiddleware

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

# rate limiting
app.add_middleware(RateLimitMiddleware)


@app.get("/health")
async def health_check() -> Dict[str, Any]:
    """Health check endpoint."""
    try:
        # Check Redis connection
        redis_client = redis.Redis(
            host=os.environ.get("REDIS_HOST", "localhost"),
            port=int(os.environ.get("REDIS_PORT", 6379)),
        )
        redis_client.ping()

        # Check PyMOL
        with pymol_lock:
            pymol.cmd.reinitialize()

        return {
            "status": "healthy",
            "timestamp": datetime.datetime.now().isoformat(),
            "services": {"redis": "connected", "pymol": "running"},
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
