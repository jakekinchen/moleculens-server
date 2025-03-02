from fastapi import Depends, FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import routers
import os

# Import and initialize the model registry at startup
from agent_management.model_config import register_models

# Register all models
register_models()

app = FastAPI(
    title="AI Backend",
    version="0.0.1",
)
app.mount("/static", StaticFiles(directory="static"), name="static")

# Determine if we're in development or production mode
is_development = os.environ.get("ENVIRONMENT", "development").lower() == "development"

# Configure CORS based on environment
if is_development:
    # Development mode - use specific origins
    origins = [
        "http://localhost:3000",  # React development server
        "http://localhost:8000",  # Backend server (for serving frontend in production)
    ]
    allow_credentials = False
else:
    # Production mode - use wildcard
    origins = ["*"]
    allow_credentials = False  # Cannot use credentials with wildcard origin

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=allow_credentials,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
)

# user management related endpoints
app.include_router(routers.prompt.router)
app.include_router(routers.geometry.router)
