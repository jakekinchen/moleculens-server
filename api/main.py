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
    allow_credentials = True
else:
    # Production mode - include specific origins instead of wildcard
    origins = [
        "https://sci-viz-ai.vercel.app",  # Vercel deployment URL
        "https://moleculens.com",         # our vercel app
        "https://www.moleculens.com",     # our www subdomain
        "https://meshmo.com",             # meshmo.com domain (primary)
        "https://www.meshmo.com"          # www.meshmo.com domain
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

# user management related endpoints
app.include_router(routers.prompt.router)
app.include_router(routers.geometry.router)
