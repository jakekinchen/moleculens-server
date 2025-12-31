"""FastAPI application entry point."""

import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from moleculens import __version__
from moleculens.api.routes import router
from moleculens.core import get_logger, settings, setup_logging
from moleculens.db import init_db

logger = get_logger(__name__)


def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    setup_logging(settings.log_level)

    app = FastAPI(
        title="Moleculens API",
        description="Psi4 orbital/density overlays for Three.js",
        version=__version__,
        docs_url="/docs",
        redoc_url="/redoc",
    )

    # Configure CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.cors_origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Include routes
    app.include_router(router)

    @app.on_event("startup")
    async def startup_event() -> None:
        """Initialize database on startup."""
        logger.info(
            "Starting Moleculens API",
            version=__version__,
            log_level=settings.log_level,
        )
        init_db()

    return app


# Create app instance for uvicorn
app = create_app()


def run() -> None:
    """Run the API server (entry point for moleculens-api command)."""
    uvicorn.run(
        "moleculens.api.main:app",
        host="0.0.0.0",
        port=8000,
        reload=False,
        log_level=settings.log_level.lower(),
    )


if __name__ == "__main__":
    run()
