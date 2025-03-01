# Standard library imports
import os
import sys
import logging
import json
from typing import Dict, Any, List, Tuple, Optional, Union, cast

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Determine environment
is_workers_env = False
try:
    from js import Response as JsResponse
    from js import URL, Request as JsRequest
    is_workers_env = True
    logger.info("Running in Cloudflare Workers environment")
except ImportError:
    logger.info("Running in local development environment")

# Conditional imports for local development
if not is_workers_env:
    try:
        from dotenv import load_dotenv
        load_dotenv()
    except ImportError:
        logger.warning("python-dotenv not available, skipping .env loading")

# Core FastAPI imports - with fallbacks
try:
    from fastapi import FastAPI, Depends
    from fastapi.middleware.cors import CORSMiddleware
    if not is_workers_env:
        from fastapi.staticfiles import StaticFiles
except ImportError:
    logger.error("FastAPI not available")
    # Provide minimal fallbacks for Workers environment
    class FastAPI:
        def __init__(self, **kwargs):
            self.routes = []
            self.middleware = []
        def add_middleware(self, *args, **kwargs): 
            pass
        def include_router(self, router): 
            self.routes.extend(getattr(router, 'routes', []))
    class CORSMiddleware:
        def __init__(self, **kwargs): 
            pass

# Create the FastAPI app
app = FastAPI(
    title="Sci-Vis AI Backend",
    version="0.1.0",
    docs_url="/docs" if not is_workers_env else None,
)

# Only mount static files when not in Workers environment and StaticFiles is available
if not is_workers_env and 'StaticFiles' in locals():
    app.mount("/static", StaticFiles(directory="static"), name="static")

# Configure CORS
origins = ["*"]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Import routers - try/except for each import
try:
    import routers
    app.include_router(routers.prompt.router)
    app.include_router(routers.geometry.router)
except ImportError as e:
    logger.error(f"Error importing routers: {str(e)}")

# Simple response for fallback if FastAPI isn't available
def simple_json_response(data: Dict[str, Any], status: int = 200) -> Dict[str, Any]:
    return {
        "body": json.dumps(data).encode(),
        "status": status,
        "headers": {"Content-Type": "application/json"}
    }

# Cloudflare Workers entry point
async def on_fetch(request, env=None):
    """Entry point for Cloudflare Workers"""
    try:
        # Set environment variables from Workers env
        if env:
            for key in dir(env):
                # Skip methods and private attributes
                if not key.startswith('_') and not callable(getattr(env, key)):
                    os.environ[key] = str(getattr(env, key))
        
        # Get request information
        method = request.method
        url = URL.new(request.url)
        path = url.pathname
        
        # Log incoming request
        logger.info(f"Received {method} request to {path}")
        
        # Basic endpoint handling if FastAPI isn't fully available
        if path == "/" or path == "/health":
            return JsResponse.new(
                json.dumps({"status": "ok", "environment": "workers"}),
                {"status": 200, "headers": {"Content-Type": "application/json"}}
            )
            
        try:
            # Convert headers to dict
            headers = dict(request.headers)
            
            # Create ASGI scope
            scope = {
                "type": "http",
                "asgi": {"version": "3.0", "spec_version": "2.0"},
                "http_version": "1.1",
                "method": method,
                "scheme": url.protocol.replace(":", ""),
                "path": path,
                "raw_path": path.encode(),
                "query_string": url.search.encode(),
                "headers": [(k.lower().encode(), v.encode()) for k, v in headers.items()],
                "client": ("0.0.0.0", 0),
                "server": (url.hostname, int(url.port) if url.port else 443),
            }
            
            # Process request through FastAPI
            body = await request.text() if method in ["POST", "PUT", "PATCH"] else ""
            response = await app(scope, lambda: [body.encode()], lambda x: None)
            
            # Extract response details
            status = response.status_code
            response_headers = dict(response.headers.items())
            response_body = response.body
            
            # Create and return JS Response
            return JsResponse.new(
                response_body,
                {
                    "status": status,
                    "headers": response_headers
                }
            )
        except Exception as e:
            logger.error(f"Error processing request with FastAPI: {str(e)}")
            return JsResponse.new(
                json.dumps({"error": str(e), "path": path}),
                {
                    "status": 500,
                    "headers": {"Content-Type": "application/json"}
                }
            )
            
    except Exception as e:
        # Return error response for any top-level errors
        logger.error(f"Critical error processing request: {str(e)}")
        return JsResponse.new(
            json.dumps({"critical_error": str(e)}),
            {
                "status": 500,
                "headers": {"Content-Type": "application/json"}
            }
        )
