"""
Minimal entry point for Cloudflare Workers Python Runtime
This file provides a simple health check endpoint in case of import issues
"""

import os
import json
import logging
from js import Response

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Attempt to load the main app
try:
    from main import on_fetch
    logger.info("Successfully imported main application")
    HAS_MAIN_APP = True
except Exception as e:
    logger.error(f"Failed to import main application: {str(e)}")
    HAS_MAIN_APP = False

async def on_fetch(request, env=None):
    """Entry point for Cloudflare Workers"""
    # If main app failed to import, provide minimal health check
    if not HAS_MAIN_APP:
        return Response.new(
            json.dumps({
                "status": "error",
                "message": "Application not loaded correctly",
                "environment": "workers"
            }),
            {
                "status": 500,
                "headers": {"Content-Type": "application/json"}
            }
        )
    
    # Pass to main application
    try:
        from main import on_fetch as main_on_fetch
        return await main_on_fetch(request, env)
    except Exception as e:
        logger.error(f"Error in main application: {str(e)}")
        return Response.new(
            json.dumps({
                "status": "error",
                "message": str(e),
                "environment": "workers"
            }),
            {
                "status": 500,
                "headers": {"Content-Type": "application/json"}
            }
        )