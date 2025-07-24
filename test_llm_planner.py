#!/usr/bin/env python3
"""
Test script to isolate and test the LLM call used in the graphic endpoint.
This replicates the exact same call that fails in the /graphic/make endpoint.
"""

import logging
import os
import sys

# Add the api directory to the path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))

from api.diagram.planner_llm.service import plan_simple

# Set up logging to see what's happening
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_llm_planner():
    """Test the exact same LLM call that's failing in the graphic endpoint."""

    # This is the exact same brief from the curl command
    brief = "Create an infographic about the Calvin cycle with RuBisCO, CO2, RuBP, and glucose"

    try:
        logger.info(f"Testing LLM planner with brief: {brief}")

        # This is the exact same call made in the /graphic/make endpoint
        yaml_spec = plan_simple(brief=brief, width=960, height=640)

        logger.info("LLM call succeeded!")
        logger.info(f"Generated YAML spec length: {len(yaml_spec)} characters")
        logger.info("First 500 characters of YAML spec:")
        logger.info(yaml_spec[:500])

        # Save the result to a file for inspection
        with open("test_llm_output.yaml", "w") as f:
            f.write(yaml_spec)

        logger.info("YAML spec saved to test_llm_output.yaml")
        return True

    except Exception as e:
        logger.error(f"LLM call failed with error: {str(e)}")
        logger.error(f"Error type: {type(e).__name__}")
        import traceback

        logger.error(f"Full traceback:\n{traceback.format_exc()}")
        return False


if __name__ == "__main__":
    success = test_llm_planner()
    if success:
        print("✅ LLM test passed!")
        sys.exit(0)
    else:
        print("❌ LLM test failed!")
        sys.exit(1)
