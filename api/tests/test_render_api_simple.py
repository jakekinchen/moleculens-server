#!/usr/bin/env python3
"""Test the enhanced render API with transparent background support."""

import json
import os
import sys
from pathlib import Path

from fastapi import FastAPI
from fastapi.testclient import TestClient

# Add path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))


def test_render_api_with_transparency():
    """Test the render API with new transparent background options."""
    print("🧪 Testing Enhanced Render API with Transparency...\n")

    # Create a minimal FastAPI app with stubbed PyMOL
    app = FastAPI()

    # Mock the render router with stubbed dependencies
    import tempfile
    from unittest.mock import MagicMock, patch

    # Create a temporary directory for outputs
    temp_dir = tempfile.mkdtemp()
    print(f"Using temp directory: {temp_dir}")

    # Import and setup the router with mocks
    try:
        # Mock PyMOL before importing the routes
        with patch.dict(
            "sys.modules",
            {
                "pymol": MagicMock(),
            },
        ):
            import routers.render.routes as render_routes

            # Mock the pymol_cmd object
            mock_pymol_cmd = MagicMock()
            render_routes.pymol_cmd = mock_pymol_cmd

            # Mock the translator
            def mock_translate(description):
                return [
                    "fetch 1ubq, async=0",
                    "hide everything",
                    "show cartoon",
                    "color chainbow, all",
                ]

            with patch("routers.render.routes.pymol_translator") as mock_translator:
                mock_translator.translate = mock_translate

                # Mock cache and security
                with patch("routers.render.routes.cache") as mock_cache, patch(
                    "routers.render.routes.security"
                ) as mock_security:

                    mock_cache.get.return_value = None  # No cached result
                    mock_cache.set.return_value = None
                    mock_cache.CACHE_DIR = Path(temp_dir)
                    mock_security.validate_commands.return_value = None

                    # Include the router
                    app.include_router(render_routes.router)

                    # Test with transparent background
                    with TestClient(app) as client:
                        print("1. Testing basic transparent background rendering...")

                        request_data = {
                            "description": "transparent background molecule visualization",
                            "format": "image",
                            "transparent_background": True,
                            "ray_trace": True,
                            "resolution": [1280, 720],
                            "dpi": 300,
                            "antialias": True,
                            "ray_shadow": True,
                        }

                        response = client.post("/render", json=request_data)

                        print(f"Response status: {response.status_code}")
                        if response.status_code != 200:
                            print(f"Response text: {response.text}")

                        assert (
                            response.status_code == 200
                        ), f"Request failed: {response.text}"

                        # Verify PyMOL commands were called with transparency settings
                        mock_pymol_cmd.set.assert_any_call("ray_opaque_background", 0)
                        mock_pymol_cmd.set.assert_any_call("antialias", 1)
                        mock_pymol_cmd.set.assert_any_call("ray_shadow", 1)
                        mock_pymol_cmd.ray.assert_called_with(1280, 720)

                        print("✅ Transparent background rendering working!")

                        print("\n2. Testing publication quality rendering...")

                        request_data = {
                            "description": "publication quality protein structure",
                            "format": "image",
                            "transparent_background": False,
                            "ray_trace": True,
                            "resolution": [2560, 1440],
                            "dpi": 400,
                            "ray_trace_mode": "poster",
                            "antialias": True,
                            "ray_shadow": True,
                            "depth_cue": True,
                        }

                        response = client.post("/render", json=request_data)
                        assert response.status_code == 200

                        # Verify high-quality settings
                        mock_pymol_cmd.set.assert_any_call(
                            "ray_trace_mode", 1
                        )  # poster mode
                        mock_pymol_cmd.ray.assert_called_with(2560, 1440)

                        print("✅ Publication quality rendering working!")

                        print("\n3. Testing different ray-trace modes...")

                        modes = {"cartoon_outline": 3, "bw": 2, "poster": 1}

                        for mode_name, mode_value in modes.items():
                            request_data = {
                                "description": f"test {mode_name} mode",
                                "ray_trace_mode": mode_name,
                                "ray_trace": True,
                            }

                            response = client.post("/render", json=request_data)
                            assert response.status_code == 200

                            mock_pymol_cmd.set.assert_any_call(
                                "ray_trace_mode", mode_value
                            )
                            print(f"  ✅ {mode_name} mode working!")

                        print("\n4. Testing background color setting...")

                        request_data = {
                            "description": "test background color",
                            "background_color": "black",
                            "transparent_background": False,
                        }

                        response = client.post("/render", json=request_data)
                        assert response.status_code == 200

                        mock_pymol_cmd.do.assert_any_call("bg_color black")
                        print("✅ Background color setting working!")

    except Exception as e:
        print(f"❌ Test failed: {str(e)}")
        import traceback

        traceback.print_exc()
        raise


def test_request_model_validation():
    """Test that the enhanced RenderRequest model validates correctly."""
    print("\n🔍 Testing RenderRequest model validation...")

    from routers.render.routes import RenderRequest

    # Test default values
    request = RenderRequest(description="test molecule")
    assert request.transparent_background is False
    assert request.ray_trace is True
    assert request.resolution == (1920, 1080)
    assert request.dpi == 300
    assert request.ray_trace_mode == "default"
    print("✅ Default values correct!")

    # Test custom values
    request = RenderRequest(
        description="custom molecule",
        transparent_background=True,
        ray_trace=True,
        resolution=(2560, 1440),
        dpi=400,
        ray_trace_mode="poster",
        antialias=True,
        ray_shadow=True,
        depth_cue=True,
        background_color="black",
    )

    assert request.transparent_background is True
    assert request.resolution == (2560, 1440)
    assert request.dpi == 400
    assert request.ray_trace_mode == "poster"
    assert request.background_color == "black"
    print("✅ Custom values validated correctly!")

    # Test validation of ray_trace_mode enum
    try:
        RenderRequest(description="test", ray_trace_mode="invalid_mode")
        assert False, "Should have failed validation"
    except ValueError:
        print("✅ Enum validation working!")


def demonstrate_api_usage():
    """Demonstrate how to use the new API features."""
    print("\n📚 API Usage Examples:")

    examples = [
        {
            "name": "Basic Transparent Background",
            "request": {
                "description": "show caffeine molecule",
                "transparent_background": True,
                "ray_trace": True,
                "dpi": 300,
            },
        },
        {
            "name": "Publication Quality",
            "request": {
                "description": "protein binding site analysis",
                "ray_trace": True,
                "resolution": [2560, 1440],
                "dpi": 400,
                "ray_trace_mode": "poster",
                "antialias": True,
                "ray_shadow": True,
                "depth_cue": True,
            },
        },
        {
            "name": "Presentation Mode",
            "request": {
                "description": "enzyme mechanism overview",
                "transparent_background": True,
                "ray_trace": True,
                "ray_trace_mode": "cartoon_outline",
                "resolution": [1920, 1080],
                "dpi": 200,
                "background_color": "white",
            },
        },
        {
            "name": "Fast Preview",
            "request": {
                "description": "quick molecule preview",
                "transparent_background": True,
                "ray_trace": False,  # Fast rendering
                "dpi": 150,
            },
        },
    ]

    for example in examples:
        print(f"\n{example['name']}:")
        print(f"  POST /render")
        print(f"  {json.dumps(example['request'], indent=2)}")


if __name__ == "__main__":
    print("🧪 Enhanced Render API Testing\n")

    try:
        test_request_model_validation()
        test_render_api_with_transparency()
        demonstrate_api_usage()

        print("\n" + "=" * 70)
        print("🎉 ALL API TESTS PASSED!")
        print("✅ Enhanced RenderRequest model working!")
        print("✅ Transparent background rendering functional!")
        print("✅ Ray-tracing modes implemented!")
        print("✅ Publication quality settings available!")
        print("✅ Background color control working!")

    except Exception as e:
        print(f"\n❌ API test failed: {str(e)}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
