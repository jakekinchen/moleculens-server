import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient


@pytest.mark.integration
def test_transparent_background_rendering():
    """Test transparent PNG generation with alpha channel support."""

    # Import the router module directly to avoid circular imports
    routes_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../routers/render/routes.py")
    )
    spec = importlib.util.spec_from_file_location("routes", routes_path)
    routes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(routes)

    app = FastAPI()
    app.include_router(routes.router)

    with TestClient(app) as client:
        # Test transparent background request
        request_data = {
            "description": "Create a transparent background visualization of caffeine molecule",
            "format": "image",
            "transparent_background": True,
            "ray_trace": True,
            "resolution": [1280, 720],
            "dpi": 300,
            "antialias": True,
            "ray_shadow": True,
            "background_color": "white",
        }

        response = client.post("/render", json=request_data)

        # Verify response
        assert response.status_code == 200, f"Request failed: {response.text}"
        assert (
            len(response.content) > 1000
        ), f"Response too small: {len(response.content)} bytes"

        # Save test output for inspection
        test_output_dir = "/tmp/transparent_test_output"
        os.makedirs(test_output_dir, exist_ok=True)

        output_path = os.path.join(test_output_dir, "transparent_molecule.png")
        with open(output_path, "wb") as f:
            f.write(response.content)

        print(f"✅ Transparent background test passed - saved to {output_path}")
        return output_path


@pytest.mark.integration
def test_publication_quality_rendering():
    """Test high-quality ray-traced output for publications."""

    routes_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../routers/render/routes.py")
    )
    spec = importlib.util.spec_from_file_location("routes", routes_path)
    routes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(routes)

    app = FastAPI()
    app.include_router(routes.router)

    with TestClient(app) as client:
        request_data = {
            "description": "Generate publication quality image of protein 1ubq with binding site highlighted",
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
        assert len(response.content) > 50000, "High-quality image should be larger"

        # Save test output
        test_output_dir = "/tmp/transparent_test_output"
        os.makedirs(test_output_dir, exist_ok=True)

        output_path = os.path.join(test_output_dir, "publication_quality.png")
        with open(output_path, "wb") as f:
            f.write(response.content)

        print(f"✅ Publication quality test passed - saved to {output_path}")
        return output_path


@pytest.mark.integration
def test_ray_trace_modes():
    """Test different ray-tracing modes for visual variety."""

    routes_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../routers/render/routes.py")
    )
    spec = importlib.util.spec_from_file_location("routes", routes_path)
    routes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(routes)

    app = FastAPI()
    app.include_router(routes.router)

    test_output_dir = "/tmp/transparent_test_output"
    os.makedirs(test_output_dir, exist_ok=True)

    ray_trace_modes = ["default", "cartoon_outline", "bw", "poster"]

    with TestClient(app) as client:
        for mode in ray_trace_modes:
            request_data = {
                "description": f"Show protein structure 1ubq in {mode} style",
                "format": "image",
                "transparent_background": True,  # All with transparent background
                "ray_trace": True,
                "resolution": [1024, 768],
                "dpi": 300,
                "ray_trace_mode": mode,
                "antialias": True,
            }

            response = client.post("/render", json=request_data)
            assert response.status_code == 200

            output_path = os.path.join(test_output_dir, f"ray_mode_{mode}.png")
            with open(output_path, "wb") as f:
                f.write(response.content)

            print(f"✅ Ray-trace mode '{mode}' test passed - saved to {output_path}")


@pytest.mark.integration
def test_advanced_rendering_options():
    """Test various combinations of advanced rendering options."""

    routes_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../routers/render/routes.py")
    )
    spec = importlib.util.spec_from_file_location("routes", routes_path)
    routes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(routes)

    app = FastAPI()
    app.include_router(routes.router)

    test_output_dir = "/tmp/transparent_test_output"
    os.makedirs(test_output_dir, exist_ok=True)

    test_cases = [
        {
            "name": "basic_transparent",
            "description": "Simple transparent molecule",
            "transparent_background": True,
            "ray_trace": False,  # Fast rendering
            "dpi": 150,
        },
        {
            "name": "high_quality_opaque",
            "description": "High quality opaque background",
            "transparent_background": False,
            "ray_trace": True,
            "ray_shadow": True,
            "antialias": True,
            "dpi": 300,
        },
        {
            "name": "presentation_transparent",
            "description": "Presentation quality with transparency",
            "transparent_background": True,
            "ray_trace": True,
            "ray_trace_mode": "poster",
            "dpi": 200,
            "resolution": [1920, 1080],
        },
    ]

    with TestClient(app) as client:
        for case in test_cases:
            case_name = case.pop("name")

            # Default values
            request_data = {
                "description": case.get("description", "molecule visualization"),
                "format": "image",
                "transparent_background": False,
                "ray_trace": True,
                "resolution": [1280, 720],
                "dpi": 200,
                "ray_trace_mode": "default",
                "antialias": True,
                "ray_shadow": True,
                "depth_cue": True,
                **case,  # Override with test case specific values
            }

            response = client.post("/render", json=request_data)
            assert response.status_code == 200

            output_path = os.path.join(test_output_dir, f"advanced_{case_name}.png")
            with open(output_path, "wb") as f:
                f.write(response.content)

            print(
                f"✅ Advanced rendering '{case_name}' test passed - saved to {output_path}"
            )


def test_template_functions():
    """Test the new template functions directly."""
    import sys

    sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

    from agent_management import pymol_templates

    # Test transparent molecule scene
    commands = pymol_templates.transparent_molecule_scene("1ubq", "cartoon")
    assert "set ray_opaque_background, 0" in commands
    assert "set antialias, 1" in commands
    print("✅ Transparent molecule template function working")

    # Test publication quality scene
    commands = pymol_templates.publication_quality_scene(
        "1ubq", highlight_selection="resi 50-60", transparent=True
    )
    assert "set ray_trace_mode, 1" in commands
    assert "set ray_opaque_background, 0" in commands
    print("✅ Publication quality template function working")

    # Test annotated molecule scene
    annotations = [
        {
            "type": "distance",
            "name": "bond1",
            "atom1": "resi 1 and name CA",
            "atom2": "resi 2 and name CA",
        },
        {"type": "label", "selection": "resi 10", "text": "Active Site"},
    ]
    commands = pymol_templates.annotated_molecule_scene("1ubq", annotations)
    assert any("distance bond1" in cmd for cmd in commands)
    assert any("label resi 10" in cmd for cmd in commands)
    print("✅ Annotated molecule template function working")

    # Test transparent binding site scene
    commands = pymol_templates.transparent_binding_site_scene(
        "1ubq", "resi 100-120", transparent_bg=True
    )
    assert "set ray_opaque_background, 0" in commands
    assert "show surface, binding_site around 4" in commands
    print("✅ Transparent binding site template function working")


if __name__ == "__main__":
    """Run tests directly to see output."""
    print("🧪 Starting transparent background rendering tests...\n")

    try:
        print("=" * 60)
        print("Testing template functions...")
        test_template_functions()

        print("\n" + "=" * 60)
        print("Testing transparent background rendering...")
        test_transparent_background_rendering()

        print("\n" + "=" * 60)
        print("Testing publication quality rendering...")
        test_publication_quality_rendering()

        print("\n" + "=" * 60)
        print("Testing ray-trace modes...")
        test_ray_trace_modes()

        print("\n" + "=" * 60)
        print("Testing advanced rendering options...")
        test_advanced_rendering_options()

        print("\n" + "=" * 60)
        print(
            "🎉 ALL TESTS PASSED! Check /tmp/transparent_test_output/ for generated images."
        )

    except Exception as e:
        print(f"❌ Test failed: {str(e)}")
        import traceback

        traceback.print_exc()
