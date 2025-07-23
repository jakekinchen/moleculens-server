#!/usr/bin/env python3
"""Simple test for the new PyMOL template functions."""

import os
import sys

# Add api directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))


def test_template_functions():
    """Test the new template functions directly."""
    print("🧪 Testing PyMOL template functions...")

    # Import the templates module directly
    from agent_management.pymol_templates import (
        annotated_molecule_scene,
        publication_quality_scene,
        transparent_binding_site_scene,
        transparent_molecule_scene,
    )

    print("\n1. Testing transparent_molecule_scene...")
    commands = transparent_molecule_scene("1ubq", "cartoon")
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify key transparency commands
    assert "set ray_opaque_background, 0" in commands
    assert "set antialias, 1" in commands
    print("✅ Transparent molecule template working correctly!")

    print("\n2. Testing publication_quality_scene...")
    commands = publication_quality_scene(
        "1ubq", highlight_selection="resi 50-60", transparent=True
    )
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify quality settings
    assert "set ray_trace_mode, 1" in commands
    assert "set ray_opaque_background, 0" in commands
    assert "set antialias, 2" in commands
    print("✅ Publication quality template working correctly!")

    print("\n3. Testing annotated_molecule_scene...")
    annotations = [
        {
            "type": "distance",
            "name": "bond1",
            "atom1": "resi 1 and name CA",
            "atom2": "resi 2 and name CA",
        },
        {"type": "label", "selection": "resi 10", "text": "Active Site"},
        {
            "type": "angle",
            "name": "angle1",
            "atom1": "resi 1 and name CA",
            "atom2": "resi 2 and name CA",
            "atom3": "resi 3 and name CA",
        },
    ]
    commands = annotated_molecule_scene("1ubq", annotations)
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify annotations
    assert any("distance bond1" in cmd for cmd in commands)
    assert any("label resi 10" in cmd for cmd in commands)
    assert any("angle angle1" in cmd for cmd in commands)
    print("✅ Annotated molecule template working correctly!")

    print("\n4. Testing transparent_binding_site_scene...")
    commands = transparent_binding_site_scene(
        "1ubq", "resi 100-120", transparent_bg=True
    )
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify binding site and transparency
    assert "set ray_opaque_background, 0" in commands
    assert "show surface, binding_site around 4" in commands
    assert "select binding_site, (resi 100-120)" in commands
    print("✅ Transparent binding site template working correctly!")


def test_security_updates():
    """Test that security module allows new commands."""
    print("\n🔒 Testing security module updates...")

    from utils.security import ALLOWED_COMMANDS, validate_commands

    # Test that new commands are in the whitelist
    new_commands = ["ray", "antialias", "depth_cue", "distance", "angle", "h_add"]

    for cmd in new_commands:
        assert cmd in ALLOWED_COMMANDS, f"Command '{cmd}' not in ALLOWED_COMMANDS"
        print(f"  ✅ {cmd} allowed")

    # Test validation of ray-tracing commands
    test_commands = [
        "set ray_opaque_background, 0",
        "set antialias, 1",
        "ray 1920, 1080",
        "distance bond1, resi 1 and name CA, resi 2 and name CA",
    ]

    try:
        validate_commands(test_commands)
        print("✅ All ray-tracing commands validated successfully!")
    except ValueError as e:
        print(f"❌ Security validation failed: {e}")
        raise


def test_scene_spec_updates():
    """Test SceneSpec model with new operations."""
    print("\n📋 Testing SceneSpec model updates...")

    from agent_management.scene_spec import RenderingOptions, SceneSpec

    # Test new operation types
    new_ops = [
        "transparent_molecule",
        "publication_quality",
        "annotated_molecule",
        "transparent_binding_site",
    ]

    for op in new_ops:
        try:
            SceneSpec(
                op=op,
                structure_id="1ubq",
                rendering=RenderingOptions(
                    transparent_background=True, ray_trace=True, dpi=300
                ),
            )
            print(f"  ✅ {op} operation spec created successfully")
        except Exception as e:
            print(f"  ❌ Failed to create {op} spec: {e}")
            raise

    # Test RenderingOptions
    rendering_opts = RenderingOptions(
        transparent_background=True,
        ray_trace=True,
        resolution=(2560, 1440),
        dpi=400,
        ray_trace_mode="poster",
        antialias=True,
        ray_shadow=True,
        depth_cue=True,
    )
    print(f"✅ RenderingOptions created: {rendering_opts}")


if __name__ == "__main__":
    print("🧪 Starting PyMOL Template and Configuration Tests...\n")

    try:
        test_template_functions()
        test_security_updates()
        test_scene_spec_updates()

        print("\n" + "=" * 60)
        print("🎉 ALL TEMPLATE TESTS PASSED!")
        print("✅ Transparent background rendering is ready!")
        print("✅ New PyMOL templates are working correctly!")
        print("✅ Security validation updated!")
        print("✅ Scene specifications enhanced!")

    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
