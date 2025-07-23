#!/usr/bin/env python3
"""Direct test of PyMOL templates without circular imports."""

import os
import sys

# Add the direct path to pymol_templates
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "agent_management"))


def test_pymol_templates_direct():
    """Test PyMOL templates by importing the module directly."""
    print("🧪 Testing PyMOL templates directly...")

    # Direct import to avoid circular dependencies
    import pymol_templates

    print("\n1. Testing transparent_molecule_scene...")
    commands = pymol_templates.transparent_molecule_scene("1ubq", "cartoon")
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify key transparency commands
    assert "set ray_opaque_background, 0" in commands, "Missing transparency setting"
    assert "set antialias, 1" in commands, "Missing antialiasing"
    assert "fetch 1ubq, async=0" in commands, "Missing fetch command"
    print("✅ Transparent molecule template working correctly!")

    print("\n2. Testing publication_quality_scene...")
    commands = pymol_templates.publication_quality_scene(
        "1ubq", highlight_selection="resi 50-60", transparent=True
    )
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify quality settings
    assert "set ray_trace_mode, 1" in commands, "Missing ray trace mode"
    assert "set ray_opaque_background, 0" in commands, "Missing transparency"
    assert "set antialias, 2" in commands, "Missing maximum antialiasing"
    assert "select highlight, (resi 50-60)" in commands, "Missing highlight selection"
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
    commands = pymol_templates.annotated_molecule_scene("1ubq", annotations)
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify annotations
    assert any(
        "distance bond1" in cmd for cmd in commands
    ), "Missing distance annotation"
    assert any("label resi 10" in cmd for cmd in commands), "Missing label annotation"
    assert any("angle angle1" in cmd for cmd in commands), "Missing angle annotation"
    print("✅ Annotated molecule template working correctly!")

    print("\n4. Testing transparent_binding_site_scene...")
    commands = pymol_templates.transparent_binding_site_scene(
        "1ubq", "resi 100-120", transparent_bg=True
    )
    print(f"Generated {len(commands)} commands:")
    for cmd in commands:
        print(f"  - {cmd}")

    # Verify binding site and transparency
    assert "set ray_opaque_background, 0" in commands, "Missing transparency"
    assert "show surface, binding_site around 4" in commands, "Missing surface display"
    assert (
        "select binding_site, (resi 100-120)" in commands
    ), "Missing binding site selection"
    print("✅ Transparent binding site template working correctly!")

    print("\n5. Testing existing template functions still work...")
    # Test existing functions to ensure we didn't break anything
    commands = pymol_templates.overview_scene("1ubq")
    assert "fetch 1ubq, async=0" in commands, "Overview scene broken"
    assert "show cartoon" in commands, "Overview scene missing cartoon"
    print("✅ Existing overview template still working!")

    commands = pymol_templates.binding_site_scene("1ubq", "resi 50-60")
    assert "select binding_site, (resi 50-60)" in commands, "Binding site scene broken"
    print("✅ Existing binding site template still working!")


def test_security_allowlist():
    """Test security module updates."""
    print("\n🔒 Testing security module...")

    # Add the utils path
    sys.path.append(os.path.join(os.path.dirname(__file__), "..", "utils"))
    import security

    # Test that new commands are in the whitelist
    new_commands = [
        "ray",
        "antialias",
        "depth_cue",
        "distance",
        "angle",
        "h_add",
        "find_pairs",
    ]

    for cmd in new_commands:
        assert (
            cmd in security.ALLOWED_COMMANDS
        ), f"Command '{cmd}' not in ALLOWED_COMMANDS"
        print(f"  ✅ {cmd} allowed")

    # Test validation of ray-tracing commands
    test_commands = [
        "set ray_opaque_background, 0",
        "set antialias, 1",
        "ray 1920, 1080",
        "distance bond1, resi 1 and name CA, resi 2 and name CA",
    ]

    try:
        security.validate_commands(test_commands)
        print("✅ All ray-tracing commands validated successfully!")
    except ValueError as e:
        print(f"❌ Security validation failed: {e}")
        raise


def demonstrate_transparent_rendering():
    """Demonstrate the new transparent rendering capabilities."""
    print("\n🎨 Demonstrating transparent rendering capabilities...")

    import pymol_templates

    # Example 1: Basic transparent molecule
    print("\nExample 1: Basic Transparent Molecule")
    commands = pymol_templates.transparent_molecule_scene("1ubq", "sticks")
    print("PyMOL commands for transparent molecule:")
    for cmd in commands:
        print(f"  {cmd}")

    # Example 2: Publication quality with highlighting
    print("\nExample 2: Publication Quality with Highlighting")
    commands = pymol_templates.publication_quality_scene(
        "1ubq", highlight_selection="resi 50-70", transparent=True
    )
    print("PyMOL commands for publication quality:")
    for cmd in commands:
        print(f"  {cmd}")

    # Example 3: Annotated structure
    print("\nExample 3: Annotated Structure with Measurements")
    annotations = [
        {
            "type": "distance",
            "name": "active_site_width",
            "atom1": "resi 50 and name CA",
            "atom2": "resi 60 and name CA",
        },
        {"type": "label", "selection": "resi 55", "text": "Catalytic Residue"},
    ]
    commands = pymol_templates.annotated_molecule_scene("1ubq", annotations)
    print("PyMOL commands for annotated structure:")
    for cmd in commands:
        print(f"  {cmd}")


if __name__ == "__main__":
    print("🧪 Direct PyMOL Template Testing (No Circular Imports)\n")

    try:
        test_pymol_templates_direct()
        test_security_allowlist()
        demonstrate_transparent_rendering()

        print("\n" + "=" * 70)
        print("🎉 ALL TESTS PASSED!")
        print("✅ Transparent background rendering templates working!")
        print("✅ Publication quality rendering available!")
        print("✅ Annotation and measurement support added!")
        print("✅ Security validation updated for new commands!")
        print("✅ Backward compatibility maintained!")

        print("\n📋 Features implemented:")
        print("  • Transparent background support (ray_opaque_background, 0)")
        print("  • High-quality ray-tracing with antialiasing")
        print("  • Publication-grade rendering modes")
        print("  • Molecular annotations (distances, angles, labels)")
        print("  • Advanced PyMOL controls (shadows, depth cues)")
        print("  • Security whitelist updated for new commands")

    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
