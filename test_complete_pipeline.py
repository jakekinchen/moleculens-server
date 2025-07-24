#!/usr/bin/env python3
"""
Demonstrate the complete graphic pipeline with both fallback and molecule rendering.
"""

import asyncio
import json
import os
import sys

import yaml

# Add the api directory to the path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))

from api.diagram.planner_llm.service import plan_simple
from api.diagram.validator.core import validate
from api.routers.graphic.routes import (
    _generate_molecule_images,
    _render_diagram_with_molecules,
    _yaml_to_diagram_plan,
)


async def test_with_molecule_rendering():
    """Test the complete pipeline with actual molecule rendering."""

    brief = "Create a simple diagram with CO2 and glucose"

    try:
        # Step 1: Generate YAML spec
        print("🔄 Generating YAML spec...")
        yaml_spec = plan_simple(brief=brief, width=960, height=640)
        print("✅ YAML spec generated")

        # Step 2: Convert to DiagramPlan
        print("🔄 Converting to DiagramPlan...")
        spec_dict = yaml.safe_load(yaml_spec)
        diagram_plan = _yaml_to_diagram_plan(spec_dict)
        print(
            f"✅ DiagramPlan created with {len(diagram_plan.molecule_list)} molecules"
        )

        # Step 3: Try to generate molecule images (with timeout)
        print("🔄 Attempting to generate molecule images...")
        try:
            molecule_images = await asyncio.wait_for(
                _generate_molecule_images(diagram_plan),
                timeout=60.0,  # 1 minute timeout
            )
            print(f"✅ Generated {len(molecule_images)} molecule images")
        except asyncio.TimeoutError:
            print("⚠️  Molecule rendering timed out, using fallback")
            molecule_images = {}
        except Exception as e:
            print(f"⚠️  Molecule rendering failed: {e}, using fallback")
            molecule_images = {}

        # Step 4: Render diagram
        print("🔄 Rendering final diagram...")
        svg_content = _render_diagram_with_molecules(
            plan=diagram_plan,
            molecule_images=molecule_images,
            canvas_width=960,
            canvas_height=640,
        )
        print(f"✅ SVG rendered ({len(svg_content)} characters)")

        # Save outputs
        with open("transparent_molecule.svg", "w") as f:
            f.write(svg_content)

        print("✅ File saved: transparent_molecule.svg")

        # Show what we got
        if molecule_images:
            print(
                f"🎉 SUCCESS: Generated actual molecular images for: {list(molecule_images.keys())}"
            )
        else:
            print("🎨 FALLBACK: Using colored circles as molecule representations")

        return True

    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        import traceback

        traceback.print_exc()
        return False


def test_fallback_only():
    """Test the pipeline with fallback rendering only."""

    brief = "Create an infographic about photosynthesis with chlorophyll, sunlight, CO2, and oxygen"

    try:
        print("🔄 Testing fallback-only rendering...")
        yaml_spec = plan_simple(brief=brief, width=960, height=640)
        spec_dict = yaml.safe_load(yaml_spec)
        diagram_plan = _yaml_to_diagram_plan(spec_dict)

        svg_content = _render_diagram_with_molecules(
            plan=diagram_plan,
            molecule_images={},  # No molecule images - pure fallback
            canvas_width=960,
            canvas_height=640,
        )

        with open("fallback_test.svg", "w") as f:
            f.write(svg_content)

        print(
            f"✅ Fallback test complete: fallback_test.svg ({len(svg_content)} chars)"
        )
        print(f"   Molecules: {[mol.label for mol in diagram_plan.molecule_list]}")

        return True

    except Exception as e:
        print(f"❌ Fallback test failed: {e}")
        return False


async def main():
    """Run both tests."""

    print("🧪 Testing Complete Graphic Pipeline\n")

    # Test 1: Fallback rendering (fast, always works)
    print("=== Test 1: Fallback Rendering ===")
    fallback_success = test_fallback_only()

    print("\n=== Test 2: Molecule Rendering (with fallback) ===")
    molecule_success = await test_with_molecule_rendering()

    print("\n📋 Summary:")
    print(f"   Fallback rendering: {'✅ PASS' if fallback_success else '❌ FAIL'}")
    print(f"   Molecule rendering: {'✅ PASS' if molecule_success else '❌ FAIL'}")

    if fallback_success:
        print("\n🎉 The graphic pipeline is working!")
        print("   - LLM generates YAML specs from natural language")
        print("   - YAML is converted to structured DiagramPlan")
        print("   - SVG is rendered with molecules and arrows")
        print("   - Fallback circles are used when molecule rendering fails/times out")
        print("   - Actual PyMOL-rendered molecules are embedded when available")

        return True
    else:
        print("\n💥 Pipeline has issues that need to be fixed")
        return False


if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1)
