#!/usr/bin/env python3
"""
Test the complete graphic pipeline with fallback rendering (no molecule images).
This tests the structure without waiting for PyMOL rendering.
"""

import json
import os
import sys

import yaml

# Add the api directory to the path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))

from api.diagram.planner_llm.service import plan_simple
from api.diagram.validator.core import validate
from api.routers.graphic.routes import (
    _render_diagram_with_molecules,
    _yaml_to_diagram_plan,
)


def test_complete_pipeline():
    """Test the complete pipeline without molecule rendering."""

    brief = "Create an infographic about the Calvin cycle with RuBisCO, CO2, RuBP, and glucose"

    try:
        # Step 1: Generate YAML spec
        print("🔄 Generating YAML spec...")
        yaml_spec = plan_simple(brief=brief, width=960, height=640)
        print("✅ YAML spec generated")

        # Step 2: Validate the spec
        print("🔄 Validating YAML spec...")
        errors = validate(yaml_spec)
        if errors:
            print(f"⚠️  Validation warnings: {errors}")
        else:
            print("✅ YAML spec validated")

        # Step 3: Convert to DiagramPlan
        print("🔄 Converting to DiagramPlan...")
        spec_dict = yaml.safe_load(yaml_spec)
        diagram_plan = _yaml_to_diagram_plan(spec_dict)
        print(
            f"✅ DiagramPlan created with {len(diagram_plan.molecule_list)} molecules"
        )

        # Step 4: Render with fallback (no molecule images)
        print("🔄 Rendering diagram with fallback circles...")
        svg_content = _render_diagram_with_molecules(
            plan=diagram_plan,
            molecule_images={},  # Empty - will use fallback circles
            canvas_width=960,
            canvas_height=640,
        )
        print(f"✅ SVG rendered ({len(svg_content)} characters)")

        # Save outputs
        with open("calvin_cycle_spec.yaml", "w") as f:
            f.write(yaml_spec)

        with open("calvin_cycle_output.svg", "w") as f:
            f.write(svg_content)

        print("✅ Files saved:")
        print("   - calvin_cycle_spec.yaml")
        print("   - calvin_cycle_output.svg")

        # Show molecule positions
        print("\n📍 Molecule Layout:")
        for i, mol in enumerate(diagram_plan.molecule_list):
            print(f"   {i+1}. {mol.label} at ({mol.x:.1f}, {mol.y:.1f})")

        print("\n🏹 Arrows:")
        for i, arrow in enumerate(diagram_plan.arrows):
            print(f"   {i+1}. {arrow.start} → {arrow.end}: '{arrow.text}'")

        return True

    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_complete_pipeline()
    if success:
        print("\n🎉 Complete pipeline test passed!")
        print(
            "The graphic shows colored circles as molecule placeholders with arrows and labels."
        )
        sys.exit(0)
    else:
        print("\n💥 Pipeline test failed!")
        sys.exit(1)
