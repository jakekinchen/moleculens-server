#!/usr/bin/env python3
"""
Test script to demonstrate the diagram data model working correctly.
This shows how the YAML spec is converted to DiagramPlan and then rendered.
"""

import json
import os
import sys

import yaml

# Add the api directory to the path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))

from api.diagram.models import DiagramPlan
from api.diagram.renderer.diagram import render_diagram
from api.routers.graphic.routes import _yaml_to_diagram_plan


def test_data_model():
    """Test the data model conversion and rendering."""

    # Load the successful YAML spec from our test
    with open("test_llm_output.yaml", "r") as f:
        yaml_content = f.read()

    print("✅ Loaded YAML spec from test_llm_output.yaml")

    # Parse YAML
    spec_dict = yaml.safe_load(yaml_content)
    print("✅ Parsed YAML successfully")

    # Convert to DiagramPlan
    try:
        diagram_plan = _yaml_to_diagram_plan(spec_dict)
        print("✅ Converted YAML to DiagramPlan successfully")
        print(f"   Plan: {diagram_plan.plan}")
        print(f"   Molecules: {len(diagram_plan.molecule_list)}")
        print(f"   Arrows: {len(diagram_plan.arrows)}")
        print(f"   Canvas: {diagram_plan.canvas_width}x{diagram_plan.canvas_height}")
    except Exception as e:
        print(f"❌ Failed to convert YAML to DiagramPlan: {e}")
        return False

    # Render the diagram
    try:
        svg_content = render_diagram(
            plan=diagram_plan,
            molecule_data={},
            canvas_width=diagram_plan.canvas_width,
            canvas_height=diagram_plan.canvas_height,
        )
        print("✅ Rendered diagram successfully")
        print(f"   SVG length: {len(svg_content)} characters")

        # Save the SVG for inspection
        with open("test_diagram_output.svg", "w") as f:
            f.write(svg_content)
        print("✅ Saved SVG to test_diagram_output.svg")

    except Exception as e:
        print(f"❌ Failed to render diagram: {e}")
        return False

    # Test the molecule placements
    print("\n📍 Molecule Placements:")
    for i, mol in enumerate(diagram_plan.molecule_list):
        print(f"   {i+1}. {mol.label} at ({mol.x:.1f}, {mol.y:.1f})")

    # Test the arrows
    print("\n🏹 Arrows:")
    for i, arrow in enumerate(diagram_plan.arrows):
        print(f"   {i+1}. {arrow.start} → {arrow.end}: '{arrow.text}'")

    return True


if __name__ == "__main__":
    success = test_data_model()
    if success:
        print("\n🎉 Data model test passed! The diagram pipeline is working correctly.")
        sys.exit(0)
    else:
        print("\n💥 Data model test failed!")
        sys.exit(1)
