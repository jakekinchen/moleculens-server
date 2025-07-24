#!/usr/bin/env python3
"""
Final demonstration of the complete PNG pipeline working end-to-end.
"""

import asyncio
import base64
import json
import os
import sys
from pathlib import Path

import yaml

# Add the api directory to the path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))

from api.diagram.planner_llm.service import plan_simple
from api.routers.graphic.routes import (
    _render_diagram_with_molecules,
    _svg_to_png,
    _yaml_to_diagram_plan,
)


def test_pipeline_without_molecules():
    """Test the pipeline with fallback molecule rendering (no network calls)."""

    brief = (
        "Create a diagram showing photosynthesis with CO2, water, glucose, and oxygen"
    )

    try:
        print("🚀 Testing Complete Pipeline (Fallback Mode)")
        print("=" * 50)

        # Step 1: Generate YAML spec from natural language
        print("🔄 Step 1: Generating YAML spec...")
        yaml_spec = plan_simple(brief=brief, width=800, height=600)
        print(f"✅ YAML generated ({len(yaml_spec)} chars)")

        # Step 2: Convert YAML to DiagramPlan
        print("🔄 Step 2: Converting to DiagramPlan...")
        spec_dict = yaml.safe_load(yaml_spec)
        diagram_plan = _yaml_to_diagram_plan(spec_dict)
        molecules = [mol.molecule for mol in diagram_plan.molecule_list]
        print(f"✅ DiagramPlan created with molecules: {molecules}")

        # Step 3: Render SVG with fallback circles (no molecule images)
        print("🔄 Step 3: Rendering SVG with fallback molecules...")
        svg_content = _render_diagram_with_molecules(
            plan=diagram_plan,
            molecule_images={},  # Empty - will use fallback circles
            canvas_width=800,
            canvas_height=600,
        )
        print(f"✅ SVG rendered ({len(svg_content)} chars)")

        # Step 4: Convert SVG to PNG
        print("🔄 Step 4: Converting SVG to PNG...")
        png_base64 = _svg_to_png(svg_content, 800, 600)

        if png_base64:
            # Save final PNG
            png_data = base64.b64decode(png_base64)
            with open("final_pipeline_demo.png", "wb") as f:
                f.write(png_data)
            print(f"✅ PNG created ({len(png_data)} bytes)")
            print("   💾 Saved: final_pipeline_demo.png")

            # Save SVG for comparison
            with open("final_pipeline_demo.svg", "w") as f:
                f.write(svg_content)
            print("   💾 Saved: final_pipeline_demo.svg")

            # Save YAML for inspection
            with open("final_pipeline_demo.yaml", "w") as f:
                f.write(yaml_spec)
            print("   💾 Saved: final_pipeline_demo.yaml")

            return True
        else:
            print("❌ PNG conversion failed")
            return False

    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        import traceback

        traceback.print_exc()
        return False


def main():
    """Run the demonstration."""

    print("🧪 Final Pipeline Demonstration")
    print("=" * 60)
    print("This test demonstrates the complete working pipeline:")
    print("  1. Natural language → YAML specification")
    print("  2. YAML → structured diagram plan")
    print("  3. Diagram plan → SVG with molecules")
    print("  4. SVG → PNG conversion")
    print()

    success = test_pipeline_without_molecules()

    print("\n📊 Results:")
    print("=" * 30)

    if success:
        print("🎉 COMPLETE PIPELINE SUCCESS!")
        print()
        print("✨ What was accomplished:")
        print("   ✅ Natural language processing (LLM)")
        print("   ✅ YAML specification generation")
        print("   ✅ Structured diagram planning")
        print("   ✅ SVG rendering with molecules")
        print("   ✅ PNG conversion (fallback method)")
        print()
        print("📁 Files created:")
        print("   - final_pipeline_demo.yaml (YAML spec)")
        print("   - final_pipeline_demo.svg (SVG diagram)")
        print("   - final_pipeline_demo.png (Final PNG output)")
        print()
        print("🔧 Next steps for production:")
        print("   - Install proper SVG→PNG converter (cairosvg + Cairo)")
        print("   - Enable /render endpoint for transparent molecules")
        print("   - Add caching and error handling")
        print("   - Implement rate limiting")

        return True
    else:
        print("❌ Pipeline has issues that need fixing")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
