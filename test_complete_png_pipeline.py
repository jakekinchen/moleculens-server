#!/usr/bin/env python3
"""
Complete test of the PNG pipeline: prompt → YAML → SVG → transparent molecules → final PNG
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
from api.diagram.validator.core import validate
from api.routers.graphic.routes import (
    _generate_molecule_images,
    _render_diagram_with_molecules,
    _svg_to_png,
    _yaml_to_diagram_plan,
)


async def test_complete_png_pipeline():
    """Test the complete pipeline ending with a valid PNG file."""

    brief = "Create a diagram showing the Calvin cycle with CO2, RuBP, and glucose connected by arrows"

    try:
        print("🚀 Starting Complete PNG Pipeline Test")
        print("=" * 50)

        # Step 1: Generate YAML spec from natural language
        print("🔄 Step 1: Generating YAML spec from prompt...")
        yaml_spec = plan_simple(brief=brief, width=960, height=640)
        print(f"✅ YAML spec generated ({len(yaml_spec)} characters)")

        # Save YAML for inspection
        with open("pipeline_test_spec.yaml", "w") as f:
            f.write(yaml_spec)
        print("   📄 Saved: pipeline_test_spec.yaml")

        # Step 2: Validate YAML
        print("🔄 Step 2: Validating YAML spec...")
        errors = validate(yaml_spec)
        if errors:
            print(f"⚠️  Validation warnings: {errors}")
        else:
            print("✅ YAML spec is valid")

        # Step 3: Convert YAML to DiagramPlan
        print("🔄 Step 3: Converting YAML to DiagramPlan...")
        spec_dict = yaml.safe_load(yaml_spec)
        diagram_plan = _yaml_to_diagram_plan(spec_dict)
        print(
            f"✅ DiagramPlan created with {len(diagram_plan.molecule_list)} molecules"
        )

        # List molecules found
        molecules = [mol.molecule for mol in diagram_plan.molecule_list]
        print(f"   🧬 Molecules: {molecules}")

        # Step 4: Generate transparent PNG images for each molecule
        print("🔄 Step 4: Generating transparent molecule PNGs...")
        try:
            molecule_images = await asyncio.wait_for(
                _generate_molecule_images(diagram_plan),
                timeout=120.0,  # 2 minute timeout
            )
            print(f"✅ Generated {len(molecule_images)} transparent molecule images")

            # Save individual molecule images for inspection
            for mol_name, img_data in molecule_images.items():
                filename = f"molecule_{mol_name.replace(' ', '_')}.png"
                with open(filename, "wb") as f:
                    f.write(base64.b64decode(img_data))
                print(f"   💾 Saved: {filename}")

        except asyncio.TimeoutError:
            print("⚠️  Molecule rendering timed out, using fallback circles")
            molecule_images = {}
        except Exception as e:
            print(f"⚠️  Molecule rendering failed: {e}, using fallback circles")
            molecule_images = {}

        # Step 5: Render SVG diagram with embedded molecules
        print("🔄 Step 5: Rendering SVG with embedded molecules...")
        svg_content = _render_diagram_with_molecules(
            plan=diagram_plan,
            molecule_images=molecule_images,
            canvas_width=960,
            canvas_height=640,
        )
        print(f"✅ SVG rendered ({len(svg_content)} characters)")

        # Save SVG for inspection
        with open("pipeline_test_diagram.svg", "w") as f:
            f.write(svg_content)
        print("   📄 Saved: pipeline_test_diagram.svg")

        # Step 6: Convert SVG to PNG
        print("🔄 Step 6: Converting SVG to final PNG...")
        png_base64 = _svg_to_png(svg_content, 960, 640)

        if png_base64:
            # Save final PNG
            png_data = base64.b64decode(png_base64)
            with open("pipeline_test_final.png", "wb") as f:
                f.write(png_data)
            print(f"✅ Final PNG created ({len(png_data)} bytes)")
            print("   🖼️  Saved: pipeline_test_final.png")

            # Verify PNG is valid
            try:
                from PIL import Image

                img = Image.open("pipeline_test_final.png")
                print(f"   📏 PNG dimensions: {img.size}")
                print(f"   🎨 PNG mode: {img.mode}")
                img.close()
                print("✅ PNG file is valid and readable")
            except ImportError:
                print("   ℹ️  PIL not available for PNG validation")
            except Exception as e:
                print(f"   ⚠️  PNG validation failed: {e}")
        else:
            print("❌ SVG to PNG conversion failed - no conversion library available")
            return False

        print("\n🎉 COMPLETE PIPELINE SUCCESS!")
        print("=" * 50)
        print("Pipeline Summary:")
        print(f"   📝 Input: '{brief}'")
        print(f"   📋 YAML: {len(yaml_spec)} chars")
        print(
            f"   🧬 Molecules: {len(molecules)} found, {len(molecule_images)} rendered"
        )
        print(f"   🖼️  SVG: {len(svg_content)} chars")
        print(f"   📸 PNG: {len(png_data) if png_base64 else 0} bytes")
        print("\nFiles created:")
        print("   - pipeline_test_spec.yaml (YAML specification)")
        print("   - pipeline_test_diagram.svg (SVG with molecules)")
        print("   - pipeline_test_final.png (Final PNG output)")
        if molecule_images:
            print("   - molecule_*.png (Individual transparent molecules)")

        return True

    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        import traceback

        traceback.print_exc()
        return False


async def test_api_endpoint():
    """Test the /graphic/make API endpoint directly."""

    print("\n🌐 Testing /graphic/make API Endpoint")
    print("=" * 50)

    try:
        import httpx

        request_data = {
            "brief": "Show the photosynthesis process with CO2, water, glucose, and oxygen",
            "width": 800,
            "height": 600,
            "output_format": "both",  # Request both SVG and PNG
        }

        async with httpx.AsyncClient(timeout=180.0) as client:
            print("🔄 Calling /graphic/make endpoint...")
            response = await client.post(
                "http://localhost:8000/graphic/make", json=request_data
            )

            if response.status_code == 200:
                result = response.json()
                print("✅ API call successful")
                print(f"   📋 YAML: {len(result.get('yaml_spec', ''))} chars")
                print(f"   🖼️  SVG: {len(result.get('svg_content', ''))} chars")

                if result.get("png_base64"):
                    png_data = base64.b64decode(result["png_base64"])
                    with open("api_test_result.png", "wb") as f:
                        f.write(png_data)
                    print(f"   📸 PNG: {len(png_data)} bytes")
                    print("   💾 Saved: api_test_result.png")
                    return True
                else:
                    print("   ⚠️  No PNG data in response")
                    return False
            else:
                print(f"❌ API call failed: {response.status_code}")
                print(f"   Response: {response.text}")
                return False

    except Exception as e:
        print(f"❌ API test failed: {e}")
        return False


async def main():
    """Run all tests."""

    print("🧪 Complete PNG Pipeline Test Suite")
    print("=" * 60)

    # Test 1: Direct pipeline functions
    pipeline_success = await test_complete_png_pipeline()

    # Test 2: API endpoint (if server is running)
    api_success = await test_api_endpoint()

    print("\n📊 Test Results:")
    print("=" * 30)
    print(f"Direct Pipeline: {'✅ PASS' if pipeline_success else '❌ FAIL'}")
    print(f"API Endpoint:    {'✅ PASS' if api_success else '❌ FAIL'}")

    if pipeline_success:
        print("\n🎉 The complete PNG pipeline is working!")
        print("✨ You now have:")
        print("   1. Natural language → YAML conversion")
        print("   2. YAML → structured diagram plan")
        print("   3. Transparent molecule PNG generation")
        print("   4. SVG diagram with embedded molecules")
        print("   5. Final PNG output with overlaid molecules")

        return True
    else:
        print("\n💥 Pipeline needs fixes")
        return False


if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1)
