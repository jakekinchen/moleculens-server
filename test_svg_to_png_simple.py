#!/usr/bin/env python3
"""
Simple test of SVG to PNG conversion functionality.
"""

import base64
import os
import sys

# Add the api directory to the path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))

from api.routers.graphic.routes import _svg_to_png


def test_svg_to_png():
    """Test SVG to PNG conversion with a simple SVG."""

    # Create a simple test SVG
    test_svg = """<?xml version="1.0" encoding="UTF-8"?>
<svg width="400" height="300" xmlns="http://www.w3.org/2000/svg">
  <rect width="400" height="300" fill="white"/>
  <circle cx="100" cy="100" r="50" fill="red" stroke="black" stroke-width="2"/>
  <circle cx="300" cy="100" r="50" fill="blue" stroke="black" stroke-width="2"/>
  <line x1="150" y1="100" x2="250" y2="100" stroke="black" stroke-width="3" marker-end="url(#arrowhead)"/>
  <text x="200" y="200" text-anchor="middle" font-family="Arial" font-size="16">Test Diagram</text>
  <defs>
    <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto">
      <polygon points="0 0, 10 3.5, 0 7" fill="black"/>
    </marker>
  </defs>
</svg>"""

    print("🧪 Testing SVG to PNG Conversion")
    print("=" * 40)

    try:
        print("🔄 Converting SVG to PNG...")
        png_base64 = _svg_to_png(test_svg, 400, 300)

        if png_base64:
            # Save the PNG
            png_data = base64.b64decode(png_base64)
            with open("svg_to_png_test.png", "wb") as f:
                f.write(png_data)

            print(f"✅ PNG conversion successful ({len(png_data)} bytes)")
            print("   💾 Saved: svg_to_png_test.png")

            # Try to validate with PIL if available
            try:
                from PIL import Image

                img = Image.open("svg_to_png_test.png")
                print(f"   📏 PNG dimensions: {img.size}")
                print(f"   🎨 PNG mode: {img.mode}")
                img.close()
                print("✅ PNG file is valid")
                return True
            except ImportError:
                print("   ℹ️  PIL not available for validation")
                return True
            except Exception as e:
                print(f"   ⚠️  PNG validation failed: {e}")
                return False
        else:
            print("❌ SVG to PNG conversion failed")
            print("   💡 Try installing: pip install cairosvg")
            return False

    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_svg_to_png()
    if success:
        print("\n🎉 SVG to PNG conversion is working!")
    else:
        print("\n💥 SVG to PNG conversion needs fixes")
    sys.exit(0 if success else 1)
