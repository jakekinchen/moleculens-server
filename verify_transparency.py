#!/usr/bin/env python3
"""Verify PNG transparency by checking alpha channel values."""

import sys
from pathlib import Path

try:
    import numpy as np
    from PIL import Image
except ImportError:
    print("Installing required packages...")
    import subprocess

    subprocess.check_call([sys.executable, "-m", "pip", "install", "Pillow", "numpy"])
    import numpy as np
    from PIL import Image


def check_transparency(png_path):
    """Check if a PNG file has transparent pixels."""
    try:
        img = Image.open(png_path)

        # Convert to RGBA if not already
        if img.mode != "RGBA":
            img = img.convert("RGBA")

        # Get alpha channel
        alpha_channel = np.array(img)[:, :, 3]

        # Check for transparency
        min_alpha = alpha_channel.min()
        max_alpha = alpha_channel.max()
        transparent_pixels = np.sum(alpha_channel < 255)
        total_pixels = alpha_channel.size

        print(f"\n=== {png_path} ===")
        print(f"Image size: {img.size}")
        print(f"Mode: {img.mode}")
        print(f"Alpha range: {min_alpha} - {max_alpha}")
        print(f"Transparent pixels: {transparent_pixels:,} / {total_pixels:,}")
        print(f"Transparency percentage: {(transparent_pixels/total_pixels)*100:.2f}%")

        if transparent_pixels > 0:
            print("✅ HAS TRANSPARENCY")
        else:
            print("❌ NO TRANSPARENCY")

        return transparent_pixels > 0

    except Exception as e:
        print(f"❌ Error checking {png_path}: {e}")
        return False


def main():
    """Check transparency of generated PNG files."""
    print("🔍 Verifying PNG Transparency")

    png_files = ["overview_test.png", "overview_opaque.png", "transparent_molecule.png"]

    results = {}
    for png_file in png_files:
        if Path(png_file).exists():
            results[png_file] = check_transparency(png_file)
        else:
            print(f"\n❌ {png_file} not found")
            results[png_file] = None

    print("\n" + "=" * 50)
    print("📊 SUMMARY")
    print("=" * 50)

    for png_file, has_transparency in results.items():
        if has_transparency is None:
            status = "NOT FOUND"
        elif has_transparency:
            status = "✅ TRANSPARENT"
        else:
            status = "❌ OPAQUE"
        print(f"{png_file:25} {status}")

    # Check if we have working transparency
    transparent_files = [f for f, t in results.items() if t is True]
    if transparent_files:
        print(f"\n🎉 SUCCESS: {len(transparent_files)} file(s) have transparency!")
        print("The PyMOL transparent background rendering is working correctly.")
    else:
        print(
            "\n⚠️  No transparent files found. Check PyMOL ray_opaque_background setting."
        )


if __name__ == "__main__":
    main()
