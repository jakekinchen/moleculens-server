#!/usr/bin/env python3
"""Test script for transparent background rendering."""

import json
import sys
from pathlib import Path

import requests


def test_transparent_rendering():
    """Test the transparent background rendering functionality."""

    # Test data
    test_cases = [
        {
            "description": "Show me a transparent molecule of caffeine with no background",
            "format": "image",
            "transparent_background": True,
            "ray_trace": True,
            "resolution": [800, 600],
            "dpi": 150,
            "background_color": "white",
        },
        {
            "description": "Create a publication quality image of 1ubq with transparent background",
            "format": "image",
            "transparent_background": True,
            "ray_trace": True,
            "ray_trace_mode": "default",
            "resolution": [1920, 1080],
            "dpi": 300,
            "antialias": True,
            "ray_shadow": True,
            "depth_cue": True,
        },
    ]

    base_url = "http://localhost:8000"

    for i, test_case in enumerate(test_cases):
        print(f"\n=== Test Case {i+1} ===")
        print(f"Description: {test_case['description']}")
        print(f"Transparent: {test_case['transparent_background']}")

        try:
            response = requests.post(f"{base_url}/render", json=test_case, timeout=60)

            if response.status_code == 200:
                # Save the response
                output_file = f"transparent_test_{i+1}.png"
                with open(output_file, "wb") as f:
                    f.write(response.content)
                print(f"✅ Success! Saved to {output_file}")

                # Check file size (transparent PNGs should have alpha channel)
                file_size = Path(output_file).stat().st_size
                print(f"File size: {file_size} bytes")

            else:
                print(f"❌ Error: {response.status_code}")
                print(f"Response: {response.text}")

        except requests.exceptions.RequestException as e:
            print(f"❌ Request failed: {e}")
        except Exception as e:
            print(f"❌ Unexpected error: {e}")


if __name__ == "__main__":
    print("Testing transparent background rendering...")
    test_transparent_rendering()
    print("\nTest completed!")
