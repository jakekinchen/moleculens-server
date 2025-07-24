#!/usr/bin/env python3
"""Test script for Calvin cycle diagram generation."""

import json
import os
import time

import requests


def test_calvin_cycle_diagram():
    """Test the Calvin cycle diagram generation through the new graphic endpoints."""

    # Base URL for the API
    base_url = "http://localhost:8000"

    # Test data
    calvin_cycle_brief = """
    Create an infographic about the Calvin cycle with several different 3D or 2D proteins/molecules.
    Include key molecules like:
    - RuBisCO (ribulose-1,5-bisphosphate carboxylase/oxygenase) enzyme
    - CO2 (carbon dioxide)
    - RuBP (ribulose-1,5-bisphosphate)
    - 3-PGA (3-phosphoglycerate)
    - G3P (glyceraldehyde-3-phosphate)
    - Glucose

    Show the three main phases:
    1. Carbon fixation
    2. Reduction
    3. Regeneration

    Use a clean scientific style with proper molecular representations.
    """

    print("🧪 Testing Calvin Cycle Diagram Generation")
    print("=" * 50)

    # Test 1: Generate YAML plan
    print("\n1. Testing /graphic/plan endpoint...")
    try:
        plan_response = requests.post(
            f"{base_url}/graphic/plan",
            json={
                "brief": calvin_cycle_brief,
                "width": 1200,
                "height": 800,
                "model_name": "gpt-3.5-turbo",
            },
            timeout=30,
        )

        if plan_response.status_code == 200:
            plan_data = plan_response.json()
            yaml_spec = plan_data.get("yaml_spec", "")
            print("✅ Plan generation successful!")
            print(f"📄 Generated YAML spec ({len(yaml_spec)} characters)")
            print("First 300 characters:")
            print(yaml_spec[:300] + "..." if len(yaml_spec) > 300 else yaml_spec)
        else:
            print(f"❌ Plan generation failed: {plan_response.status_code}")
            print(f"Error: {plan_response.text}")
            return False

    except Exception as e:
        print(f"❌ Plan generation error: {e}")
        return False

    # Test 2: Validate the generated YAML
    print("\n2. Testing /graphic/validate endpoint...")
    try:
        validate_response = requests.post(
            f"{base_url}/graphic/validate", json={"yaml_spec": yaml_spec}, timeout=10
        )

        if validate_response.status_code == 200:
            validate_data = validate_response.json()
            is_valid = validate_data.get("valid", False)
            errors = validate_data.get("errors", [])

            if is_valid:
                print("✅ YAML validation successful!")
            else:
                print(f"⚠️ YAML validation issues: {errors}")
        else:
            print(f"❌ Validation failed: {validate_response.status_code}")

    except Exception as e:
        print(f"❌ Validation error: {e}")

    # Test 3: Full pipeline with /graphic/make
    print("\n3. Testing /graphic/make endpoint (full pipeline)...")
    try:
        make_response = requests.post(
            f"{base_url}/graphic/make",
            json={
                "brief": calvin_cycle_brief,
                "width": 1200,
                "height": 800,
                "model_name": "gpt-3.5-turbo",
            },
            timeout=60,
        )

        if make_response.status_code == 200:
            make_data = make_response.json()
            svg_content = make_data.get("svg_content", "")
            yaml_spec_full = make_data.get("yaml_spec", "")

            print("✅ Full pipeline successful!")
            print(f"📊 Generated SVG ({len(svg_content)} characters)")
            print(f"📄 Generated YAML ({len(yaml_spec_full)} characters)")

            # Save outputs for inspection
            with open("calvin_cycle_output.svg", "w") as f:
                f.write(svg_content)
            with open("calvin_cycle_spec.yaml", "w") as f:
                f.write(yaml_spec_full)

            print(
                "💾 Saved outputs to calvin_cycle_output.svg and calvin_cycle_spec.yaml"
            )

            # Check for molecular content
            if any(
                molecule in svg_content.lower()
                for molecule in ["rubisco", "co2", "rubp", "glucose"]
            ):
                print("✅ SVG contains expected molecular content!")
            else:
                print("⚠️ SVG may not contain expected molecular content")

            return True

        else:
            print(f"❌ Full pipeline failed: {make_response.status_code}")
            print(f"Error: {make_response.text}")
            return False

    except Exception as e:
        print(f"❌ Full pipeline error: {e}")
        return False


def wait_for_server(base_url="http://localhost:8000", timeout=60):
    """Wait for the server to be ready."""
    print("⏳ Waiting for server to be ready...")
    start_time = time.time()

    while time.time() - start_time < timeout:
        try:
            response = requests.get(f"{base_url}/health", timeout=5)
            if response.status_code == 200:
                print("✅ Server is ready!")
                return True
        except:
            pass
        time.sleep(2)

    print("❌ Server not ready within timeout")
    return False


if __name__ == "__main__":
    # Check if server is running
    if wait_for_server():
        success = test_calvin_cycle_diagram()
        if success:
            print("\n🎉 All tests passed! Calvin cycle diagram generation is working.")
        else:
            print("\n❌ Some tests failed. Check the output above.")
    else:
        print("\n❌ Server is not running. Please start the Docker container first.")
        print(
            "Run: docker run -p 8000:8000 -e OPENAI_API_KEY=your_key moleculens-server"
        )
