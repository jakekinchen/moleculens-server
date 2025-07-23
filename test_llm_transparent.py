#!/usr/bin/env python3
"""Test LLM understanding of transparent background requests."""

import os
import sys

sys.path.append("api")

from agent_management.pymol_translator import _spec_from_prompt


def test_llm_transparent_understanding():
    """Test if the LLM can understand and generate transparent background requests."""

    test_prompts = [
        "Show me caffeine with a transparent background for a presentation",
        "Create a molecule visualization with no background so I can overlay it",
        "Generate a transparent PNG of 1ubq protein structure",
        "Make a publication quality image with transparent background",
        "Show binding site of 1abc with transparent background for infographic use",
    ]

    print("Testing LLM understanding of transparent background requests...\n")

    for i, prompt in enumerate(test_prompts, 1):
        print(f"=== Test {i} ===")
        print(f"Prompt: {prompt}")

        try:
            spec = _spec_from_prompt(prompt)
            print(f"Operation: {spec.op}")
            print(f"Structure ID: {spec.structure_id}")
            print(f"Transparent Background: {spec.rendering.transparent_background}")
            print(f"Ray Trace: {spec.rendering.ray_trace}")
            print(f"Resolution: {spec.rendering.resolution}")
            print(f"DPI: {spec.rendering.dpi}")
            print(f"Background Color: {spec.rendering.background_color}")

            if spec.rendering.transparent_background:
                print("✅ LLM correctly identified transparent background request")
            else:
                print("⚠️  LLM did not identify transparent background request")

        except Exception as e:
            print(f"❌ Error: {e}")

        print()


if __name__ == "__main__":
    test_llm_transparent_understanding()
