#!/usr/bin/env python3
"""
Isolated test of LLM service with o3-mini model and structured outputs.
"""

import os
import sys

sys.path.append("api")

from typing import List

from pydantic import BaseModel

from api.agent_management.llm_service import (
    LLMModelConfig,
    LLMService,
    ProviderType,
    StructuredLLMRequest,
)


class PyMOLCommand(BaseModel):
    """Single PyMOL command with description."""

    command: str
    description: str


class PyMOLScript(BaseModel):
    """Complete PyMOL script for molecule visualization."""

    molecule_name: str
    commands: List[PyMOLCommand]
    expected_output: str


def test_llm_service():
    """Test LLM service with o3-mini model."""
    print("🧪 Testing LLM Service with o3-mini model...")

    # Create LLM service with o3-mini configuration
    config = LLMModelConfig(provider=ProviderType.OPENAI, model_name="o3-mini")

    service = LLMService(config)
    print(f"✅ LLM Service created with model: {config.model_name}")

    # Test structured output generation
    test_molecules = ["glucose", "caffeine", "water", "CO2"]

    for molecule in test_molecules:
        print(f"\n🔬 Testing molecule: {molecule}")

        prompt = f"""Generate PyMOL commands to create and visualize the molecule '{molecule}'.

The commands should:
1. Create or load the molecule structure
2. Set appropriate representation (sticks, spheres, etc.)
3. Set colors and styling
4. Orient the view properly
5. Save as PNG with transparent background

Return a complete PyMOL script with individual commands and descriptions."""

        try:
            request = StructuredLLMRequest(
                user_prompt=prompt,
                system_prompt="You are a PyMOL expert. Generate precise PyMOL commands for molecular visualization.",
                response_model=PyMOLScript,
                llm_config=config,
                temperature=0.1,  # Low temperature for consistency
            )

            result = service.generate_structured(request)

            print(f"  ✅ Generated script for {result.molecule_name}")
            print(f"  📝 Commands ({len(result.commands)}):")
            for i, cmd in enumerate(result.commands[:3]):  # Show first 3 commands
                print(f"    {i+1}. {cmd.command}")
                print(f"       → {cmd.description}")
            if len(result.commands) > 3:
                print(f"    ... and {len(result.commands) - 3} more commands")
            print(f"  🎯 Expected: {result.expected_output}")

        except Exception as e:
            print(f"  ❌ Error for {molecule}: {str(e)}")
            return False

    print("\n🎉 All LLM tests passed!")
    return True


if __name__ == "__main__":
    success = test_llm_service()
    sys.exit(0 if success else 1)
