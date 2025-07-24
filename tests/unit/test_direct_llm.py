#!/usr/bin/env python3
"""
Test direct LLM utility without importing the complex agent system.
"""

import os
import sys
from typing import Optional, Type, TypeVar

from pydantic import BaseModel

# Add the api directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))

from agent_management.llm_service import (
    LLMModelConfig,
    ProviderType,
    StructuredLLMRequest,
)

# Import only what we need directly
from agent_management.providers.openai_provider import OpenAIProvider

T = TypeVar("T", bound=BaseModel)


def generate_structured_direct(
    user_prompt: str,
    response_model: Type[T],
    system_prompt: Optional[str] = None,
    model_name: str = "o3-mini",
    temperature: Optional[float] = None,
    max_tokens: Optional[int] = None,
    api_key: Optional[str] = None,
) -> T:
    """Generate structured output using OpenAI provider directly."""
    # Create provider
    provider = OpenAIProvider(api_key=api_key)

    # Create config
    config = LLMModelConfig(
        provider=ProviderType.OPENAI, model_name=model_name, api_key=api_key
    )

    # Create request
    request = StructuredLLMRequest(
        user_prompt=user_prompt,
        system_prompt=system_prompt,
        response_model=response_model,
        llm_config=config,
        temperature=temperature,
        max_tokens=max_tokens,
    )

    # Generate structured response
    return provider.generate_structured(request)


# Test models
class MoleculeInfo(BaseModel):
    name: str
    formula: str
    description: str


class PyMOLCommand(BaseModel):
    molecule: str
    commands: str


def test_molecule_info():
    """Test generating molecule information."""
    print("🧪 Testing molecule information generation...")

    try:
        result = generate_structured_direct(
            user_prompt="Describe the glucose molecule",
            response_model=MoleculeInfo,
            system_prompt="You are a chemistry expert. Provide accurate molecular information.",
            model_name="o3-mini",
        )

        print(f"✅ Molecule: {result.name}")
        print(f"✅ Formula: {result.formula}")
        print(f"✅ Description: {result.description[:100]}...")
        return True

    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False


def test_pymol_commands():
    """Test generating PyMOL commands for different molecules."""
    print("\n🧪 Testing PyMOL command generation...")

    molecules = ["glucose", "caffeine", "water", "CO2"]

    for molecule in molecules:
        try:
            result = generate_structured_direct(
                user_prompt=f"Generate PyMOL commands to create and visualize {molecule} molecule",
                response_model=PyMOLCommand,
                system_prompt="You are a PyMOL expert. Generate precise commands for molecular visualization.",
                model_name="o3-mini",
            )

            print(f"✅ {molecule:8} -> {result.commands[:50]}...")

        except Exception as e:
            print(f"❌ {molecule:8} -> Error: {str(e)}")
            return False

    return True


def test_yaml_generation():
    """Test generating YAML for diagram planning."""
    print("\n🧪 Testing YAML diagram generation...")

    class YAMLResponse(BaseModel):
        content: str

    try:
        result = generate_structured_direct(
            user_prompt="Create a molecular diagram showing glucose and CO2 interaction",
            response_model=YAMLResponse,
            system_prompt="""Generate a YAML specification for a molecular diagram.
Include meta information, canvas size, and diagram elements with molecules and arrows.""",
            model_name="o3-mini",
            max_tokens=1000,
        )

        print(f"✅ Generated YAML ({len(result.content)} chars):")
        print(
            result.content[:200] + "..."
            if len(result.content) > 200
            else result.content
        )
        return True

    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False


if __name__ == "__main__":
    print("🚀 Testing Direct LLM Utility with o3-mini")
    print("=" * 50)

    success = True
    success &= test_molecule_info()
    success &= test_pymol_commands()
    success &= test_yaml_generation()

    print("\n" + "=" * 50)
    if success:
        print("🎉 All tests passed! Direct LLM utility is working correctly.")
    else:
        print("❌ Some tests failed.")

    sys.exit(0 if success else 1)
