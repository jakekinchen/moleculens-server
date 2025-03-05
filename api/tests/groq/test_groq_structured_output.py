#!/usr/bin/env python3
"""
Test script to verify structured output generation with the Qwen 2.5 32B model.
"""

import os
import sys
from typing import List, Dict, Any, Optional
from pydantic import BaseModel, Field

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from agent_management.model_config import ModelRegistry, ModelCategory
from agent_management.llm_service import ProviderType, LLMService, LLMRequest, LLMModelConfig, StructuredLLMRequest

# Define some Pydantic models for structured output testing
class UserInfo(BaseModel):
    """Simple user information model."""
    name: str = Field(description="User's full name")
    age: int = Field(description="User's age in years")
    email: str = Field(description="User's email address")
    is_active: bool = Field(description="Whether the user is active")

class Ingredient(BaseModel):
    """Ingredient model for recipes."""
    name: str = Field(description="Name of the ingredient")
    quantity: float = Field(description="Quantity of the ingredient")
    unit: str = Field(description="Unit of measurement")

class Recipe(BaseModel):
    """Recipe model with nested ingredients."""
    title: str = Field(description="Title of the recipe")
    description: str = Field(description="Description of the recipe")
    prep_time: int = Field(description="Preparation time in minutes")
    cook_time: int = Field(description="Cooking time in minutes")
    servings: int = Field(description="Number of servings")
    ingredients: List[Ingredient] = Field(description="List of ingredients")
    instructions: List[str] = Field(description="List of instructions")

def test_structured_output():
    """Test structured output generation with the Qwen 2.5 32B model."""
    print(f"\n=== Testing Qwen 2.5 32B structured output generation ===")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("❌ GROQ_API_KEY environment variable is not set. Skipping test.")
        return False
    
    try:
        # Get model info
        model_name = "qwen-2.5-32b"
        model_info = ModelRegistry.create_instance(model_name)
        print(f"Model display name: {model_info.display_name}")
        print(f"Context length: {model_info.context_length}")
        categories = [cat.value for cat in model_info.categories]
        print(f"Categories: {', '.join(categories)}")
        
        # Create LLM config
        llm_config = LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name=model_info.name,
            api_key=os.getenv("GROQ_API_KEY")
        )
        
        # Create LLM service
        llm_service = LLMService(config=llm_config)
        
        # Test cases for structured output
        test_cases = [
            {
                "name": "Simple UserInfo model",
                "prompt": "Generate information for a fictional user named John Smith who is 35 years old, has email john.smith@example.com, and is an active user.",
                "model": UserInfo
            },
            {
                "name": "Complex Recipe model",
                "prompt": "Generate a recipe for chocolate chip cookies.",
                "model": Recipe
            }
        ]
        
        # Run test cases
        for i, test_case in enumerate(test_cases):
            print(f"\n--- Test Case {i+1}: {test_case['name']} ---")
            
            print(f"Prompt: {test_case['prompt']}")
            print(f"Model: {test_case['model'].__name__}")
            print("\nSending request to Groq API...")
            
            # Create a structured request
            request = StructuredLLMRequest(
                user_prompt=test_case['prompt'],
                system_prompt=f"Generate a valid {test_case['model'].__name__} object based on the prompt.",
                response_model=test_case['model'],
                llm_config=llm_config,
                max_tokens=1000
            )
            
            # Generate structured response
            response = llm_service.generate_structured(request)
            
            print("\nResponse:")
            print(response.model_dump_json(indent=2))
        
        print("\n✅ All structured output tests completed!")
        return True
    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        return False

if __name__ == "__main__":
    test_structured_output() 