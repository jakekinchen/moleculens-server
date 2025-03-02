"""
Test script for structured output with the Groq provider
"""

import os
import sys
from typing import List, Optional
from pydantic import BaseModel
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.llm_service import LLMService, StructuredLLMRequest, LLMModelConfig, ProviderType

# Load environment variables
load_dotenv()

# Define a structured output model
class Recipe(BaseModel):
    name: str
    ingredients: List[str]
    steps: List[str]
    prep_time_minutes: Optional[int] = None
    cook_time_minutes: Optional[int] = None

def test_groq_structured():
    """Test structured output with the Groq provider"""
    print("Testing Structured Output with Groq Provider...")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("Error: GROQ_API_KEY environment variable is not set")
        return False
    
    # Check if groq package is installed
    try:
        import groq
        print("Groq package is installed")
    except ImportError:
        print("Error: groq package is not installed")
        print("Install with: pip install groq")
        return False
    
    try:
        # Create an LLM service with Groq provider
        llm_service = LLMService(
            config=LLMModelConfig(
                provider=ProviderType.GROQ,
                model_name="llama3-70b-8192"
            )
        )
        
        # Create a structured request
        request = StructuredLLMRequest(
            user_prompt="Give me a recipe for a simple chocolate chip cookie.",
            system_prompt="You are a helpful cooking assistant that provides recipes in a structured format.",
            response_model=Recipe,
            llm_config=LLMModelConfig(
                provider=ProviderType.GROQ,
                model_name="llama3-70b-8192"
            )
        )
        
        # Generate a structured response
        recipe = llm_service.generate_structured(request)
        
        # Print the response
        print("\nGroq LLama-3 Structured Response:")
        print("-" * 50)
        print(f"Recipe: {recipe.name}")
        print("\nIngredients:")
        for ingredient in recipe.ingredients:
            print(f"- {ingredient}")
        print("\nSteps:")
        for i, step in enumerate(recipe.steps, 1):
            print(f"{i}. {step}")
        print(f"\nPrep Time: {recipe.prep_time_minutes} minutes")
        print(f"Cook Time: {recipe.cook_time_minutes} minutes")
        print("-" * 50)
        
        return True
    except Exception as e:
        print(f"Error testing structured output with Groq provider: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_groq_structured()
    if success:
        print("\n✅ Structured output with Groq provider test passed")
    else:
        print("\n❌ Structured output with Groq provider test failed") 