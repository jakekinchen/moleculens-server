#!/usr/bin/env python3
"""
Test script to verify the code generation capabilities of the DeepSeek R1 Distill Llama 70B model.
"""

import os
import sys
from typing import List, Dict, Any

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from agent_management.model_config import ModelRegistry, ModelCategory
from agent_management.llm_service import ProviderType, LLMService, LLMRequest, LLMModelConfig

def test_code_generation():
    """Test the code generation capabilities of the DeepSeek R1 Distill Llama 70B model."""
    print(f"\n=== Testing DeepSeek R1 Distill Llama 70B code generation ===")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("❌ GROQ_API_KEY environment variable is not set. Skipping test.")
        return False
    
    try:
        # Get model info
        model_name = "deepseek-r1-distill-llama-70b"
        model_info = ModelRegistry.create_instance(model_name)
        print(f"Model display name: {model_info.display_name}")
        print(f"Context length: {model_info.context_length}")
        categories = [cat.value for cat in model_info.categories]
        print(f"Categories: {', '.join(categories)}")
        
        # Check if the model has code capability
        if ModelCategory.CODE not in model_info.categories:
            print(f"❌ Model {model_name} does not have code generation capability.")
            return False
        
        # Create LLM config
        llm_config = LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name=model_info.name,
            api_key=os.getenv("GROQ_API_KEY")
        )
        
        # Create LLM service
        llm_service = LLMService(config=llm_config)
        
        # Test cases for code generation
        test_cases = [
            {
                "name": "Simple Python function",
                "prompt": "Write a Python function to calculate the Fibonacci sequence up to n terms."
            },
            {
                "name": "Data processing",
                "prompt": "Write a Python function that takes a list of dictionaries with 'name' and 'age' keys, and returns a new list sorted by age in descending order."
            },
            {
                "name": "API endpoint",
                "prompt": "Create a FastAPI endpoint that accepts a POST request with a JSON body containing 'text' field and returns sentiment analysis (positive, negative, or neutral)."
            }
        ]
        
        # Run test cases
        for i, test_case in enumerate(test_cases):
            print(f"\n--- Test Case {i+1}: {test_case['name']} ---")
            
            # Create request
            request = LLMRequest(
                llm_config=llm_config,
                user_prompt=test_case['prompt'],
                max_tokens=500
            )
            
            print(f"Prompt: {test_case['prompt']}")
            print("\nSending request to Groq API...")
            
            # Generate response
            response = llm_service.generate(request)
            
            print("\nResponse:")
            print(response.content[:1000] + "..." if len(response.content) > 1000 else response.content)
            print(f"\nUsage: {response.usage}")
        
        print("\n✅ All code generation tests completed!")
        return True
    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        return False

if __name__ == "__main__":
    test_code_generation() 