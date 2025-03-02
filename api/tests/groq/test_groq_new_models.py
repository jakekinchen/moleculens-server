#!/usr/bin/env python3
"""
Test script to verify that the new models can be used with the Groq provider.
"""

import os
import sys
from typing import List, Dict, Any

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from agent_management.model_config import ModelRegistry, ModelCategory
from agent_management.llm_service import ProviderType, LLMService, LLMRequest, LLMModelConfig

def test_groq_model(model_name: str):
    """Test a specific Groq model."""
    print(f"\n=== Testing model: {model_name} ===")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("❌ GROQ_API_KEY environment variable is not set. Skipping test.")
        return False
    
    try:
        # Get model info
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
        
        # Create request
        request = LLMRequest(
            llm_config=llm_config,
            user_prompt="Hello, can you introduce yourself and tell me what capabilities you have?",
            max_tokens=100
        )
        
        print("\nSending request to Groq API...")
        
        # Generate response
        response = llm_service.generate(request)
        
        print("\nResponse:")
        print(response.content[:500] + "..." if len(response.content) > 500 else response.content)
        print(f"\nUsage: {response.usage}")
        
        print("\n✅ Test passed!")
        return True
    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        return False

def test_all_new_models():
    """Test all new Groq models."""
    # Models to test
    models_to_test = [
        "deepseek-r1-distill-llama-70b",
        "llama-3.2-90b-vision-preview",
        "qwen-2.5-32b",
        "llama-3.3-70b-specdec"
    ]
    
    results = {}
    
    for model_name in models_to_test:
        results[model_name] = test_groq_model(model_name)
    
    # Print summary
    print("\n=== Test Summary ===")
    for model_name, result in results.items():
        print(f"{model_name}: {'✅ Passed' if result else '❌ Failed'}")
    
    # Return True if all tests passed, False otherwise
    return all(results.values())

if __name__ == "__main__":
    # Check if a specific model is specified
    if len(sys.argv) > 1:
        test_groq_model(sys.argv[1])
    else:
        test_all_new_models() 