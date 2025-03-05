#!/usr/bin/env python3
"""
Test script to verify the reasoning capabilities of the Llama 3.3 70B SpecDec model.
"""

import os
import sys
from typing import List, Dict, Any

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from agent_management.model_config import ModelRegistry, ModelCategory
from agent_management.llm_service import ProviderType, LLMService, LLMRequest, LLMModelConfig

def test_reasoning_capabilities():
    """Test the reasoning capabilities of the Llama 3.3 70B SpecDec model."""
    print(f"\n=== Testing Llama 3.3 70B SpecDec reasoning capabilities ===")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("❌ GROQ_API_KEY environment variable is not set. Skipping test.")
        return False
    
    try:
        # Get model info
        model_name = "llama-3.3-70b-specdec"
        model_info = ModelRegistry.create_instance(model_name)
        print(f"Model display name: {model_info.display_name}")
        print(f"Context length: {model_info.context_length}")
        categories = [cat.value for cat in model_info.categories]
        print(f"Categories: {', '.join(categories)}")
        
        # Check if the model has reasoning capability
        if ModelCategory.REASONING not in model_info.categories:
            print(f"❌ Model {model_name} does not have reasoning capability.")
            return False
        
        # Create LLM config
        llm_config = LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name=model_info.name,
            api_key=os.getenv("GROQ_API_KEY")
        )
        
        # Create LLM service
        llm_service = LLMService(config=llm_config)
        
        # Test cases for reasoning
        test_cases = [
            {
                "name": "Logical reasoning",
                "prompt": "If all A are B, and all B are C, what can we conclude about A and C?"
            },
            {
                "name": "Mathematical reasoning",
                "prompt": "A train travels at 60 mph. How far will it travel in 2.5 hours? Explain your reasoning step by step."
            },
            {
                "name": "Ethical dilemma",
                "prompt": "Consider the trolley problem: A trolley is headed toward five people who cannot move. You can pull a lever to redirect the trolley to a track where it will hit one person instead. What are the ethical considerations in this scenario? Analyze this from multiple perspectives."
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
        
        print("\n✅ All reasoning tests completed!")
        return True
    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        return False

if __name__ == "__main__":
    test_reasoning_capabilities() 