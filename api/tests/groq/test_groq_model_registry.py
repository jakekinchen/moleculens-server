#!/usr/bin/env python3
"""
Test script to verify the models we want to add to the Groq provider.
"""

import os
import sys
from typing import List, Dict, Any

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from agent_management.model_config import ModelRegistry, ModelCategory
from agent_management.llm_service import ProviderType

def test_groq_model_registry():
    """Test the Groq model registry to ensure all models are registered correctly."""
    print("Testing Groq model registry...")
    
    # Get all registered models
    all_models = ModelRegistry.list_models()
    
    # Filter for Groq models
    groq_models = []
    for model_name in all_models:
        model_info = ModelRegistry.create_instance(model_name)
        if model_info.provider == ProviderType.GROQ:
            groq_models.append(model_info)
    
    # Print all Groq models
    print(f"\nFound {len(groq_models)} Groq models:")
    for model in groq_models:
        categories = [cat.value for cat in model.categories]
        print(f"- {model.name} (Display: {model.display_name})")
        print(f"  Context length: {model.context_length}")
        print(f"  Categories: {', '.join(categories)}")
        print(f"  Default: {model.is_default}")
        print(f"  API params: {model.api_params}")
        print()
    
    # Check if the models we want to add are already registered
    models_to_check = [
        "deepseek-r1-distill-llama-70b",
        "llama-3.2-90b-vision-preview",
        "qwen-2.5-32b",
        "llama-3.3-70b-specdec"
    ]
    
    print("\nChecking for specific models:")
    all_found = True
    for model_name in models_to_check:
        found = False
        for model in groq_models:
            if model.name.lower() == model_name.lower():
                found = True
                print(f"✅ {model_name} is already registered")
                break
        
        if not found:
            all_found = False
            print(f"❌ {model_name} is not registered")
    
    print("\nTest completed.")
    return all_found

if __name__ == "__main__":
    test_groq_model_registry() 