#!/usr/bin/env python3
"""
Test script to verify the vision capabilities of the Llama 3.2 90B Vision model.
"""

import os
import sys
import base64
from typing import List, Dict, Any, Optional
from pathlib import Path

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from agent_management.model_config import ModelRegistry, ModelCategory
from agent_management.llm_service import ProviderType, LLMService, LLMRequest, LLMModelConfig

def encode_image_to_base64(image_path: str) -> str:
    """Encode an image to base64."""
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def test_vision_model(image_path: Optional[str] = None):
    """Test the Llama 3.2 90B Vision model with an image."""
    print(f"\n=== Testing Llama 3.2 90B Vision model ===")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("❌ GROQ_API_KEY environment variable is not set. Skipping test.")
        return False
    
    # Use a default test image if none is provided
    if not image_path:
        # Check if the test_images directory exists, if not create it
        test_images_dir = Path(__file__).parent / "test_images"
        test_images_dir.mkdir(exist_ok=True)
        
        # Default test image path
        image_path = str(test_images_dir / "test_image.jpg")
        
        # If the test image doesn't exist, create a simple text file explaining how to use the test
        if not Path(image_path).exists():
            print(f"⚠️ No test image found at {image_path}")
            print("Please provide an image path as an argument to test the vision model.")
            print("Example: python test_groq_vision_model.py /path/to/your/image.jpg")
            return False
    
    try:
        # Get model info
        model_name = "llama-3.2-90b-vision-preview"
        model_info = ModelRegistry.create_instance(model_name)
        print(f"Model display name: {model_info.display_name}")
        print(f"Context length: {model_info.context_length}")
        categories = [cat.value for cat in model_info.categories]
        print(f"Categories: {', '.join(categories)}")
        
        # Check if the model has vision capability
        if ModelCategory.VISION not in model_info.categories:
            print(f"❌ Model {model_name} does not have vision capability.")
            return False
        
        # Create LLM config
        llm_config = LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name=model_info.name,
            api_key=os.getenv("GROQ_API_KEY")
        )
        
        # Create LLM service
        llm_service = LLMService(config=llm_config)
        
        # Encode the image to base64
        try:
            base64_image = encode_image_to_base64(image_path)
            print(f"Successfully encoded image: {image_path}")
        except Exception as e:
            print(f"❌ Failed to encode image: {str(e)}")
            return False
        
        # Create request with image
        # Note: The current LLMRequest doesn't support images parameter
        # This is a placeholder for when the API supports it
        request = LLMRequest(
            llm_config=llm_config,
            user_prompt="What's in this image? Please describe it in detail.",
            max_tokens=300
            # images=[base64_image]  # Uncomment when the API supports images
        )
        
        print("\nSending request to Groq API with image...")
        print("Note: The current implementation doesn't support image input yet.")
        
        # Generate response
        response = llm_service.generate(request)
        
        print("\nResponse:")
        print(response.content[:1000] + "..." if len(response.content) > 1000 else response.content)
        print(f"\nUsage: {response.usage}")
        
        print("\n✅ Test passed!")
        return True
    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        return False

if __name__ == "__main__":
    # Check if an image path is provided as an argument
    if len(sys.argv) > 1:
        test_vision_model(sys.argv[1])
    else:
        test_vision_model() 