#!/usr/bin/env python3
"""
Test script to verify that the DeepSeek extraction works with the caption agent.
"""

import os
import sys
import json

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.model_config import get_llm_service
from api.agent_management.agents.caption_agent import CaptionAgent

def test_caption_agent_with_deepseek():
    """Test the caption agent with a DeepSeek model."""
    print("\n=== Testing Caption Agent with DeepSeek Model ===")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("GROQ_API_KEY environment variable not set, skipping test")
        return
    
    try:
        # Get the DeepSeek model service
        llm_service = get_llm_service("deepseek-r1-distill-llama-70b")
        
        # Create a caption agent
        caption_agent = CaptionAgent(llm_service)
        
        # Test with a simple prompt
        prompt = "A rotating planet Earth showing continents and oceans"
        
        print(f"Generating captions for prompt: '{prompt}'")
        
        # Get the caption snippet
        caption_snippet = caption_agent.get_caption_snippet(prompt)
        
        # Print the caption snippet
        print("\nGenerated Caption Snippet:")
        print("=" * 80)
        print(caption_snippet)
        print("=" * 80)
        
        # Check if the extraction worked
        if "<think" in caption_snippet or "<thinking" in caption_snippet:
            print("WARNING: Thinking tags found in the output. Extraction may have failed.")
        else:
            print("SUCCESS: No thinking tags found in the output. Extraction worked correctly.")
        
        # Try to extract the captions array from the JavaScript code
        try:
            # Find the captions array in the JavaScript code
            captions_start = caption_snippet.find("const captions = ") + len("const captions = ")
            captions_end = caption_snippet.find(";", captions_start)
            captions_json = caption_snippet[captions_start:captions_end].strip()
            
            # Parse the JSON
            captions = json.loads(captions_json)
            
            # Print the captions
            print("\nExtracted Captions:")
            for caption in captions:
                print(f"Time: {caption['time']}s, Text: '{caption['text']}'")
                
            print(f"\nTotal captions: {len(captions)}")
        except Exception as e:
            print(f"Error extracting captions from JavaScript: {str(e)}")
        
    except Exception as e:
        print(f"Error testing caption agent with DeepSeek model: {str(e)}")

if __name__ == "__main__":
    test_caption_agent_with_deepseek() 