#!/usr/bin/env python3
"""
Test script to verify that the DeepSeek extraction works with the geometry agent.
"""

import os
import sys

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.model_config import get_llm_service
from api.agent_management.agents.geometry_agent import GeometryAgent

def test_geometry_agent_with_deepseek():
    """Test the geometry agent with a DeepSeek model."""
    print("\n=== Testing Geometry Agent with DeepSeek Model ===")
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("GROQ_API_KEY environment variable not set, skipping test")
        return
    
    try:
        # Get the DeepSeek model service
        llm_service = get_llm_service("deepseek-r1-distill-llama-70b")
        
        # Create a geometry agent
        geometry_agent = GeometryAgent(llm_service)
        
        # Test with a simple prompt
        prompt = "Create a simple red sphere"
        
        print(f"Generating geometry for prompt: '{prompt}'")
        
        # Get the geometry snippet
        geometry_snippet = geometry_agent.get_geometry_snippet(prompt)
        
        # Print the geometry snippet
        print("\nGenerated Geometry Snippet:")
        print("=" * 80)
        print(geometry_snippet)
        print("=" * 80)
        
        # Check if the extraction worked
        if "<think" in geometry_snippet or "<thinking" in geometry_snippet:
            print("WARNING: Thinking tags found in the output. Extraction may have failed.")
        else:
            print("SUCCESS: No thinking tags found in the output. Extraction worked correctly.")
        
    except Exception as e:
        print(f"Error testing geometry agent with DeepSeek model: {str(e)}")

if __name__ == "__main__":
    test_geometry_agent_with_deepseek() 