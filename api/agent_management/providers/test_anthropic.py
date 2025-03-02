"""
Test script for the Anthropic provider using the geometry agent.
"""

import os
from dotenv import load_dotenv
from agent_management.agent_factory import AgentFactory
from agent_management.llm_service import LLMModelConfig, ProviderType

def test_geometry_agent():
    """Test the geometry agent with the Anthropic provider"""
    
    # Load environment variables
    load_dotenv()
    
    # Create a geometry agent with Anthropic's Claude model
    agent = AgentFactory.create_geometry_agent("claude-3-7-sonnet-latest")
    
    # Test with a simple molecule
    prompt = "Create a water molecule (H2O) with appropriate bond angles"
    
    try:
        # Generate geometry code
        geometry_code = agent.get_geometry_snippet(prompt)
        
        # Print the result
        print("\nGenerated Geometry Code:")
        print("------------------------")
        print(geometry_code)
        print("------------------------")
        
        return True
    except Exception as e:
        print(f"Error testing geometry agent: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_geometry_agent()
    print(f"\nTest {'succeeded' if success else 'failed'}") 