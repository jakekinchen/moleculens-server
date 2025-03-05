"""
Test script for using the LLMService with the Groq provider
"""

import os
import sys
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.llm_service import LLMService, LLMRequest, LLMModelConfig, ProviderType

# Load environment variables
load_dotenv()

def test_groq_service():
    """Test the LLMService with the Groq provider"""
    print("Testing LLMService with Groq Provider...")
    
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
        
        # Create a simple request
        request = LLMRequest(
            user_prompt="What are the key differences between supervised and unsupervised learning in machine learning?",
            system_prompt="You are a helpful AI assistant specializing in machine learning concepts."
        )
        
        # Generate a response
        response = llm_service.generate(request)
        
        # Print the response
        print("\nGroq LLama-3 Response via LLMService:")
        print("-" * 50)
        print(response.content)
        print("-" * 50)
        print(f"Model: {response.model}")
        print(f"Tokens: {response.usage}")
        
        return True
    except Exception as e:
        print(f"Error testing LLMService with Groq provider: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_groq_service()
    if success:
        print("\n✅ LLMService with Groq provider test passed")
    else:
        print("\n❌ LLMService with Groq provider test failed") 