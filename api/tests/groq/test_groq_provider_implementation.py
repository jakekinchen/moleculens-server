"""
Test script to verify the Groq provider implementation
"""

import os
import sys
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.providers.groq_provider import GroqProvider
from api.agent_management.llm_service import LLMRequest, LLMModelConfig, ProviderType

# Load environment variables
load_dotenv()

def test_groq_provider_implementation():
    """Test that the Groq provider can be instantiated and used"""
    print("\n" + "="*50)
    print("Testing Groq Provider Implementation")
    print("="*50)
    
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
        # Test 1: Create a Groq provider instance
        print("\nTest 1: Creating Groq provider instance...")
        provider = GroqProvider()
        print("✅ Successfully created Groq provider instance")
        
        # Test 2: Create a simple request
        print("\nTest 2: Creating a simple request...")
        request = LLMRequest(
            user_prompt="Hello, how are you?",
            llm_config=LLMModelConfig(
                provider=ProviderType.GROQ,
                model_name="llama3-8b-8192"
            )
        )
        print("✅ Successfully created request")
        
        # Test 3: Generate a response
        print("\nTest 3: Generating a response...")
        response = provider.generate(request)
        print("✅ Successfully generated response")
        print(f"\nResponse content (first 100 chars): {response.content[:100]}...")
        print(f"Model used: {response.model}")
        print(f"Token usage: {response.usage}")
        
        return True
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_groq_provider_implementation()
    if success:
        print("\n✅ All tests passed - Groq provider implementation is working correctly")
    else:
        print("\n❌ Tests failed - Groq provider implementation has issues") 