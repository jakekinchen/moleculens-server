"""
Test script to verify the Groq provider registration
"""

import os
import sys
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.llm_service import LLMService, LLMModelConfig, ProviderType

# Load environment variables
load_dotenv()

def test_groq_provider_registration():
    """Test that the Groq provider is properly registered in the provider system"""
    print("\n" + "="*50)
    print("Testing Groq Provider Registration")
    print("="*50)
    
    # Check if GROQ_API_KEY is set
    if not os.getenv("GROQ_API_KEY"):
        print("Error: GROQ_API_KEY environment variable is not set")
        return False
    
    try:
        # Test 1: Check that GroqProvider is imported in __init__.py
        print("\nTest 1: Checking provider import in __init__.py...")
        try:
            from api.agent_management.providers import GroqProvider
            print("✅ GroqProvider is properly imported in __init__.py")
        except ImportError:
            print("❌ GroqProvider is not properly imported in __init__.py")
            return False
        
        # Test 2: Check that LLMService can create a Groq provider
        print("\nTest 2: Checking LLMService can create a Groq provider...")
        config = LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name="llama3-8b-8192"
        )
        
        llm_service = LLMService(config=config)
        
        # Check that the provider is a GroqProvider
        from api.agent_management.providers.groq_provider import GroqProvider
        if isinstance(llm_service._provider, GroqProvider):
            print("✅ LLMService successfully created a GroqProvider instance")
        else:
            print(f"❌ LLMService created a {type(llm_service._provider).__name__} instead of GroqProvider")
            return False
        
        # Test 3: Check that the provider can generate a response
        print("\nTest 3: Checking the provider can generate a response...")
        response = llm_service.generate("Hello, how are you?")
        if response and response.content:
            print("✅ Provider successfully generated a response")
            print(f"\nResponse content (first 100 chars): {response.content[:100]}...")
        else:
            print("❌ Provider failed to generate a response")
            return False
        
        return True
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_groq_provider_registration()
    if success:
        print("\n✅ All tests passed - Groq provider registration is working correctly")
    else:
        print("\n❌ Tests failed - Groq provider registration has issues") 