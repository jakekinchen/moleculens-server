"""
Test script to verify the Groq provider API key management
"""

import os
import sys
import tempfile
from dotenv import load_dotenv

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.providers.groq_provider import GroqProvider
from api.agent_management.llm_service import LLMRequest, LLMModelConfig, ProviderType

# Load environment variables
load_dotenv()

def test_groq_api_key_management():
    """Test that the Groq provider properly handles API key management"""
    print("\n" + "="*50)
    print("Testing Groq Provider API Key Management")
    print("="*50)
    
    # Store the original API key
    original_api_key = os.getenv("GROQ_API_KEY")
    if not original_api_key:
        print("Error: GROQ_API_KEY environment variable is not set")
        return False
    
    try:
        # Test 1: Create provider with API key from environment variable
        print("\nTest 1: Creating provider with API key from environment variable...")
        provider_env = GroqProvider()
        
        # Create a simple request
        request = LLMRequest(
            user_prompt="Hello, this is a test.",
            llm_config=LLMModelConfig(
                provider=ProviderType.GROQ,
                model_name="llama3-8b-8192"
            )
        )
        
        # Generate a response
        response_env = provider_env.generate(request)
        if response_env and response_env.content:
            print("✅ Provider with environment variable API key successfully generated a response")
        else:
            print("❌ Provider with environment variable API key failed to generate a response")
            return False
        
        # Test 2: Create provider with explicit API key
        print("\nTest 2: Creating provider with explicit API key...")
        provider_explicit = GroqProvider(api_key=original_api_key)
        
        # Generate a response
        response_explicit = provider_explicit.generate(request)
        if response_explicit and response_explicit.content:
            print("✅ Provider with explicit API key successfully generated a response")
        else:
            print("❌ Provider with explicit API key failed to generate a response")
            return False
        
        # Test 3: Test error handling with missing API key
        print("\nTest 3: Testing error handling with missing API key...")
        
        # Temporarily unset the environment variable
        os.environ["GROQ_API_KEY"] = ""
        
        try:
            # This should raise an error
            provider_missing = GroqProvider()
            print("❌ Provider did not raise an error with missing API key")
            return False
        except ValueError as e:
            if "Groq API key is not provided and GROQ_API_KEY environment variable is not set" in str(e):
                print("✅ Provider correctly raised an error with missing API key")
            else:
                print(f"❌ Provider raised an unexpected error: {str(e)}")
                return False
        finally:
            # Restore the original API key
            if original_api_key:
                os.environ["GROQ_API_KEY"] = original_api_key
        
        return True
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False
    finally:
        # Ensure the original API key is restored
        if original_api_key:
            os.environ["GROQ_API_KEY"] = original_api_key

if __name__ == "__main__":
    success = test_groq_api_key_management()
    if success:
        print("\n✅ All tests passed - Groq provider API key management is working correctly")
    else:
        print("\n❌ Tests failed - Groq provider API key management has issues") 