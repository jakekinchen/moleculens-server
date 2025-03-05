"""
Test script for the Groq provider implementation
"""

import os
from dotenv import load_dotenv
from api.agent_management.llm_service import LLMService, LLMRequest, LLMModelConfig, ProviderType

# Load environment variables
load_dotenv()

def test_groq_provider():
    """Test the Groq provider with a simple text generation task"""
    print("Testing Groq Provider...")
    
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
    
    # Create an LLM service with Groq provider
    llm_service = LLMService(
        config=LLMModelConfig(
            provider=ProviderType.GROQ,
            model_name="llama3-70b-8192",
            api_key=os.getenv("GROQ_API_KEY")
        )
    )
    
    # Create a simple request
    request = LLMRequest(
        user_prompt="Explain the difference between a design pattern and an algorithm in 2-3 sentences.",
        system_prompt="You are a helpful assistant that provides clear, concise explanations."
    )
    
    try:
        # Generate a response
        response = llm_service.generate(request)
        
        # Print the response
        print("\nGroq LLama-3 Response:")
        print("-" * 50)
        print(response.content)
        print("-" * 50)
        print(f"Model: {response.model}")
        print(f"Tokens: {response.usage}")
        
        return True
    except Exception as e:
        print(f"Error testing Groq provider: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_groq_provider()
    if success:
        print("\n✅ Groq provider test passed")
    else:
        print("\n❌ Groq provider test failed")