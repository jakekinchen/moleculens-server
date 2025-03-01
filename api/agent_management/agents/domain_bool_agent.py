import os
import base64
from dotenv import load_dotenv
from pydantic import BaseModel
from typing import Optional

# Import our refactored LLM service
from llm_service import (
    LLMService, 
    LLMRequest, 
    StructuredLLMRequest, 
    LLMModelConfig, 
    ProviderType
)



def test_openai_structured():
    """Test OpenAI's o3-mini model with structured boolean output"""
    print("\n===== Testing OpenAI    with structured boolean output =====")
    
    llm_service = LLMService()
    
    # Create a request for a simple true/false classification
    request = StructuredLLMRequest(
        user_prompt="""Is the following query scientific or not? '{prmopt}' 
        
Respond with a JSON object that has this exact structure:
{
  "is_true": false,
  "confidence": 0.9
}

Where:
- is_true: boolean (true or false)
- confidence: number between 0 and 1
- reasoning: string explaining your reasoning""",
        system_prompt="You are a fact-checking AI. Determine if statements are true or false. Always respond with JSON exactly matching the requested format.",
        llm_config=LLMModelConfig(
            provider=ProviderType.OPENAI,
            model_name="o3-mini"  # Using a model that supports temperature parameter
        ),
        response_model=BooleanResponse,
    )
    
    try:
        # Get the structured response
        response = llm_service.generate_structured(request)
        
        # Print the results
        print(f"Statement: 'Python is a compiled programming language.'")
        print(f"Response: {response.is_true}")
        print(f"Confidence: {response.confidence}")
        print(f"Reasoning: {response.reasoning}")
        
        return True
    except Exception as e:
        print(f"Error testing OpenAI structured output: {e}")
        return False
