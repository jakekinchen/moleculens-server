from fastapi import Depends, Query, HTTPException
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
from pydantic import BaseModel
from typing import Literal
import os

# Define the allowed model names and their configurations
class ModelNameEnum(str):
    # OpenAI models
    O3_MINI = "o3-mini"
    O1 = "o1"
    GPT_4_5 = "gpt-4.5-preview"
    GPT_4_0 = "gpt-4o"
    GPT_4 = "gpt-4"
    GPT_4_TURBO = "gpt-4-turbo"
    GPT_3_5_TURBO = "gpt-3.5-turbo"
    GPT_3_5 = "gpt-3.5"
    
    # Anthropic models
    A_3_7 = "claude-3-7-sonnet-latest"
    A_3_5 = "claude-3-5-sonnet-latest"
    CLAUDE_3_OPUS = "claude-3-opus-20240229"
    CLAUDE_3_SONNET = "claude-3-sonnet-20240229"
    CLAUDE_3_HAIKU = "claude-3-haiku-20240307"
    
    # Groq models
    GROQ_LLAMA_3 = "llama3-70b-8192"
    GROQ_QWEN = "qwen-2.5-coder-32b"
    LLAMA_2 = "llama-2"
    LLAMA2_70B = "llama2-70b-4096"
    MIXTRAL = "mixtral-8x7b-32768"
    QWEN_72B = "qwen-72b-max"

# Mapping of models to providers and environment variables
MODEL_PROVIDER_MAPPING = {
    # OpenAI models
    ModelNameEnum.O3_MINI: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.O1: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_4_5: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_4_0: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_4: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_4_TURBO: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_3_5_TURBO: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_3_5: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    
    # Anthropic models
    ModelNameEnum.A_3_7: (ProviderType.ANTHROPIC, "ANTHROPIC_API_KEY"),
    ModelNameEnum.A_3_5: (ProviderType.ANTHROPIC, "ANTHROPIC_API_KEY"),
    ModelNameEnum.CLAUDE_3_OPUS: (ProviderType.ANTHROPIC, "ANTHROPIC_API_KEY"),
    ModelNameEnum.CLAUDE_3_SONNET: (ProviderType.ANTHROPIC, "ANTHROPIC_API_KEY"),
    ModelNameEnum.CLAUDE_3_HAIKU: (ProviderType.ANTHROPIC, "ANTHROPIC_API_KEY"),
    
    # Groq models
    ModelNameEnum.GROQ_LLAMA_3: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.GROQ_QWEN: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.LLAMA_2: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.LLAMA2_70B: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.MIXTRAL: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.QWEN_72B: (ProviderType.GROQ, "GROQ_API_KEY"),
}

# List of all supported model names
ALLOWED_MODELS = list(MODEL_PROVIDER_MAPPING.keys())

class ModelValidator(BaseModel):
    model_name: str
    
    # Add configuration to avoid model_name warning
    model_config = {
        "protected_namespaces": ()
    }
    
    def validate_model(self):
        if self.model_name not in ALLOWED_MODELS:
            raise ValueError(f"Invalid model name '{self.model_name}'. Allowed values: {', '.join(ALLOWED_MODELS)}")
        return True

def use_llm(model_name: str = Query("o3-mini", alias="model")):
    """
    Dependency that initializes LLMService based on the selected model.
    Ensures model_name is valid and selects the appropriate provider and API key.
    """
    # Validate the model name
    validator = ModelValidator(model_name=model_name)
    try:
        validator.validate_model()
    except ValueError:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid model name '{model_name}'. Allowed values: {', '.join(ALLOWED_MODELS)}"
        )
    
    # Get provider and API key environment variable
    provider, api_key_env = MODEL_PROVIDER_MAPPING[model_name]
    api_key = os.getenv(api_key_env)
    
    # Check if we have the API key
    if not api_key:
        raise HTTPException(
            status_code=500,
            detail=f"Missing API key for provider {provider.value}. Please set the {api_key_env} environment variable."
        )
    
    # Create config and service
    llm_config = LLMModelConfig(
        provider=provider,
        model_name=model_name,
        api_key=api_key
    )
    
    try:
        # Create and return the LLM service
        return LLMService(llm_config)
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error initializing LLM service: {str(e)}"
        )