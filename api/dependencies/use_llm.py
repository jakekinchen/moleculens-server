from fastapi import Depends, Query, HTTPException
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
from pydantic import BaseModel
from typing import Literal
import os

# Define the allowed model names and their providers
class ModelConfig:
    # OpenAI models
    O3_MINI = ("o3-mini", ProviderType.OPENAI)
    GPT_4 = ("gpt-4", ProviderType.OPENAI)
    GPT_4_TURBO = ("gpt-4-turbo", ProviderType.OPENAI)
    GPT_3_5_TURBO = ("gpt-3.5-turbo", ProviderType.OPENAI)
    
    # Anthropic models
    CLAUDE_3_OPUS = ("claude-3-opus-20240229", ProviderType.ANTHROPIC)
    CLAUDE_3_SONNET = ("claude-3-sonnet-20240229", ProviderType.ANTHROPIC)
    CLAUDE_3_HAIKU = ("claude-3-haiku-20240307", ProviderType.ANTHROPIC)
    
    # Groq models
    LLAMA_3 = ("llama3-70b-8192", ProviderType.GROQ)
    LLAMA_2 = ("llama2-70b-4096", ProviderType.GROQ)
    MIXTRAL = ("mixtral-8x7b-32768", ProviderType.GROQ)
    QWEN = ("qwen-72b-max", ProviderType.GROQ)

# Create a dictionary mapping model names to their provider types
MODEL_PROVIDER_MAP = {
    # OpenAI models
    "o3-mini": ProviderType.OPENAI,
    "gpt-4": ProviderType.OPENAI,
    "gpt-4-turbo": ProviderType.OPENAI,
    "gpt-3.5-turbo": ProviderType.OPENAI,
    
    # Anthropic models
    "claude-3-opus-20240229": ProviderType.ANTHROPIC,
    "claude-3-sonnet-20240229": ProviderType.ANTHROPIC, 
    "claude-3-haiku-20240307": ProviderType.ANTHROPIC,
    
    # Groq models
    "llama3-70b-8192": ProviderType.GROQ,
    "llama2-70b-4096": ProviderType.GROQ,
    "mixtral-8x7b-32768": ProviderType.GROQ,
    "qwen-72b-max": ProviderType.GROQ
}

# List of all supported model names
ALLOWED_MODELS = list(MODEL_PROVIDER_MAP.keys())

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
    Dependency that initializes LLMService based on the model query parameter.
    Selects the appropriate provider based on the model name.
    """
    # Validate the model_name
    validator = ModelValidator(model_name=model_name)
    try:
        validator.validate_model()
    except ValueError as e:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid model name '{model_name}'. Allowed values: {', '.join(ALLOWED_MODELS)}"
        )

    # Get the provider type for this model
    provider_type = MODEL_PROVIDER_MAP[model_name]
    
    # Select the appropriate API key based on provider
    api_key = None
    if provider_type == ProviderType.OPENAI:
        api_key = os.getenv("OPENAI_API_KEY")
    elif provider_type == ProviderType.ANTHROPIC:
        api_key = os.getenv("ANTHROPIC_API_KEY")
    elif provider_type == ProviderType.GROQ:
        api_key = os.getenv("GROQ_API_KEY")
    
    # Check if we have the API key
    if not api_key:
        raise HTTPException(
            status_code=500,
            detail=f"Missing API key for provider {provider_type.value}. Please set the appropriate environment variable."
        )

    # Create config with correct provider and model
    llm_config = LLMModelConfig(
        provider=provider_type,
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
