from fastapi import Depends, Query, HTTPException
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
from pydantic import BaseModel
from typing import Literal
import os

# Define the allowed model names and their corresponding providers
class ModelNameEnum(str):
    O3_MINI = "o3-mini"
    O1 = "o1"
    GPT_4_5 = "gpt-4.5-preview"
    GPT_4_0 = "gpt-4o"
    A_3_7 = "claude-3-7-sonnet-latest"
    A_3_5 = "claude-3-5-sonnet-latest"
    GROQ_LLAMA_3 = "llama3-70b-8192"
    GROQ_QWEN = "qwen-2.5-coder-32b"
    LLAMA_2 = "llama-2"
    GPT_3_5 = "gpt-3.5"

# Mapping of models to providers and environment variables
MODEL_PROVIDER_MAPPING = {
    ModelNameEnum.O3_MINI: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.O1: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_4_5: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.GPT_4_0: (ProviderType.OPENAI, "OPENAI_API_KEY"),
    ModelNameEnum.A_3_7: (ProviderType.ANTHROPIC, "ANTHROPIC_API_KEY"),
    ModelNameEnum.A_3_5: (ProviderType.ANTHROPIC, "ANTHROPIC_API_KEY"),
    ModelNameEnum.GROQ_LLAMA_3: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.GROQ_QWEN: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.LLAMA_2: (ProviderType.GROQ, "GROQ_API_KEY"),
    ModelNameEnum.GPT_3_5: (ProviderType.OPENAI, "OPENAI_API_KEY"),
}

class ModelValidator(BaseModel):
    model_name: Literal[
        ModelNameEnum.O3_MINI,
        ModelNameEnum.O1,
        ModelNameEnum.GPT_4_5,
        ModelNameEnum.GPT_4_0,
        ModelNameEnum.A_3_7,
        ModelNameEnum.A_3_5,
        ModelNameEnum.GROQ_LLAMA_3,
        ModelNameEnum.GROQ_QWEN,
        ModelNameEnum.LLAMA_2,
        ModelNameEnum.GPT_3_5
    ]
    
    model_config = {"protected_namespaces": ()}

def use_llm(model_name: str = Query("o3-mini", alias="model")):
    """
    Dependency that initializes LLMService based on the selected model.
    Ensures model_name is valid and selects the appropriate provider and API key.
    """
    try:
        ModelValidator(model_name=model_name)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid model name: {model_name}")
    
    provider, api_key_env = MODEL_PROVIDER_MAPPING[model_name]
    api_key = os.getenv(api_key_env)
    
    if not api_key:
        raise HTTPException(status_code=500, detail=f"Missing API key for provider: {provider}")
    
    llm_config = LLMModelConfig(model_name=model_name, provider=provider, api_key=api_key)
    return LLMService(llm_config)
