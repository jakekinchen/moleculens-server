from fastapi import Depends, Query, HTTPException
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
from pydantic import BaseModel
from typing import Literal
import os

# Define the allowed model names
class ModelNameEnum(str):
    O3_MINI = "o3-mini"
    GPT_4 = "gpt-4"
    GPT_3_5 = "gpt-3.5-turbo"
    LLAMA_2 = "llama-2"

class ModelValidator(BaseModel):
    model_name: Literal[
        ModelNameEnum.O3_MINI, 
        ModelNameEnum.GPT_4, 
        ModelNameEnum.GPT_3_5, 
        ModelNameEnum.LLAMA_2
    ]

def use_llm(model_name: str = Query("o3-mini", alias="model")):
    """
    Dependency that initializes LLMService based on the model query parameter.
    Ensures model_name is within the approved list.
    """
    # Validate the model_name using Pydantic
    try:
        ModelValidator(model_name=model_name)
    except ValueError:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid model name '{model_name}'. Allowed values: {', '.join(ModelNameEnum.__dict__.values())}"
        )

    llm_config = LLMModelConfig(
        provider=ProviderType.OPENAI,
        model_name=model_name,
        api_key=os.getenv("OPENAI_API_KEY")
    )
    return LLMService(llm_config)
