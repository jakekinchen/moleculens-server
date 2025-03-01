from fastapi import Depends, Query, HTTPException
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
from pydantic import BaseModel
from typing import Literal
import os


# openai: o1, o3-mini, gpt-4.5-preview, gpt-4o
# anthropic: claude-3-7-sonnet-latest, claude-3-5-sonnet-latest
# groq: llama3-70b-8192, qwen-2.5-coder-32b


from pydantic import BaseModel
from typing import Literal

# Define the allowed model names
class ModelNameEnum(str):
    O3_MINI = "o3-mini"
    O1 = "o1"
    GPT_4_5 = "gpt-4.5-preview"
    GPT_4 = "gpt-4"
    GPT_4_0 = "gpt-4o"
    A_3_7 = "claude-3-7-sonnet-latest"
    A_3_5 = "claude-3-5-sonnet-latest"
    GROQ_LLAMA_3 = "llama3-70b-8192"
    GROQ_QWEN = "qwen-2.5-coder-32b"
    LLAMA_2 = "llama-2"
    GPT_3_5 = "gpt-3.5"  # Adding this since it was referenced but missing

class ModelValidator(BaseModel):
    model_name: Literal[
        ModelNameEnum.O3_MINI,
        ModelNameEnum.O1,
        ModelNameEnum.GPT_4_5,
        ModelNameEnum.GPT_4,
        ModelNameEnum.GPT_4_0,
        ModelNameEnum.A_3_7,
        ModelNameEnum.A_3_5,
        ModelNameEnum.GROQ_LLAMA_3,
        ModelNameEnum.GROQ_QWEN,
        ModelNameEnum.LLAMA_2,
        ModelNameEnum.GPT_3_5  # Ensuring consistency
    ]
    
    # Add configuration to avoid model_name warning
    model_config = {
        "protected_namespaces": ()
    }


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
