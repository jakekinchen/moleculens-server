"""Model configuration and registry."""

from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel

from api.agent_management.llm_service import LLMService
from api.agent_management.model_registry import ModelRegistry


class ProviderType(str, Enum):
    """Supported LLM providers."""

    OPENAI = "openai"


class ModelCategory(str, Enum):
    """Categories of models based on their capabilities."""

    GENERAL = "general"
    CODE = "code"
    REASONING = "reasoning"
    VISION = "vision"


class LLMModelConfig(BaseModel):
    """Configuration for an LLM model."""

    provider: ProviderType
    model_name: str
    api_key: Optional[str] = None
    api_base: Optional[str] = None
    api_version: Optional[str] = None
    organization: Optional[str] = None
    category: Optional[ModelCategory] = None
    max_tokens: Optional[int] = None
    temperature: Optional[float] = None
    top_p: Optional[float] = None
    frequency_penalty: Optional[float] = None
    presence_penalty: Optional[float] = None
    stop: Optional[List[str]] = None
    logit_bias: Optional[Dict[str, float]] = None
    user: Optional[str] = None
    default_system_prompt: Optional[str] = None
    default_user_prompt: Optional[str] = None
    default_assistant_prompt: Optional[str] = None
    default_function_prompt: Optional[str] = None
    default_tool_prompt: Optional[str] = None
    default_response_format: Optional[Dict[str, Any]] = None
    default_tools: Optional[List[Dict[str, Any]]] = None
    default_functions: Optional[List[Dict[str, Any]]] = None
    default_tool_choice: Optional[str] = None
    default_function_call: Optional[str] = None
    default_seed: Optional[int] = None
    default_response_format_type: Optional[str] = None
    default_response_format_schema: Optional[Dict[str, Any]] = None
    default_max_retries: Optional[int] = None
    default_timeout: Optional[float] = None
    default_request_timeout: Optional[float] = None
    default_pool_timeout: Optional[float] = None
    default_api_type: Optional[str] = None
    default_api_version: Optional[str] = None
    default_engine: Optional[str] = None
    default_deployment_id: Optional[str] = None
    default_organization: Optional[str] = None
    default_api_base: Optional[str] = None
    default_api_key: Optional[str] = None
    default_model: Optional[str] = None
    default_temperature: Optional[float] = None
    default_max_tokens: Optional[int] = None
    default_top_p: Optional[float] = None
    default_frequency_penalty: Optional[float] = None
    default_presence_penalty: Optional[float] = None
    default_stop: Optional[List[str]] = None
    default_logit_bias: Optional[Dict[str, float]] = None
    default_user: Optional[str] = None
    default_n: Optional[int] = None
    default_stream: Optional[bool] = None
    default_echo: Optional[bool] = None
    default_best_of: Optional[int] = None
    default_logprobs: Optional[int] = None
    default_suffix: Optional[str] = None
    default_prompt: Optional[str] = None
    default_context: Optional[str] = None
    default_examples: Optional[List[Dict[str, Any]]] = None
    default_labels: Optional[List[str]] = None
    default_temperature_multiplier: Optional[float] = None
    default_top_p_multiplier: Optional[float] = None
    default_frequency_penalty_multiplier: Optional[float] = None
    default_presence_penalty_multiplier: Optional[float] = None
    default_max_tokens_multiplier: Optional[float] = None
    default_stop_multiplier: Optional[float] = None
    default_logit_bias_multiplier: Optional[float] = None
    default_user_multiplier: Optional[float] = None
    default_n_multiplier: Optional[float] = None
    default_best_of_multiplier: Optional[float] = None
    default_logprobs_multiplier: Optional[float] = None
    default_suffix_multiplier: Optional[float] = None
    default_prompt_multiplier: Optional[float] = None
    default_context_multiplier: Optional[float] = None
    default_examples_multiplier: Optional[float] = None
    default_labels_multiplier: Optional[float] = None


def register_models() -> None:
    """Register all available models."""
    # Register OpenAI models
    ModelRegistry.register(
        "gpt-3.5-turbo",
        LLMModelConfig,
        lambda: LLMModelConfig(
            provider=ProviderType.OPENAI,
            model_name="gpt-3.5-turbo",
            category=ModelCategory.GENERAL,
        ),
    )

    ModelRegistry.register(
        "gpt-4",
        LLMModelConfig,
        lambda: LLMModelConfig(
            provider=ProviderType.OPENAI,
            model_name="gpt-4",
            category=ModelCategory.GENERAL,
        ),
    )


def get_default_model() -> str:
    """Get the default model name."""
    return "gpt-3.5-turbo"


def get_models_by_category(category: ModelCategory) -> List[str]:
    """Get all models in a category."""
    models = []
    for model_name in ModelRegistry.list_models():
        model = ModelRegistry.create_instance(model_name)
        if model.category == category:
            models.append(model_name)
    return models


def get_default_model_for_use_case(use_case: str) -> str:
    """Get the default model for a specific use case."""
    # For now, return the general default model
    # This can be expanded to have use-case specific models
    return get_default_model()


def get_llm_service(model_name: Optional[str] = None) -> LLMService:
    """Get an LLM service instance for a model."""
    import os

    from api.agent_management.llm_service import LLMModelConfig, ProviderType

    if model_name is None:
        model_name = get_default_model()
    # Use OPENAI_API_KEY from environment
    api_key = os.environ.get("OPENAI_API_KEY")
    config = LLMModelConfig(
        provider=ProviderType.OPENAI, model_name=model_name, api_key=api_key
    )
    return LLMService(config)
