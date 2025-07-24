"""
Direct LLM utility using OpenAI provider without the complex agent factory pattern.
"""

from typing import Optional, Type, TypeVar

from pydantic import BaseModel

from api.agent_management.llm_service import (
    LLMModelConfig,
    ProviderType,
    StructuredLLMRequest,
)
from api.agent_management.providers.openai_provider import OpenAIProvider

T = TypeVar("T", bound=BaseModel)


def generate_structured(
    user_prompt: str,
    response_model: Type[T],
    system_prompt: Optional[str] = None,
    model_name: str = "gpt-4o-mini",
    temperature: Optional[float] = None,
    max_tokens: Optional[int] = None,
    api_key: Optional[str] = None,
) -> T:
    """
    Generate structured output using OpenAI provider directly.

    Args:
        user_prompt: The user's input prompt
        response_model: Pydantic model class for structured output
        system_prompt: Optional system prompt
        model_name: OpenAI model name (default: o3-mini)
        temperature: Temperature for generation (ignored for o* models)
        max_tokens: Maximum tokens to generate
        api_key: Optional API key (uses env var if not provided)

    Returns:
        Instance of response_model with generated data

    Raises:
        Exception: If generation fails
    """
    # Create provider
    provider = OpenAIProvider(api_key=api_key)

    # Create config
    config = LLMModelConfig(
        provider=ProviderType.OPENAI, model_name=model_name, api_key=api_key
    )

    # Create request
    request = StructuredLLMRequest(
        user_prompt=user_prompt,
        system_prompt=system_prompt,
        response_model=response_model,
        llm_config=config,
        temperature=temperature,
        max_tokens=max_tokens,
    )

    # Generate structured response
    return provider.generate_structured(request)


def generate_text(
    user_prompt: str,
    system_prompt: Optional[str] = None,
    model_name: str = "gpt-4o-mini",
    temperature: Optional[float] = None,
    max_tokens: Optional[int] = None,
    api_key: Optional[str] = None,
) -> str:
    """
    Generate simple text response using OpenAI provider directly.

    Args:
        user_prompt: The user's input prompt
        system_prompt: Optional system prompt
        model_name: OpenAI model name (default: o3-mini)
        temperature: Temperature for generation (ignored for o* models)
        max_tokens: Maximum tokens to generate
        api_key: Optional API key (uses env var if not provided)

    Returns:
        Generated text response

    Raises:
        Exception: If generation fails
    """

    # Simple text response model
    class TextResponse(BaseModel):
        content: str

    result = generate_structured(
        user_prompt=user_prompt,
        response_model=TextResponse,
        system_prompt=system_prompt,
        model_name=model_name,
        temperature=temperature,
        max_tokens=max_tokens,
        api_key=api_key,
    )

    return result.content
