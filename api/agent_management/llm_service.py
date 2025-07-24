"""
LLM Service Module - Handles interactions with LLM provider OpenAI
"""

from enum import Enum
from typing import Any, Dict, Generic, Optional, Type, TypeVar

# Load environment variables from .env file
from dotenv import load_dotenv
from pydantic import BaseModel

load_dotenv()

T = TypeVar("T", bound=BaseModel)


class ProviderType(str, Enum):
    """Supported LLM providers."""

    OPENAI = "openai"


class LLMModelConfig(BaseModel):
    """Configuration for an LLM model."""

    provider: ProviderType = ProviderType.OPENAI
    model_name: str
    api_key: Optional[str] = None

    # Address warning about model_name field
    model_config = {
        "use_enum_values": True,
        "protected_namespaces": (),  # Remove protection for the 'model_' namespace
    }


class LLMService(Generic[T]):
    """Main service class for interacting with LLM providers."""

    @classmethod
    def __class_getitem__(cls, item):
        # Support subscription syntax for tests
        return cls

    def __init__(self, config: Optional[LLMModelConfig] = None):
        # Provide a sensible default configuration
        if config is None:
            config = LLMModelConfig(provider=ProviderType.OPENAI, model_name="o3-mini")
        self.config = config
        self._provider = self._create_provider(config)

    def _create_provider(self, config: LLMModelConfig) -> Any:
        """Create the appropriate provider based on configuration."""
        if config.provider == ProviderType.OPENAI:
            # Use the openai_provider module directly (provides generate_structured)
            from .providers import openai_provider

            return openai_provider
        else:
            raise ValueError(f"Unsupported provider type: {config.provider}")

    def generate_structured(
        self,
        user_prompt: str,
        response_model: Type[T],
        system_prompt: str = "You are a helpful assistant.",
        model_name: Optional[str] = None,
    ) -> T:
        """Generate a structured response from the LLM."""
        return self._provider.generate_structured(
            user_prompt=user_prompt,
            response_model=response_model,
            system_prompt=system_prompt,
            model_name=model_name or self.config.model_name,
        )
