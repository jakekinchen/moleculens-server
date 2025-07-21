"""
LLM Service Module - Handles interactions with LLM provider OpenAI
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Generic, Literal, Optional, Type, TypeVar, Union

from dotenv import load_dotenv
from pydantic import BaseModel

# Load environment variables from .env file
load_dotenv()

# Type for the response model
T = TypeVar("T", bound=BaseModel)


class ProviderType(str, Enum):
    """Supported LLM providers."""

    OPENAI = "openai"


class MessageRole(str, Enum):
    """Standardized message roles across providers."""

    SYSTEM = "system"
    USER = "user"
    ASSISTANT = "assistant"


class Message(BaseModel):
    """Standardized message format."""

    role: MessageRole
    content: str


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


class LLMRequest(BaseModel):
    """Request to generate text from an LLM."""

    user_prompt: str
    llm_config: Optional[LLMModelConfig] = None
    system_prompt: Optional[str] = None
    assistant_prompt: Optional[str] = None
    function_prompt: Optional[str] = None
    tool_prompt: Optional[str] = None
    response_format: Optional[Dict[str, Any]] = None
    tools: Optional[list[Dict[str, Any]]] = None
    functions: Optional[list[Dict[str, Any]]] = None
    tool_choice: Optional[str] = None
    function_call: Optional[str] = None
    seed: Optional[int] = None
    response_format_type: Optional[str] = None
    response_format_schema: Optional[Dict[str, Any]] = None
    max_retries: Optional[int] = None
    timeout: Optional[float] = None
    request_timeout: Optional[float] = None
    pool_timeout: Optional[float] = None
    api_type: Optional[str] = None
    api_version: Optional[str] = None
    api_base: Optional[str] = None
    api_key: Optional[str] = None
    organization: Optional[str] = None
    engine: Optional[str] = None
    deployment_id: Optional[str] = None
    model: Optional[str] = None
    temperature: Optional[float] = None
    max_tokens: Optional[int] = None
    top_p: Optional[float] = None
    frequency_penalty: Optional[float] = None
    presence_penalty: Optional[float] = None
    stop: Optional[list[str]] = None
    logit_bias: Optional[Dict[str, float]] = None
    user: Optional[str] = None
    n: Optional[int] = None
    stream: Optional[bool] = None
    echo: Optional[bool] = None
    best_of: Optional[int] = None
    logprobs: Optional[int] = None
    suffix: Optional[str] = None
    prompt: Optional[str] = None
    context: Optional[str] = None
    examples: Optional[list[Dict[str, Any]]] = None
    labels: Optional[list[str]] = None
    temperature_multiplier: Optional[float] = None
    top_p_multiplier: Optional[float] = None
    frequency_penalty_multiplier: Optional[float] = None
    presence_penalty_multiplier: Optional[float] = None
    max_tokens_multiplier: Optional[float] = None
    stop_multiplier: Optional[float] = None
    logit_bias_multiplier: Optional[float] = None
    user_multiplier: Optional[float] = None
    n_multiplier: Optional[float] = None
    best_of_multiplier: Optional[float] = None
    logprobs_multiplier: Optional[float] = None
    suffix_multiplier: Optional[float] = None
    prompt_multiplier: Optional[float] = None
    context_multiplier: Optional[float] = None
    examples_multiplier: Optional[float] = None
    labels_multiplier: Optional[float] = None


class ResponseFormat(BaseModel):
    """Response format configuration."""

    type: Literal["text", "json_object"] = "text"


class StructuredLLMRequest(LLMRequest, Generic[T]):
    """Request structure for structured output."""

    response_model: Type[T]
    response_format: ResponseFormat = ResponseFormat(type="json_object")


@dataclass
class LLMResponse:
    """Standardized response from any LLM provider."""

    content: str
    model: str
    usage: Dict[str, int]


class LLMProvider(ABC):
    """Abstract base class for LLM providers."""

    @abstractmethod
    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response from the LLM."""

    @abstractmethod
    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response from the LLM."""


class LLMService:
    """Main service class for interacting with LLM providers."""

    def __init__(self, config: Optional[LLMModelConfig] = None):
        # Provide a sensible default configuration so callers (including tests)
        # can simply instantiate `LLMService()` without arguments.
        if config is None:
            config = LLMModelConfig(
                provider=ProviderType.OPENAI, model_name="gpt-3.5-turbo"
            )
        self.config = config
        self._provider = self._create_provider(config)

    def _create_provider(self, config: LLMModelConfig) -> "LLMProvider":
        """Create the appropriate provider based on configuration."""
        if config.provider == ProviderType.OPENAI:
            from .providers.openai_provider import OpenAIProvider

            return OpenAIProvider(api_key=config.api_key)
        else:
            raise ValueError(f"Unsupported provider type: {config.provider}")

    def generate(self, request: Union[str, "LLMRequest"]) -> LLMResponse:
        """Generate a response from the LLM."""
        if isinstance(request, str):
            request = LLMRequest(user_prompt=request, llm_config=self.config)
        elif request.llm_config is None:
            request.llm_config = self.config

        return self._provider.generate(request)

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response from the LLM."""
        if request.llm_config is None:
            request.llm_config = self.config
        return self._provider.generate_structured(request)


# Example usage with ThreeGroup from models.py
from api.agent_management.models import ThreeGroup


def generate_three_group(
    service: LLMService, description: str, provider: ProviderType, model_name: str
) -> ThreeGroup:
    """Generate a Three.js group from a description using specified provider.

    Args:
        service: LLMService instance
        description: Natural language description of the 3D object to generate
        provider: Which LLM provider to use
        model_name: The specific model name to use

    Returns:
        ThreeGroup: A validated Three.js group containing the described object
    """
    prompt = f"""Create a Three.js group that represents the following object: {description}

The group should:
- Have a descriptive name
- Contain one or more meshes that form the object
- Use appropriate materials and geometries
- Position and scale components correctly

You must return a valid ThreeGroup structure with all required fields. The response must follow this schema:

- name (string): A descriptive name for the group
- position (Vector3): x, y, z coordinates for group position
- rotation (Vector3): x, y, z rotation values in radians
- scale (Vector3): x, y, z scale factors (default 1.0)
- children (array of Mesh objects): Each mesh must have:
  - name (string): Unique identifier
  - geometry:
    - type (string): Must end in 'Geometry'
    - parameters (dict): Size/shape parameters
  - material:
    - type (string): Must end in 'Material'
    - color (int): Hex color value
    - opacity (float): 0.0 to 1.0
    - transparent (bool): Whether material is transparent
    - shininess (float, optional): For PhongMaterial"""

    request = StructuredLLMRequest(
        user_prompt=prompt,
        system_prompt="You are a 3D modeling expert. Generate precise Three.js structures based on descriptions.",
        response_model=ThreeGroup,
        llm_config=LLMModelConfig(provider=provider, model_name=model_name),
        temperature=0.2,  # Low temperature for more deterministic output
    )

    return service.generate_structured(request)


# Similarly expose OpenAIProvider directly for extensions that expect it
try:
    from .providers.openai_provider import OpenAIProvider  # noqa: F401
except Exception:

    class OpenAIProvider:  # type: ignore
        def __init__(self, *args, **kwargs):
            raise ImportError(
                "OpenAIProvider requires the `openai` package, which is not installed."
            )
