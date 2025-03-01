"""
LLM Service Module - Handles interactions with various LLM providers (OpenAI, Anthropic, Groq).
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, TypeVar, Type, Any, Union, Generic, Literal
import os
import json
from enum import Enum
from dataclasses import dataclass
from dotenv import load_dotenv
from pydantic import BaseModel, Field

# Load environment variables from .env file
load_dotenv()

# Type for the response model
T = TypeVar('T', bound=BaseModel)

class MessageRole(str, Enum):
    """Standardized message roles across providers"""
    SYSTEM = "system"
    USER = "user"
    ASSISTANT = "assistant"

class Message(BaseModel):
    """Standardized message format"""
    role: MessageRole
    content: str

class ProviderType(str, Enum):
    """Supported LLM providers"""
    OPENAI = "openai"
    ANTHROPIC = "anthropic"
    GROQ = "groq"

class LLMModelConfig(BaseModel):
    """Configuration for an LLM model"""
    provider: ProviderType
    model_name: str
    api_key: Optional[str] = None
    
    # Address warning about model_name field
    model_config = {
        "use_enum_values": True,
        "protected_namespaces": ()  # Remove protection for the 'model_' namespace
    }

class LLMRequest(BaseModel):
    """Standard request structure for all LLM providers"""
    system_prompt: Optional[str] = None
    user_prompt: str
    llm_config: Optional[LLMModelConfig] = None
    temperature: float = 0.7
    max_tokens: Optional[int] = None
    top_p: Optional[float] = None
    stream: bool = False
    additional_params: Dict[str, Any] = Field(default_factory=dict)

class ResponseFormat(BaseModel):
    """Response format configuration"""
    type: Literal["text", "json_object"] = "text"

class StructuredLLMRequest(LLMRequest, Generic[T]):
    """Request structure for structured output"""
    response_model: Type[T]
    response_format: ResponseFormat = ResponseFormat(type="json_object")

@dataclass
class LLMResponse:
    """Standardized response from any LLM provider"""
    content: str
    model: str
    usage: Dict[str, int]

class LLMProvider(ABC):
    """Abstract base class for LLM providers"""
    
    @abstractmethod
    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response from the LLM"""
        pass

    @abstractmethod
    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response from the LLM"""
        pass

class LLMService:
    """Main service class for interacting with LLM providers"""
    
    def __init__(self, config: LLMModelConfig):
        self.config = config
        self._provider = self._create_provider(config)

    def _create_provider(self, config: LLMModelConfig) -> LLMProvider:
        """Create the appropriate provider based on configuration"""
        if config.provider == ProviderType.OPENAI:
            from .providers.openai_provider import OpenAIProvider
            return OpenAIProvider(config.api_key)
        elif config.provider == ProviderType.ANTHROPIC:
            from .providers.anthropic_provider import AnthropicProvider
            return AnthropicProvider(config.api_key)
        elif config.provider == ProviderType.GROQ:
            from .providers.groq_provider import GroqProvider
            return GroqProvider(config.api_key)
        else:
            raise ValueError(f"Unsupported provider type: {config.provider}")

    def generate(self, request: Union[str, LLMRequest]) -> LLMResponse:
        """Generate a response from the LLM"""
        if isinstance(request, str):
            request = LLMRequest(
                user_prompt=request,
                llm_config=self.config
            )
        elif request.llm_config is None:
            request.llm_config = self.config
            
        return self._provider.generate(request)

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response from the LLM"""
        if request.llm_config is None:
            request.llm_config = self.config
        return self._provider.generate_structured(request)

# Example usage with ThreeGroup from models.py 
from agent_management.models import ThreeGroup

def generate_three_group(service: LLMService, description: str, provider: ProviderType, model_name: str) -> ThreeGroup:
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
        llm_config=LLMModelConfig(
            provider=provider,
            model_name=model_name
        ),
        temperature=0.2  # Low temperature for more deterministic output
    )

    return service.generate_structured(request)