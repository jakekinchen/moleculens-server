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
    
    model_config = {"use_enum_values": True}

class LLMRequest(BaseModel):
    """Standard request structure for all LLM providers"""
    system_prompt: Optional[str] = None
    user_prompt: str
    llm_config: LLMModelConfig
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

class OpenAIProvider(LLMProvider):
    """OpenAI-specific implementation"""
    
    def __init__(self, api_key: str):
        import openai
        self.client = openai.OpenAI(
            api_key=api_key
        )

    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response using OpenAI's API"""
        try:
            messages = []
            if request.system_prompt:
                messages.append({"role": "system", "content": request.system_prompt})
            messages.append({"role": "user", "content": request.user_prompt})
            
            # Build parameters, excluding None values
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "temperature": request.temperature,
                "stream": request.stream,
            }
            
            # Only include non-None parameters
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
            if request.top_p is not None:
                params["top_p"] = request.top_p
                
            # Add any additional params
            params.update(request.additional_params)
            
            # Call the API
            response = self.client.chat.completions.create(**params)
            
            return LLMResponse(
                content=response.choices[0].message.content,
                model=request.llm_config.model_name,
                usage={
                    "prompt_tokens": response.usage.prompt_tokens,
                    "completion_tokens": response.usage.completion_tokens,
                    "total_tokens": response.usage.total_tokens
                }
            )
        except Exception as e:
            raise Exception(f"OpenAI API error: {str(e)}")

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using OpenAI's API"""
        try:
            messages = []
            if request.system_prompt:
                messages.append({"role": "system", "content": request.system_prompt})
            messages.append({"role": "user", "content": request.user_prompt})
            
            # Build parameters, excluding None values
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "response_format": {"type": "json_object"},
                "temperature": request.temperature,
                "stream": request.stream,
            }
            
            # Only include non-None parameters
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
            if request.top_p is not None:
                params["top_p"] = request.top_p
                
            # Add any additional params
            params.update(request.additional_params)
            
            # Call the API
            response = self.client.chat.completions.create(**params)
            
            json_response = json.loads(response.choices[0].message.content)
            return request.response_model.model_validate(json_response)
        except Exception as e:
            raise Exception(f"OpenAI structured output error: {str(e)}")

class AnthropicProvider(LLMProvider):
    """Anthropic-specific implementation"""
    
    def __init__(self, api_key: str):
        try:
            import anthropic
            import httpx
            # Create client with explicitly setting the httpx client to avoid socket_options issue
            # This helps with compatibility across different versions of anthropic and httpx
            self.client = anthropic.Anthropic(
                api_key=api_key,
                # Don't pass httpx client or additional options that may not be supported
            )
        except ImportError:
            raise ImportError("Anthropic package is not installed. Install with 'pip install anthropic'")

    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response using Anthropic's API"""
        try:
            # Build messages format for Anthropic
            messages = []
            if request.system_prompt:
                messages.append({"role": "system", "content": request.system_prompt})
            messages.append({"role": "user", "content": request.user_prompt})
            
            # Create params dict to avoid incompatible parameters
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "temperature": request.temperature,
                "max_tokens": request.max_tokens or 1024,
            }
            
            # Add any supported additional params
            supported_params = ["top_p", "top_k", "stop_sequences", "stream"]
            for param, value in request.additional_params.items():
                if param in supported_params:
                    params[param] = value
            
            response = self.client.messages.create(**params)
            
            return LLMResponse(
                content=response.content[0].text,
                model=request.llm_config.model_name,
                usage={
                    "input_tokens": response.usage.input_tokens,
                    "output_tokens": response.usage.output_tokens,
                    "total_tokens": response.usage.input_tokens + response.usage.output_tokens
                }
            )
        except Exception as e:
            raise Exception(f"Anthropic API error: {str(e)}")

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using Anthropic's API"""
        try:
            # Add specific instructions for JSON format
            json_instruction = """
            You must respond with a valid JSON object that conforms to the specified schema.
            Your entire response should be valid JSON without any additional text or explanation.
            """
            
            system_prompt = request.system_prompt or ""
            enhanced_system_prompt = system_prompt + json_instruction
            
            messages = [
                {"role": "system", "content": enhanced_system_prompt},
                {"role": "user", "content": request.user_prompt}
            ]
            
            # Create params dict to avoid incompatible parameters
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "temperature": request.temperature,
                "max_tokens": request.max_tokens or 1024,
            }
            
            # Add any supported additional params
            supported_params = ["top_p", "top_k", "stop_sequences", "stream"]
            for param, value in request.additional_params.items():
                if param in supported_params:
                    params[param] = value
            
            response = self.client.messages.create(**params)
            
            # Try to parse JSON from the response
            try:
                json_response = json.loads(response.content[0].text)
                return request.response_model.model_validate(json_response)
            except json.JSONDecodeError:
                # Attempt to extract JSON if it's wrapped in markdown or other text
                import re
                json_pattern = r'```(?:json)?\s*([\s\S]*?)\s*```'
                match = re.search(json_pattern, response.content[0].text)
                if match:
                    json_str = match.group(1)
                    json_response = json.loads(json_str)
                    return request.response_model.model_validate(json_response)
                else:
                    raise Exception("Failed to parse JSON from Anthropic response")
                
        except Exception as e:
            raise Exception(f"Anthropic structured output error: {str(e)}")

class GroqProvider(LLMProvider):
    """Groq-specific implementation"""
    
    def __init__(self, api_key: str):
        try:
            import groq
            self.client = groq.Groq(
                api_key=api_key
            )
        except ImportError:
            raise ImportError("Groq package is not installed. Install with 'pip install groq'")

    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response using Groq's API"""
        try:
            messages = []
            if request.system_prompt:
                messages.append({"role": "system", "content": request.system_prompt})
            messages.append({"role": "user", "content": request.user_prompt})
            
            # Build parameters, excluding None values
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "temperature": request.temperature,
                "stream": request.stream,
            }
            
            # Only include non-None parameters
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
            if request.top_p is not None:
                params["top_p"] = request.top_p
                
            # Add any additional params
            params.update(request.additional_params)
            
            # Call the API
            response = self.client.chat.completions.create(**params)
            
            return LLMResponse(
                content=response.choices[0].message.content,
                model=request.llm_config.model_name,
                usage={
                    "prompt_tokens": response.usage.prompt_tokens,
                    "completion_tokens": response.usage.completion_tokens,
                    "total_tokens": response.usage.total_tokens
                }
            )
        except Exception as e:
            raise Exception(f"Groq API error: {str(e)}")

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using Groq's API"""
        try:
            # Add instructions for JSON format
            json_instruction = """
            Return a valid JSON object that matches the required schema.
            Your entire response should be valid JSON without any explanations.
            """
            
            system_prompt = request.system_prompt or ""
            enhanced_system_prompt = system_prompt + json_instruction
            
            messages = [
                {"role": "system", "content": enhanced_system_prompt},
                {"role": "user", "content": request.user_prompt}
            ]
            
            # Build parameters, excluding None values
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "temperature": request.temperature,
            }
            
            # Only include non-None parameters
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
            if request.top_p is not None:
                params["top_p"] = request.top_p
                
            # Add any additional params
            params.update(request.additional_params)
            
            # Call the API
            response = self.client.chat.completions.create(**params)
            
            try:
                json_response = json.loads(response.choices[0].message.content)
                return request.response_model.model_validate(json_response)
            except json.JSONDecodeError:
                # Try to extract JSON if it's wrapped in markdown
                import re
                json_pattern = r'```(?:json)?\s*([\s\S]*?)\s*```'
                match = re.search(json_pattern, response.choices[0].message.content)
                if match:
                    json_str = match.group(1)
                    json_response = json.loads(json_str)
                    return request.response_model.model_validate(json_response)
                else:
                    raise Exception("Failed to parse JSON from Groq response")
                
        except Exception as e:
            raise Exception(f"Groq structured output error: {str(e)}")

class LLMService:
    """Main service class for handling LLM interactions"""
    
    def __init__(self):
        """Initialize the service"""
        self._providers = {}
    
    def _get_provider(self, provider_type: ProviderType, api_key: Optional[str] = None) -> LLMProvider:
        """Get or create a provider instance"""
        # Create a unique provider key that includes the API key to handle multiple API keys
        provider_key = f"{provider_type}_{api_key}"
        
        if provider_key not in self._providers:
            if provider_type == ProviderType.OPENAI:
                # Use provided API key or load from .env
                api_key = api_key or os.getenv("OPENAI_API_KEY")
                if not api_key:
                    raise ValueError("No OpenAI API key provided in request or .env file")
                self._providers[provider_key] = OpenAIProvider(api_key)
                
            elif provider_type == ProviderType.ANTHROPIC:
                api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
                if not api_key:
                    raise ValueError("No Anthropic API key provided in request or .env file")
                self._providers[provider_key] = AnthropicProvider(api_key)
                
            elif provider_type == ProviderType.GROQ:
                api_key = api_key or os.getenv("GROQ_API_KEY")
                if not api_key:
                    raise ValueError("No Groq API key provided in request or .env file")
                self._providers[provider_key] = GroqProvider(api_key)
                
            else:
                raise ValueError(f"Unsupported provider: {provider_type}")
        
        return self._providers[provider_key]
    
    def generate(self, request: LLMRequest) -> str:
        """Generate a response using the specified provider"""
        provider = self._get_provider(
            request.llm_config.provider, 
            request.llm_config.api_key
        )
        response = provider.generate(request)
        return response.content

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using the specified provider"""
        provider = self._get_provider(
            request.llm_config.provider, 
            request.llm_config.api_key
        )
        return provider.generate_structured(request)

# Example usage with ThreeGroup from models.py 
from models import ThreeGroup

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