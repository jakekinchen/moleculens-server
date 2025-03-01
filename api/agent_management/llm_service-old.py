"""
LLM Service Module - Handles interactions with various LLM providers.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, TypeVar, Type, cast
import os
import json
from dataclasses import dataclass
from dotenv import load_dotenv
from pydantic import BaseModel
from models import ThreeGroup

# Load environment variables from .env file
load_dotenv()

@dataclass
class LLMResponse:
    """Standardized response from any LLM provider"""
    content: str
    model: str
    usage: Dict[str, int]

class LLMProvider(ABC):
    """Abstract base class for LLM providers"""
    
    @abstractmethod
    def generate(self, prompt: str, **kwargs) -> LLMResponse:
        """Generate a response from the LLM"""
        pass

    @abstractmethod
    def generate_structured(self, prompt: str, response_model: Type[BaseModel], **kwargs) -> BaseModel:
        """Generate a structured response from the LLM"""
        pass

class OpenAIProvider(LLMProvider):
    """OpenAI-specific implementation"""
    
    def __init__(self, api_key: Optional[str] = None, model: str = "o3-mini"):
        import openai
        self.client = openai.OpenAI(
            api_key=api_key or os.getenv("OPENAI_API_KEY")
        )
        self.model = model

    def generate(self, prompt: str, **kwargs) -> LLMResponse:
        """Generate a response using OpenAI's API"""
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[{"role": "user", "content": prompt}],
                **kwargs
            )
            
            return LLMResponse(
                content=response.choices[0].message.content,
                model=self.model,
                usage={
                    "prompt_tokens": response.usage.prompt_tokens,
                    "completion_tokens": response.usage.completion_tokens,
                    "total_tokens": response.usage.total_tokens
                }
            )
        except Exception as e:
            raise Exception(f"OpenAI API error: {str(e)}")

    def generate_structured(self, prompt: str, response_model: Type[BaseModel], **kwargs) -> BaseModel:
        """Generate a structured response using OpenAI's API"""
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[{"role": "user", "content": prompt}],
                response_format={"type": "json_object"},
                **kwargs
            )
            
            json_response = json.loads(response.choices[0].message.content)
            return response_model.model_validate(json_response)
        except Exception as e:
            raise Exception(f"OpenAI structured output error: {str(e)}")

class LLMService:
    """Main service class for handling LLM interactions"""
    
    def __init__(self, provider: Optional[LLMProvider] = None):
        """Initialize with a specific provider, defaults to OpenAI"""
        self.provider = provider or OpenAIProvider()
    
    def generate(self, prompt: str, **kwargs) -> str:
        """Generate a response using the configured provider"""
        response = self.provider.generate(prompt, **kwargs)
        return response.content

    def generate_structured(self, prompt: str, response_model: Type[BaseModel], **kwargs) -> BaseModel:
        """Generate a structured response using the configured provider"""
        return self.provider.generate_structured(prompt, response_model, **kwargs)

    def generate_three_group(self, description: str, **kwargs) -> ThreeGroup:
        """Generate a Three.js group from a description.
        
        Args:
            description: Natural language description of the 3D object to generate
            **kwargs: Additional arguments to pass to the LLM
            
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

        result = self.generate_structured(prompt, ThreeGroup, **kwargs)
        return cast(ThreeGroup, result) 