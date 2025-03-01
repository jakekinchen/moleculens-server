"""
Extension to the LLM Service for image processing capabilities
"""

from typing import Optional, List, Dict, Any, Union
import base64
from pydantic import BaseModel, Field

from llm_service import (
    LLMProvider, 
    LLMResponse, 
    ProviderType,
    LLMRequest,
    LLMModelConfig
)

# Extended message content types for multimodal models
class TextContent(BaseModel):
    type: str = "text"
    text: str

class ImageUrlContent(BaseModel):
    type: str = "image_url"
    image_url: Dict[str, str]

# Union type for content
ContentItem = Union[TextContent, ImageUrlContent]

class MultiModalMessage(BaseModel):
    role: str
    content: List[ContentItem]

class ImageToTextRequest(LLMRequest):
    """Extended request structure for image-to-text tasks"""
    image_path: str
    image_description_prompt: str = "Describe what you see in this image in detail."
    
    def get_base64_image(self) -> Optional[str]:
        """Read an image file and convert it to base64 encoding"""
        try:
            with open(self.image_path, "rb") as image_file:
                return base64.b64encode(image_file.read()).decode('utf-8')
        except Exception as e:
            raise ValueError(f"Error reading image file {self.image_path}: {e}")

def extend_llm_service(original_module):
    """Add image processing capabilities to the LLM service"""
    
    # Extend the GroqProvider for image handling
    original_groq_provider = original_module.GroqProvider
    
    class ExtendedGroqProvider(original_groq_provider):
        """Extended Groq provider with image processing capabilities"""
        
        def process_image(self, request: ImageToTextRequest) -> LLMResponse:
            """Process an image and return a text description"""
            try:
                # Get base64 encoded image
                base64_image = request.get_base64_image()
                
                # For Groq, we need to incorporate system message into the user prompt
                # because Groq doesn't support system messages with image prompts
                user_prompt = request.image_description_prompt
                if request.system_prompt:
                    user_prompt = f"{request.system_prompt}\n\n{user_prompt}"
                
                # Format messages for vision model
                messages = [
                    {
                        "role": "user",
                        "content": [
                            {"type": "text", "text": user_prompt},
                            {
                                "type": "image_url",
                                "image_url": {
                                    "url": f"data:image/png;base64,{base64_image}"
                                }
                            }
                        ]
                    }
                ]
                
                # Call the model
                response = self.client.chat.completions.create(
                    model=request.llm_config.model_name,
                    messages=messages,
                    temperature=request.temperature,
                    max_tokens=request.max_tokens or 1024,
                    **request.additional_params
                )
                
                return LLMResponse(
                    content=response.choices[0].message.content,
                    model=request.llm_config.model_name,
                    usage={
                        "prompt_tokens": response.usage.prompt_tokens if hasattr(response.usage, "prompt_tokens") else 0,
                        "completion_tokens": response.usage.completion_tokens if hasattr(response.usage, "completion_tokens") else 0,
                        "total_tokens": response.usage.total_tokens if hasattr(response.usage, "total_tokens") else 0
                    }
                )
            except Exception as e:
                raise Exception(f"Groq image processing error: {str(e)}")
    
    # Extend the OpenAIProvider for image handling
    original_openai_provider = original_module.OpenAIProvider
    
    class ExtendedOpenAIProvider(original_openai_provider):
        """Extended OpenAI provider with image processing capabilities"""
        
        def process_image(self, request: ImageToTextRequest) -> LLMResponse:
            """Process an image and return a text description"""
            try:
                # Get base64 encoded image
                base64_image = request.get_base64_image()
                
                # Format messages for vision model (e.g. gpt-4-vision)
                messages = [
                    {
                        "role": "user",
                        "content": [
                            {"type": "text", "text": request.image_description_prompt},
                            {
                                "type": "image_url",
                                "image_url": {
                                    "url": f"data:image/png;base64,{base64_image}"
                                }
                            }
                        ]
                    }
                ]
                
                # Add system message if provided
                if request.system_prompt:
                    messages.insert(0, {
                        "role": "system",
                        "content": request.system_prompt
                    })
                
                # Call the model (must be a vision-capable model)
                response = self.client.chat.completions.create(
                    model=request.llm_config.model_name,  # Should be gpt-4-vision or similar
                    messages=messages,
                    temperature=request.temperature,
                    max_tokens=request.max_tokens or 1024,
                    **request.additional_params
                )
                
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
                raise Exception(f"OpenAI image processing error: {str(e)}")
    
    # Replace the original providers with extended ones
    original_module.GroqProvider = ExtendedGroqProvider
    original_module.OpenAIProvider = ExtendedOpenAIProvider
    
    # Extend the LLMService with image processing
    original_llm_service = original_module.LLMService
    
    class ExtendedLLMService(original_llm_service):
        """Extended LLM service with image processing capabilities"""
        
        def process_image(self, request: ImageToTextRequest) -> str:
            """Process an image and return a text description"""
            provider = self._get_provider(
                request.llm_config.provider,
                request.llm_config.api_key
            )
            
            # Check if provider has image processing capability
            if hasattr(provider, 'process_image'):
                response = provider.process_image(request)
                return response.content
            else:
                raise ValueError(f"Provider {request.model_config.provider} does not support image processing")
    
    # Replace the original LLMService with the extended one
    original_module.LLMService = ExtendedLLMService
    
    # Add the new request type to the module
    original_module.ImageToTextRequest = ImageToTextRequest
    
    return original_module