"""
OpenAI provider implementation for LLM service.
"""

import json
from typing import Dict, Any, TypeVar, Type, List, Optional, Union, cast, Literal
from openai.types.chat import ChatCompletion, ChatCompletionMessage, ChatCompletionMessageParam
from openai.types.chat.chat_completion import Choice
from openai.types import Completion
from ..llm_service import LLMProvider, LLMRequest, LLMResponse, StructuredLLMRequest, T, MessageRole
import os

class OpenAIProvider(LLMProvider):
    """OpenAI-specific implementation"""
    
    def __init__(self, api_key: Optional[str] = None):
        """Initialize the OpenAI provider with API key"""
        import openai
        self.api_key = api_key or os.environ.get("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError("OpenAI API key is required")
            
        self.client = openai.OpenAI(
            api_key=self.api_key
        )

    def _convert_messages(self, request: LLMRequest) -> List[ChatCompletionMessageParam]:
        """Convert our message format to OpenAI's format"""
        messages: List[ChatCompletionMessageParam] = []
        if request.system_prompt:
            messages.append({
                "role": "system",
                "content": request.system_prompt
            })
        messages.append({
            "role": "user",
            "content": request.user_prompt
        })
        return messages

    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response using OpenAI's API"""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)
            
            # Build parameters dictionary with proper typing
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "stream": request.stream
            }
            
            # Only add parameters if the model supports them (o* models don't support these)
            if not request.llm_config.model_name.startswith("o"):
                params["temperature"] = request.temperature
                if request.max_tokens is not None:
                    params["max_tokens"] = request.max_tokens
                if request.top_p is not None:
                    params["top_p"] = request.top_p
            
            response = self.client.chat.completions.create(**params)
            
            # Cast the response to ChatCompletion since we know it's not a stream
            completion = cast(ChatCompletion, response)
            
            if not completion.choices or not completion.choices[0].message or not completion.choices[0].message.content:
                raise ValueError("No response content received from OpenAI")
            
            return LLMResponse(
                content=completion.choices[0].message.content,
                model=request.llm_config.model_name,
                usage={
                    "prompt_tokens": completion.usage.prompt_tokens if completion.usage else 0,
                    "completion_tokens": completion.usage.completion_tokens if completion.usage else 0,
                    "total_tokens": completion.usage.total_tokens if completion.usage else 0
                }
            )
        except Exception as e:
            raise Exception(f"OpenAI API error: {str(e)}")

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using OpenAI's API"""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)
            
            # Build parameters dictionary with proper typing
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "response_format": {"type": "json_object"},
                "stream": request.stream
            }
            
            # Only add parameters if the model supports them (o* models don't support these)
            if not request.llm_config.model_name.startswith("o"):
                params["temperature"] = request.temperature
                if request.max_tokens is not None:
                    params["max_tokens"] = request.max_tokens
                if request.top_p is not None:
                    params["top_p"] = request.top_p
                
            response = self.client.chat.completions.create(**params)
            
            # Cast the response to ChatCompletion since we know it's not a stream
            completion = cast(ChatCompletion, response)
            
            if not completion.choices or not completion.choices[0].message or not completion.choices[0].message.content:
                raise ValueError("No response content received from OpenAI")
            
            content = completion.choices[0].message.content
            if not content:
                raise ValueError("Empty response content from OpenAI")
                
            json_response = json.loads(content)
            return request.response_model.model_validate(json_response)
        except Exception as e:
            raise Exception(f"OpenAI structured output error: {str(e)}") 