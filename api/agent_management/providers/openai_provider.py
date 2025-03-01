"""
OpenAI provider implementation for LLM service.
"""

import json
from typing import Dict, Any, TypeVar, Type, List, Optional, Union, cast
from openai.types.chat import ChatCompletion, ChatCompletionMessage, ChatCompletionMessageParam
from openai.types.chat.chat_completion import Choice
from openai.types import Completion
from ..llm_service import LLMProvider, LLMRequest, LLMResponse, StructuredLLMRequest, T, MessageRole
from openai.types.chat.chat_completion_message_param import ChatCompletionMessageParam
from openai.types.chat.chat_completion_message import ChatCompletionMessage
from openai.types.chat.chat_completion_message_param import ChatCompletionMessageParam
from dotenv import load_dotenv

class OpenAIProvider(LLMProvider):
    """OpenAI-specific implementation"""
    
    def __init__(self, api_key: str):
        import openai
        self.client = openai.OpenAI(
            api_key=load_dotenv("OPENAI_API_KEY")
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
            
            # Build parameters dictionary
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "stream": request.stream
            }
            
            # Only add temperature if it's supported by the model (o3 models don't support it)
            if not request.llm_config.model_name.startswith("o3-"):
                params["temperature"] = request.temperature
            
            # Only add max_tokens if it's not None (omit it entirely for o3 models)
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
                
            # Add top_p if provided (o3 models don't support it)
            if request.top_p is not None and not request.llm_config.model_name.startswith("o3-"):
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
            
            # Build parameters dictionary
            params = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "response_format": {"type": "json_object"},
                "stream": request.stream
            }
            
            # Only add temperature if it's supported by the model (o3 models don't support it)
            if not request.llm_config.model_name.startswith("o3-"):
                params["temperature"] = request.temperature
            
            # Only add max_tokens if it's not None (omit it entirely for o3 models)
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
                
            # Add top_p if provided (o3 models don't support it)
            if request.top_p is not None and not request.llm_config.model_name.startswith("o3-"):
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