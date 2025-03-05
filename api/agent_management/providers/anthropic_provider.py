"""
Anthropic provider implementation for LLM service.
"""

import json
from typing import Dict, Any, TypeVar, Type, List, Optional, Union, cast
import httpx
from anthropic import Anthropic, HUMAN_PROMPT, AI_PROMPT
from ..llm_service import LLMProvider, LLMRequest, LLMResponse, StructuredLLMRequest, T, MessageRole
import time
import logging
from typing import TypeVar, Generic, Optional, Dict, Any, List, Union, Callable
import os

class AnthropicProvider(LLMProvider):
    """Anthropic-specific implementation"""
    
    MAX_RETRIES = 3  # Maximum number of retry attempts
    INITIAL_RETRY_DELAY = 1  # Initial delay in seconds
    
    def __init__(self, api_key: Optional[str] = None):
        """Initialize the Anthropic provider with API key"""
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("Anthropic API key is required")
        
        transport = httpx.HTTPTransport(retries=3)
        client = httpx.Client(transport=transport)
        self.client = Anthropic(
            api_key=self.api_key,
            http_client=client
        )
        self.logger = logging.getLogger(__name__)
        
    def _get_token_usage_from_stream(self, stream):
        """Extract token usage from stream if available"""
        if hasattr(stream, "usage"):
            return {
                "prompt_tokens": stream.usage.input_tokens,
                "completion_tokens": stream.usage.output_tokens,
                "total_tokens": stream.usage.input_tokens + stream.usage.output_tokens
            }
        else:
            # Default usage when not available
            return {
                "prompt_tokens": 0,
                "completion_tokens": 0,
                "total_tokens": 0
            }

    def _convert_messages(self, request: LLMRequest) -> List[Dict[str, str]]:
        """Convert our message format to Anthropic's format"""
        return [{
            "role": "user",
            "content": request.user_prompt
        }]

    def _call_with_retry(self, func: Callable, *args, **kwargs) -> Any:
        """
        Call a function with retry logic and exponential backoff
        
        Args:
            func: The function to call
            *args: Arguments to pass to the function
            **kwargs: Keyword arguments to pass to the function
            
        Returns:
            The result of the function call
            
        Raises:
            Exception: If all retry attempts fail
        """
        retry_count = 0
        last_exception = None
        
        while retry_count < self.MAX_RETRIES:
            try:
                return func(*args, **kwargs)
            except Exception as e:
                last_exception = e
                error_message = str(e).lower()
                
                # Only retry on specific errors that might be temporary
                if "overloaded" in error_message or "rate limit" in error_message or "timeout" in error_message:
                    retry_count += 1
                    if retry_count < self.MAX_RETRIES:
                        # Calculate exponential backoff delay: 1s, 2s, 4s, etc.
                        delay = self.INITIAL_RETRY_DELAY * (2 ** (retry_count - 1))
                        self.logger.warning(f"Anthropic API error: {str(e)}. Retrying in {delay}s (attempt {retry_count}/{self.MAX_RETRIES})")
                        time.sleep(delay)
                    else:
                        self.logger.error(f"Anthropic API error: {str(e)}. Max retries ({self.MAX_RETRIES}) exceeded.")
                else:
                    # For other errors, don't retry
                    break
        
        # If we get here, all retries failed or the error wasn't retryable
        raise last_exception

    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response using Anthropic's API"""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)
            
            # Build parameters dictionary with proper typing
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "max_tokens": request.max_tokens or 4096  # Use a smaller default max_tokens
            }
            
            # Add system prompt as a top-level parameter if provided
            if request.system_prompt:
                params["system"] = request.system_prompt
            
            # Add optional parameters
            if request.temperature is not None:
                params["temperature"] = request.temperature
            if request.top_p is not None:
                params["top_p"] = request.top_p
            
            # Use streaming for large token requests
            if params["max_tokens"] > 4096:
                content_parts = []
                with self.client.messages.stream(**params) as stream:
                    for message in stream:
                        if message.type == "content_block":
                            content_parts.append(message.content[0].text)
                        elif message.type == "content_block_delta":
                            if hasattr(message.delta, "text"):
                                content_parts.append(message.delta.text)
                
                if not content_parts:
                    raise ValueError("No response content received from Anthropic")
                
                content = "".join(content_parts)
                return LLMResponse(
                    content=content,
                    model=request.llm_config.model_name,
                    usage=self._get_token_usage_from_stream(stream)
                )
            else:
                # Use non-streaming for shorter responses
                response = self._call_with_retry(self.client.messages.create, **params)
                
                if not response.content or not response.content[0].text:
                    raise ValueError("No response content received from Anthropic")
                
                return LLMResponse(
                    content=response.content[0].text,
                    model=request.llm_config.model_name,
                    usage=self._get_token_usage_from_stream(response)
                )
        except Exception as e:
            self.logger.error(f"Anthropic API error: {str(e)}")
            raise Exception(f"Anthropic API error: {str(e)}")

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using Anthropic's API"""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            # Add JSON structure requirement to the system prompt
            system_prompt = (request.system_prompt or "") + "\nYou must respond with valid JSON that matches the required schema. Do not include any explanatory text, only output the JSON object."
            
            messages = self._convert_messages(LLMRequest(
                system_prompt=None,  # We'll pass this as a top-level parameter
                user_prompt=request.user_prompt,
                llm_config=request.llm_config
            ))
            
            # Build parameters dictionary with proper typing
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "system": system_prompt,  # Pass system prompt as top-level parameter
                "max_tokens": request.max_tokens or 20000  # Use a smaller default max_tokens
            }
            
            # Add optional parameters
            if request.temperature is not None:
                params["temperature"] = request.temperature
            if request.top_p is not None:
                params["top_p"] = request.top_p
            
            # Use streaming for large token requests
            if params["max_tokens"] > 8192:
                content_parts = []
                with self.client.messages.stream(**params) as stream:
                    for message in stream:
                        if message.type == "content_block":
                            content_parts.append(message.content[0].text)
                        elif message.type == "content_block_delta":
                            if hasattr(message.delta, "text"):
                                content_parts.append(message.delta.text)
                
                if not content_parts:
                    raise ValueError("No response content received from Anthropic")
                
                content = "".join(content_parts)
            else:
                # Use non-streaming for shorter responses
                response = self.client.messages.create(**params)
                
                if not response.content or not response.content[0].text:
                    raise ValueError("No response content received from Anthropic")
                
                content = response.content[0].text
            
            if not content:
                raise ValueError("Empty response content from Anthropic")
            
            # Clean up the response to extract JSON
            content = content.strip()
            
            # Sometimes Anthropic might wrap the JSON in markdown code blocks
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0].strip()
            elif "```" in content:
                content = content.split("```")[1].strip()
            
            # Try to find JSON object if there's additional text
            try:
                # Find the first '{' and last '}'
                start = content.find('{')
                end = content.rfind('}')
                if start != -1 and end != -1:
                    content = content[start:end + 1]
            except:
                pass
                
            try:
                json_response = json.loads(content)
                return request.response_model.model_validate(json_response)
            except json.JSONDecodeError as json_err:
                raise ValueError(f"Failed to parse JSON response: {str(json_err)}. Response content: {content}")
        except Exception as e:
            raise Exception(f"Anthropic structured output error: {str(e)}") 