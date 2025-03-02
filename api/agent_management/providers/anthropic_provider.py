"""
Anthropic provider implementation for LLM service.
"""

import json
from typing import Dict, Any, TypeVar, Type, List, Optional, Union, cast
import httpx
from anthropic import Anthropic, HUMAN_PROMPT, AI_PROMPT
from ..llm_service import LLMProvider, LLMRequest, LLMResponse, StructuredLLMRequest, T, MessageRole

class AnthropicProvider(LLMProvider):
    """Anthropic-specific implementation"""
    
    def __init__(self, api_key: str):
        transport = httpx.HTTPTransport(retries=3)
        client = httpx.Client(transport=transport)
        self.client = Anthropic(
            api_key=api_key,
            http_client=client
        )

    def _convert_messages(self, request: LLMRequest) -> List[Dict[str, str]]:
        """Convert our message format to Anthropic's format"""
        return [{
            "role": "user",
            "content": request.user_prompt
        }]

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
                "max_tokens": request.max_tokens or 1024
            }
            
            # Add system prompt as a top-level parameter if provided
            if request.system_prompt:
                params["system"] = request.system_prompt
            
            # Add optional parameters
            if request.temperature is not None:
                params["temperature"] = request.temperature
            if request.top_p is not None:
                params["top_p"] = request.top_p
            
            response = self.client.messages.create(**params)
            
            if not response.content or len(response.content) == 0:
                raise ValueError("No response content received from Anthropic")
            
            # Get the text content from the first content block
            content = response.content[0].text
            
            # Get token usage from response
            usage = {
                "prompt_tokens": response.usage.input_tokens,
                "completion_tokens": response.usage.output_tokens,
                "total_tokens": response.usage.input_tokens + response.usage.output_tokens
            }
            
            return LLMResponse(
                content=content,
                model=request.llm_config.model_name,
                usage=usage
            )
        except Exception as e:
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
                "max_tokens": request.max_tokens or 20480
            }
            
            # Add optional parameters
            if request.temperature is not None:
                params["temperature"] = request.temperature
            if request.top_p is not None:
                params["top_p"] = request.top_p
                
            response = self.client.messages.create(**params)
            
            if not response.content or len(response.content) == 0:
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