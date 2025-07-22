"""OpenAI provider implementation for LLM service."""

from typing import Any, Dict, List, Optional, cast

from openai.types.chat import ChatCompletion, ChatCompletionMessageParam

try:
    from ...utils.openai_client import get_client
except ImportError:
    # Fallback for direct module loading
    import os
    import sys

    sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
    from utils.openai_client import get_client

from ..llm_service import LLMProvider, LLMRequest, LLMResponse, StructuredLLMRequest, T


class OpenAIProvider(LLMProvider):
    """OpenAI-specific implementation."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize the OpenAI provider with API key."""
        self.client = get_client(api_key)

    def _convert_messages(
        self, request: LLMRequest
    ) -> List[ChatCompletionMessageParam]:
        """Convert our message format to OpenAI's format."""
        messages: List[ChatCompletionMessageParam] = []
        if request.system_prompt:
            messages.append({"role": "system", "content": request.system_prompt})
        messages.append({"role": "user", "content": request.user_prompt})
        return messages

    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response using OpenAI's API."""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)

            # Build parameters dictionary with proper typing
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "stream": request.stream,
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

            if (
                not completion.choices
                or not completion.choices[0].message
                or not completion.choices[0].message.content
            ):
                raise ValueError("No response content received from OpenAI")

            return LLMResponse(
                content=completion.choices[0].message.content,
                model=request.llm_config.model_name,
                usage={
                    "prompt_tokens": (
                        completion.usage.prompt_tokens if completion.usage else 0
                    ),
                    "completion_tokens": (
                        completion.usage.completion_tokens if completion.usage else 0
                    ),
                    "total_tokens": (
                        completion.usage.total_tokens if completion.usage else 0
                    ),
                },
            )
        except Exception as e:
            raise Exception(f"OpenAI API error: {str(e)}")

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using OpenAI's Responses API with structured outputs."""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)

            # Build parameters dictionary for the Responses API
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "input": messages,
                "text_format": request.response_model,
            }

            # Only add parameters if the model supports them (o* models don't support these)
            if not request.llm_config.model_name.startswith("o"):
                if request.temperature is not None:
                    params["temperature"] = request.temperature
                if request.max_tokens is not None:
                    params["max_output_tokens"] = (
                        request.max_tokens
                    )  # Note: different parameter name
                if request.top_p is not None:
                    params["top_p"] = request.top_p

            # Use the Responses API with structured outputs
            response = self.client.responses.parse(**params)

            # Handle different response statuses
            if response.status == "incomplete":
                if response.incomplete_details.reason == "max_output_tokens":
                    raise ValueError("Response was incomplete due to max tokens limit")
                elif response.incomplete_details.reason == "content_filter":
                    raise ValueError("Response was incomplete due to content filtering")
                else:
                    raise ValueError(
                        f"Response was incomplete: {response.incomplete_details.reason}"
                    )

            # Check for refusal
            if (
                response.output
                and len(response.output) > 0
                and response.output[0].content
                and len(response.output[0].content) > 0
                and response.output[0].content[0].type == "refusal"
            ):
                raise ValueError(
                    f"Model refused to respond: {response.output[0].content[0].refusal}"
                )

            # Return the parsed object
            if response.output_parsed:
                return cast(T, response.output_parsed)
            else:
                raise ValueError("No parsed output received from the API")

        except Exception as e:
            raise Exception(
                f"OpenAI structured output error (using Responses API): {str(e)}"
            )
