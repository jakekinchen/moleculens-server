"""OpenAI provider implementation for LLM service."""

import os
from typing import Any, Dict, List, Optional, cast

from openai.types.chat import ChatCompletion, ChatCompletionMessageParam

from ..llm_service import LLMProvider, LLMRequest, LLMResponse, StructuredLLMRequest, T


class OpenAIProvider(LLMProvider):
    """OpenAI-specific implementation."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize the OpenAI provider with API key."""
        import openai

        self.api_key = api_key or os.environ.get("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError("OpenAI API key is required")

        self.client = openai.OpenAI(api_key=self.api_key)

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
        """Generate a structured response using OpenAI's API with
        client.beta.chat.completions.parse."""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)

            # Parameters for the .parse() method might differ slightly or some might be implicit.
            # Temperature, max_tokens, top_p are typically part of the create/parse call.
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "response_format": request.response_model,  # Pass Pydantic model directly
                # stream is not typically used with .parse() as it expects a full response to parse.
                # If streaming with parsing is needed, SDK might have other utilities or it needs custom handling.
            }

            # Add temperature, max_tokens, top_p if applicable and supported by .parse()
            # (Consult SDK docs for .parse() specific parameters)
            if not request.llm_config.model_name.startswith("o"):
                params["temperature"] = request.temperature
                if request.max_tokens is not None:
                    params["max_tokens"] = (
                        request.max_tokens
                    )  # Check if .parse() supports this
                if request.top_p is not None:
                    params["top_p"] = request.top_p  # Check if .parse() supports this

            # Using client.beta.chat.completions.parse
            completion_parse_result = self.client.beta.chat.completions.parse(**params)

            # The .parse() method should return a result where the parsed Pydantic object is accessible.
            # Based on search results, it might be in: completion.choices[0].message.parsed
            if (
                not completion_parse_result.choices
                or not completion_parse_result.choices[0].message
                or not hasattr(completion_parse_result.choices[0].message, "parsed")
                or completion_parse_result.choices[0].message.parsed is None
            ):
                raise ValueError(
                    "No parsed structured content received from OpenAI using .parse()"
                )

            # The .parsed attribute should already be an instance of request.response_model (T)
            parsed_object = completion_parse_result.choices[0].message.parsed

            # Ensure it's the correct type (though .parse() should handle this)
            if not isinstance(parsed_object, request.response_model):
                raise TypeError(
                    f"Parsed object is not of type {request.response_model.__name__}. "
                    f"Got {type(parsed_object).__name__}."
                )

            return cast(
                T, parsed_object
            )  # Cast for type safety, though it should already be T

        except Exception as e:
            # Catch specific OpenAI API errors if possible for better error reporting
            # For now, a general catch with a clear message.
            raise Exception(
                f"OpenAI structured output error (using .parse()): {str(e)}"
            )
