from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar, cast

from pydantic import BaseModel

if TYPE_CHECKING:
    from openai.types.chat import ChatCompletionMessageParam

try:
    from openai.types.response_format import ResponseFormatJSONSchema
except ImportError:
    # Fallback for older OpenAI versions
    ResponseFormatJSONSchema = None

from .openai_client import get_client

T = TypeVar("T", bound=BaseModel)


def generate_structured(
    user_prompt: str,
    response_model: type[T],
    *,
    model_name: str = "o3-mini",
    system_prompt: str = "You are a helpful assistant.",
    max_tokens: int = 2000,
    temperature: float = 0.3,
) -> T:
    client = get_client()  # singleton

    # For o3-mini, we need to use JSON mode with explicit instructions
    enhanced_system_prompt = (
        f"{system_prompt} Respond with valid JSON that matches this schema: {response_model.model_json_schema()}"
    )
    enhanced_user_prompt = f"{user_prompt}\n\nPlease respond with valid JSON only."

    messages: list[ChatCompletionMessageParam] = cast(
        "list[ChatCompletionMessageParam]",
        [
            {"role": "system", "content": enhanced_system_prompt},
            {"role": "user", "content": enhanced_user_prompt},
        ],
    )

    try:
        # Try structured output first (for newer models)
        if ResponseFormatJSONSchema is not None and model_name not in ["o3-mini"]:
            completion = client.chat.completions.create(
                model=model_name,
                messages=messages,
                max_tokens=max_tokens,
                temperature=temperature,
                response_format=ResponseFormatJSONSchema(
                    type="json_schema",
                    schema=response_model.model_json_schema(),
                    strict=True,
                ),
            )
        else:
            # Use JSON mode for o3-mini and fallback
            completion = client.chat.completions.create(
                model=model_name,
                messages=messages,
                max_tokens=max_tokens,
                temperature=temperature,
                response_format={"type": "json_object"},
            )
    except Exception:
        # If structured output fails, try without response_format
        completion = client.chat.completions.create(
            model=model_name,
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
        )

    content = completion.choices[0].message.content
    if content is None:
        raise ValueError("No content received from OpenAI")

    # Try to parse as JSON, if it fails, try to extract JSON from the response
    try:
        return response_model.model_validate_json(content)
    except Exception:
        # Try to find JSON in the response
        import re

        # Look for JSON-like content
        json_match = re.search(r"\{.*\}", content, re.DOTALL)
        if json_match:
            json_str = json_match.group()
            return response_model.model_validate_json(json_str)

        # If no JSON found, create a simple response
        if hasattr(response_model, "content"):
            # For YAMLResponse and similar models
            return response_model(content=content)
        else:
            raise ValueError(f"Could not parse response as JSON: {content}") from None
