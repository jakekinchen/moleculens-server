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
) -> T:
    client = get_client()  # singleton

    messages: list[ChatCompletionMessageParam] = cast(
        "list[ChatCompletionMessageParam]",
        [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
    )

    if ResponseFormatJSONSchema is not None:
        # Use structured output for newer OpenAI versions
        completion = client.chat.completions.create(
            model=model_name,
            messages=messages,
            response_format=ResponseFormatJSONSchema(
                type="json_schema",
                schema=response_model.model_json_schema(),
                strict=True,
            ),
        )
    else:
        # Fallback for older OpenAI versions
        completion = client.chat.completions.create(
            model=model_name,
            messages=messages,
            response_format={"type": "json_object"},
        )

    content = completion.choices[0].message.content
    if content is None:
        raise ValueError("No content received from OpenAI")
    return response_model.model_validate_json(content)
