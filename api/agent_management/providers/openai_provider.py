import os
from typing import List, Type, TypeVar, cast

from openai import OpenAI
from openai.types.chat import ChatCompletionMessageParam
from openai.types.response_format import ResponseFormatJSONSchema
from pydantic import BaseModel

T = TypeVar("T", bound=BaseModel)

client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))


def generate_structured(
    user_prompt: str,
    response_model: Type[T],
    *,
    model_name: str = "o3-mini",
    system_prompt: str = "You are a helpful assistant.",
) -> T:
    messages: List[ChatCompletionMessageParam] = cast(
        List[ChatCompletionMessageParam],
        [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
    )

    completion = client.chat.completions.create(
        model=model_name,
        messages=messages,
        response_format=ResponseFormatJSONSchema(  # <- typed param
            type="json_schema",
            schema=response_model.model_json_schema(),
            strict=True,
        ),
    )

    content = completion.choices[0].message.content
    if content is None:
        raise ValueError("No content received from OpenAI")
    return response_model.model_validate_json(content)
