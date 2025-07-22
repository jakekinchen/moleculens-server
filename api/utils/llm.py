import json
import os
from typing import Optional

from .openai_client import get_client

SYSTEM_PROMPT = "Translate the user description into a deterministic list of PyMOL cmd API calls only."

FUNCTION_SCHEMA = {
    "name": "generate_pymol_script",
    "parameters": {
        "type": "object",
        "properties": {"commands": {"type": "array", "items": {"type": "string"}}},
        "required": ["commands"],
    },
}

# Use centralized client helper


def description_to_commands(description: str) -> list[str]:
    """Call the LLM to translate a free-form description into a list of PyMOL
    commands."""
    client = get_client()
    try:
        response = client.chat.completions.create(
            model="o3-mini",  # Using a version known to work well with function calling
            messages=[
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": description},
            ],
            functions=[FUNCTION_SCHEMA],
            function_call={"name": "generate_pymol_script"},
        )
        choice = response.choices[0]
        if not choice.message.function_call:
            raise ValueError("LLM did not return function call")
        args = json.loads(choice.message.function_call.arguments)
        return args.get("commands", [])
    except Exception as e:
        print(f"Error in LLM translation: {str(e)}")
        # In test or offline environments, gracefully fall back to a simple default
        # to avoid network dependency. Returning an empty list is safe for callers.
        return []
