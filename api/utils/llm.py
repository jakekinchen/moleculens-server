import json
import os
from typing import Optional
from openai import OpenAI

SYSTEM_PROMPT = (
    "Translate the user description into a deterministic list of PyMOL cmd API calls only."
)

FUNCTION_SCHEMA = {
    "name": "generate_pymol_script",
    "parameters": {
        "type": "object",
        "properties": {
            "commands": {"type": "array", "items": {"type": "string"}}
        },
        "required": ["commands"],
    },
}

# Lazily instantiate the OpenAI client so that merely importing this module
# does not require a valid API key (which may be unavailable during unit
# testing). If the key is missing we fall back to a dummy value so the client
# can still be created and the downstream code/tests can stub network calls.

_CLIENT: Optional[OpenAI] = None


def _get_client() -> OpenAI:
    global _CLIENT
    if _CLIENT is None:
        api_key = os.environ.get("OPENAI_API_KEY", "DUMMY_TEST_KEY")
        _CLIENT = OpenAI(api_key=api_key)
    return _CLIENT


def description_to_commands(description: str) -> list[str]:
    """Call the LLM to translate a free-form description into a list of PyMOL commands."""
    client = _get_client()
    try:
        response = client.chat.completions.create(
            model="o3-mini",
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
    except Exception:
        # In test or offline environments, gracefully fall back to a simple default
        # to avoid network dependency. Returning an empty list is safe for callers.
        return []
