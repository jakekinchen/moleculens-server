import json
import os
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

client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))


async def description_to_commands(description: str) -> list[str]:
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
