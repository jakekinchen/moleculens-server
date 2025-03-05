# Using the Groq Provider

The Groq provider integrates Groq's API with our LLM service, allowing you to use models like Llama 3 and Mixtral in your applications.

## Setup

1. Install the required package:
   ```
   pip install groq
   ```

2. Set up your Groq API key:
   - Add `GROQ_API_KEY=your_api_key` to your `.env` file
   - Or pass it directly when creating the provider

## Available Models

The Groq provider supports the following models:

- `llama3-70b-8192`: Llama 3 70B model (8K context)
- `llama3-8b-8192`: Llama 3 8B model (8K context)
- `mixtral-8x7b-32768`: Mixtral 8x7B model (32K context)
- `gemma-7b-it`: Gemma 7B-IT model (8K context)

## Usage

### Simple text generation

```python
from agent_management.llm_service import LLMService, LLMRequest, LLMModelConfig, ProviderType

# Create an LLM service with Groq provider
llm_service = LLMService(
    config=LLMModelConfig(
        provider=ProviderType.GROQ,
        model_name="llama3-70b-8192"
    )
)

# Create a request
request = LLMRequest(
    user_prompt="Explain the concept of neural networks in simple terms.",
    system_prompt="You are a helpful assistant that explains complex topics simply."
)

# Generate a response
response = llm_service.generate(request)
print(response.content)
```

### Structured output generation

```python
from pydantic import BaseModel
from typing import List, Optional
from agent_management.llm_service import StructuredLLMRequest

# Define your output model
class Recipe(BaseModel):
    name: str
    ingredients: List[str]
    steps: List[str]
    prep_time: Optional[int] = None

# Create a structured request
request = StructuredLLMRequest(
    user_prompt="Give me a recipe for chocolate chip cookies.",
    system_prompt="You are a helpful cooking assistant.",
    response_model=Recipe,
    llm_config=LLMModelConfig(
        provider=ProviderType.GROQ,
        model_name="llama3-70b-8192"
    )
)

# Generate a structured response
recipe = llm_service.generate_structured(request)
print(f"Recipe: {recipe.name}")
print(f"Ingredients: {recipe.ingredients}")
```

## Troubleshooting

- If you get an error about the Groq package not being installed, run `pip install groq`
- If you get an authentication error, check that your API key is correct
- If you get a model not found error, check that you're using one of the supported model names
- If you need a different model, check the Groq API documentation for available models