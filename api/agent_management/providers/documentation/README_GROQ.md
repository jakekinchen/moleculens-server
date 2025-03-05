# Groq Provider for LLM Service

This document provides comprehensive instructions for using the Groq provider with the LLM service in the Scientific Visualization AI Server.

## Overview

The Groq provider allows you to use Groq's API to access models like Llama 3 and Mixtral with extremely fast inference times. Groq's infrastructure is optimized for low-latency inference, making it ideal for interactive applications.

## Setup

### 1. Install Dependencies

The Groq provider requires the `groq` Python package. It's already included in the project's requirements.txt, but you can install it manually if needed:

```bash
pip install groq
```

### 2. API Key

You need a Groq API key to use the provider. You can get one by signing up at [https://console.groq.com/](https://console.groq.com/).

Set the API key in your environment:

```bash
# Add to your .env file
GROQ_API_KEY=your_api_key_here
```

Or set it directly in your code:

```python
from api.agent_management.llm_service import LLMModelConfig, ProviderType

config = LLMModelConfig(
    provider=ProviderType.GROQ,
    model_name="llama3-70b-8192",
    api_key="your_api_key_here"
)
```

## Available Models

The Groq provider supports the following models:

| Model Name | Description | Context Length |
|------------|-------------|---------------|
| `llama3-70b-8192` | Llama 3 70B model | 8,192 tokens |
| `llama3-8b-8192` | Llama 3 8B model | 8,192 tokens |
| `mixtral-8x7b-32768` | Mixtral 8x7B model | 32,768 tokens |
| `gemma-7b-it` | Gemma 7B-IT model | 8,192 tokens |

## Usage Examples

### Basic Text Generation

```python
from api.agent_management.llm_service import LLMService, LLMRequest, LLMModelConfig, ProviderType

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

### Structured Output Generation

The Groq provider fully supports structured output generation using Pydantic models. This allows you to get responses in a structured format that can be easily processed by your application.

```python
from pydantic import BaseModel
from typing import List, Optional
from api.agent_management.llm_service import LLMService, StructuredLLMRequest, LLMModelConfig, ProviderType

# Define your output model
class Recipe(BaseModel):
    name: str
    ingredients: List[str]
    steps: List[str]
    prep_time_minutes: Optional[int] = None
    cook_time_minutes: Optional[int] = None

# Create an LLM service with Groq provider
llm_service = LLMService(
    config=LLMModelConfig(
        provider=ProviderType.GROQ,
        model_name="llama3-70b-8192"
    )
)

# Create a structured request
request = StructuredLLMRequest(
    user_prompt="Give me a recipe for a simple chocolate chip cookie.",
    system_prompt="You are a helpful cooking assistant that provides recipes in a structured format.",
    response_model=Recipe
)

# Generate a structured response
recipe = llm_service.generate_structured(request)
print(f"Recipe: {recipe.name}")
print(f"Ingredients: {recipe.ingredients}")
print(f"Steps: {recipe.steps}")
```

### Supported Response Models

The Groq provider has been tested and confirmed to work with the following response models:

1. **Simple Models**
   - `Recipe` - For generating cooking recipes
   - `BooleanResponse` - For true/false validation with confidence scores
   - `MovieReview` - For generating structured movie reviews
   - `Vector3` - For generating 3D coordinates

2. **Complex Models**
   - `ThreeGroup` - For generating Three.js 3D scene groups
   - `SceneObject` - For generating scene objects with properties and relationships
   - And other nested models used throughout the application

## Advanced Configuration

### Temperature

Control the randomness of the output by setting the temperature:

```python
request = LLMRequest(
    user_prompt="Write a creative story about a robot.",
    temperature=0.9  # Higher temperature for more creative output
)
```

### Max Tokens

Limit the length of the response:

```python
request = LLMRequest(
    user_prompt="Summarize the history of artificial intelligence.",
    max_tokens=500  # Limit response to 500 tokens
)
```

### Top-p Sampling

Control the diversity of the output:

```python
request = LLMRequest(
    user_prompt="Generate ideas for a science project.",
    top_p=0.8  # Use nucleus sampling with p=0.8
)
```

## Troubleshooting

### API Key Issues

If you get an authentication error, check that your API key is correct and properly set in your environment.

### Model Not Found

If you get a model not found error, check that you're using one of the supported model names listed above.

### Package Not Installed

If you get an error about the Groq package not being installed, run `pip install groq`.

### JSON Parsing Errors

The Groq provider includes robust JSON parsing for structured output. If you encounter JSON parsing errors:

1. Check that your Pydantic model matches the expected output structure
2. Try using a lower temperature setting (e.g., 0.2) for more deterministic output
3. Make sure your prompt clearly specifies the expected format

### Rate Limits

If you hit rate limits, consider implementing retry logic or reducing the frequency of requests.

## Performance Considerations

- Groq's API is optimized for low-latency inference, making it ideal for interactive applications.
- The Llama 3 70B model provides the best quality but may be more expensive.
- For cost-sensitive applications, consider using the Llama 3 8B model.
- For long context applications, consider using the Mixtral model with its 32K context window.

## Testing

You can test the Groq provider using the provided test scripts:

```bash
# Test basic text generation
python api/test_groq_simple.py

# Test using the LLMService
python api/test_groq_service.py

# Test structured output
python api/test_groq_structured.py

# Test all response models
python api/test_groq_all_models.py
``` 