# Groq Model Additions

This document summarizes the models added to the Groq provider.

## Models Added

The following models have been added to the Groq provider:

### DeepSeek R1 Distill Llama 70B

- **Model ID**: `deepseek-r1-distill-llama-70b`
- **Display Name**: DeepSeek R1 Distill Llama 70B
- **Context Length**: 128,000 tokens
- **Categories**: General, Reasoning, Code
- **Description**: A distilled version of DeepSeek's R1 model, fine-tuned from the Llama-3.3-70B-Instruct base model. This model leverages knowledge distillation to retain robust reasoning capabilities while enhancing efficiency. It demonstrates strong performance across various benchmarks, particularly in mathematical reasoning and coding tasks.

### Llama 3.2 90B Vision

- **Model ID**: `llama-3.2-90b-vision-preview`
- **Display Name**: Llama 3.2 90B Vision
- **Context Length**: 128,000 tokens
- **Categories**: General, Vision
- **Description**: A multimodal model that can understand and interpret visual data from images. This model can be used for tasks such as visual question answering, caption generation, and image understanding.

### Qwen 2.5 32B

- **Model ID**: `qwen-2.5-32b`
- **Display Name**: Qwen 2.5 32B
- **Context Length**: 128,000 tokens
- **Categories**: General, Code
- **Description**: Qwen 2.5 offers increased knowledge plus boosts for coding and mathematics, big improvements in instruction following, understanding structured data (e.g., tables), and generating structured outputs especially JSON. It also has more resiliency to diverse system prompts, enhancing role-play implementation and condition-setting for chatbots.

### Llama 3.3 70B SpecDec

- **Model ID**: `llama-3.3-70b-specdec`
- **Display Name**: Llama 3.3 70B SpecDec
- **Context Length**: 8,192 tokens
- **Categories**: General, Reasoning
- **Description**: A specialized version of Llama 3.3 70B that uses speculative decoding to improve inference speed. This model is particularly useful for reasoning tasks that require step-by-step analysis, logical deduction, and structured thinking.

## Usage

These models can be used with the Groq provider by specifying the model ID in the `LLMModelConfig`:

```python
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType

# Create LLM config
llm_config = LLMModelConfig(
    provider=ProviderType.GROQ,
    model_name="deepseek-r1-distill-llama-70b",  # Use any of the model IDs listed above
    api_key=os.getenv("GROQ_API_KEY")
)

# Create LLM service
llm_service = LLMService(config=llm_config)
```

## Notes

- These models are available on the Groq platform and can be used with the Groq provider.
- The Groq API key must be set in the environment variable `GROQ_API_KEY`.
- Some models may be in preview mode and may not be available for production use.
- The context length and capabilities of each model may change as the models are updated by Groq.
- For the most up-to-date information, refer to the Groq API documentation. 