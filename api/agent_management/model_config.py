"""
Model configuration module for registering LLM providers and models.
"""

from typing import Dict, Any, List, Callable, Optional
from enum import Enum, auto
from pydantic import BaseModel, Field
from agent_management.models import ModelRegistry
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType

class ModelCategory(str, Enum):
    """Categories of LLM models by capability"""
    GENERAL = "general"
    CODE = "code"
    VISION = "vision"
    REASONING = "reasoning"

class ModelInfo(BaseModel):
    """Information about a specific LLM model"""
    name: str = Field(description="Full model name including version")
    display_name: str = Field(description="User-friendly display name")
    provider: ProviderType = Field(description="The LLM provider")
    categories: List[ModelCategory] = Field(description="Model capabilities")
    context_length: int = Field(description="Maximum context length in tokens")
    is_default: bool = Field(default=False, description="Whether this is a default model")
    default_for: List[str] = Field(default_factory=list, description="List of use cases for which this model is default (e.g. 'geometry', 'molecular')")
    api_params: Dict[str, Any] = Field(default_factory=dict, description="Provider-specific API parameters")

# Register LLM models
def register_models():
    """Register all supported LLM models with the registry"""
    
    ModelRegistry.register(
        "o1",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "o1"),
            display_name=kwargs.get("display_name", "o1"),
            provider=kwargs.get("provider", ProviderType.OPENAI),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING]),
            context_length=kwargs.get("context_length", 128000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )

    ModelRegistry.register(
        "o3-mini",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "o3-mini"),
            display_name=kwargs.get("display_name", "o3-mini"),
            provider=kwargs.get("provider", ProviderType.OPENAI),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING]),
            context_length=kwargs.get("context_length", 128000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )

    ModelRegistry.register(
        "gpt-4.5-preview",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "gpt-4.5-preview"),
            display_name=kwargs.get("display_name", "GPT-4.5 Preview"),
            provider=kwargs.get("provider", ProviderType.OPENAI),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING, ModelCategory.CODE]),
            context_length=kwargs.get("context_length", 128000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    
    ModelRegistry.register(
        "gpt-4o",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "gpt-4o"),
            display_name=kwargs.get("display_name", "GPT-4o"),
            provider=kwargs.get("provider", ProviderType.OPENAI),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING, ModelCategory.VISION]),
            context_length=kwargs.get("context_length", 128000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    # Anthropic models
    ModelRegistry.register(
        "claude-3-7-sonnet-latest",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "claude-3-7-sonnet-latest"),
            display_name=kwargs.get("display_name", "Claude 3.7 Sonnet"),
            provider=kwargs.get("provider", ProviderType.ANTHROPIC),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING, ModelCategory.CODE]),
            context_length=kwargs.get("context_length", 200000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    ModelRegistry.register(
        "claude-3-5-sonnet-latest",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "claude-3-5-sonnet-latest"),
            display_name=kwargs.get("display_name", "Claude 3.5 Sonnet"),
            provider=kwargs.get("provider", ProviderType.ANTHROPIC),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING]),
            context_length=kwargs.get("context_length", 200000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    # Groq models - Use the correct model names as shown in Groq API docs
    ModelRegistry.register(
        "llama3-70b-8192",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "llama3-70b-8192"),
            display_name=kwargs.get("display_name", "Llama 3 70B"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL]),
            context_length=kwargs.get("context_length", 8192),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    ModelRegistry.register(
        "llama3-8b-8192",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "llama3-8b-8192"),
            display_name=kwargs.get("display_name", "Llama 3 8B"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL]),
            context_length=kwargs.get("context_length", 8192),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    ModelRegistry.register(
        "mixtral-8x7b-32768",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "mixtral-8x7b-32768"),
            display_name=kwargs.get("display_name", "Mixtral 8x7B"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL]),
            context_length=kwargs.get("context_length", 32768),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    ModelRegistry.register(
        "gemma-7b-it",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "gemma-7b-it"),
            display_name=kwargs.get("display_name", "Gemma 7B-IT"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL]),
            context_length=kwargs.get("context_length", 8192),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    # New Groq models
    ModelRegistry.register(
        "deepseek-r1-distill-llama-70b",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "deepseek-r1-distill-llama-70b"),
            display_name=kwargs.get("display_name", "DeepSeek R1 Distill Llama 70B"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING, ModelCategory.CODE]),
            context_length=kwargs.get("context_length", 128000),
            is_default=kwargs.get("is_default", False),
            default_for=kwargs.get("default_for", ["geometry"]),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    ModelRegistry.register(
        "llama-3.2-90b-vision-preview",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "llama-3.2-90b-vision-preview"),
            display_name=kwargs.get("display_name", "Llama 3.2 90B Vision"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.VISION]),
            context_length=kwargs.get("context_length", 128000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    ModelRegistry.register(
        "qwen-2.5-32b",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "qwen-2.5-32b"),
            display_name=kwargs.get("display_name", "Qwen 2.5 32B"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.CODE]),
            context_length=kwargs.get("context_length", 128000),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )
    
    ModelRegistry.register(
        "llama-3.3-70b-specdec",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "llama-3.3-70b-specdec"),
            display_name=kwargs.get("display_name", "Llama 3.3 70B SpecDec"),
            provider=kwargs.get("provider", ProviderType.GROQ),
            categories=kwargs.get("categories", [ModelCategory.GENERAL, ModelCategory.REASONING]),
            context_length=kwargs.get("context_length", 8192),
            is_default=kwargs.get("is_default", False),
            api_params=kwargs.get("api_params", {})
        )
    )

# Helper functions for working with models
def get_llm_service(model_name: str) -> LLMService:
    """
    Create an LLM service instance for a specific model.
    
    Args:
        model_name: The registered model name
        
    Returns:
        An initialized LLMService instance
        
    Raises:
        ValueError: If the model is not registered
    """
    model_info = ModelRegistry.create_instance(model_name)
    
    # Create LLM config from model info
    llm_config = LLMModelConfig(
        provider=model_info.provider,
        model_name=model_info.name,
        api_key=None  # API key will be loaded from environment variables
    )
    
    # Create and return the LLM service
    return LLMService(config=llm_config)

def get_default_model_for_use_case(use_case: str) -> str:
    """
    Get the name of the default model for a specific use case.
    
    Args:
        use_case: The use case to get the default model for (e.g. 'geometry', 'molecular')
        
    Returns:
        The name of the default model for the use case
        
    Raises:
        ValueError: If no default model is found for the use case
    """
    for name in ModelRegistry.list_models():
        model_info = ModelRegistry.create_instance(name)
        if use_case in model_info.default_for:
            return name
    
    # If no default is set for this use case, return the global default or first model
    return get_default_model()

def get_default_model() -> str:
    """Get the name of the default model"""
    # First check for global default
    for name in ModelRegistry.list_models():
        model_info = ModelRegistry.create_instance(name)
        if model_info.is_default:
            return name
    
    # If no default is set, return the first model
    models = ModelRegistry.list_models()
    if models:
        return models[0]
    else:
        raise ValueError("No models registered")

def get_models_by_category(category: ModelCategory) -> List[str]:
    """
    Get all models that support a specific category.
    
    Args:
        category: The model category to filter by
        
    Returns:
        A list of model names that support the category
    """
    matching_models = []
    
    for name in ModelRegistry.list_models():
        model_info = ModelRegistry.create_instance(name)
        if category in model_info.categories:
            matching_models.append(name)
    
    return matching_models

# Initialize the model registry
register_models()