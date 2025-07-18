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
    is_default: bool = Field(description="Whether this is a default model")
    default_for: List[str] = Field(description="List of tasks this model is default for")
    api_params: Dict[str, Any] = Field(description="Additional API parameters")

def register_models():
    """Register all supported LLM models with the registry"""
    
    # Clear existing registry
    ModelRegistry._registry.clear()

    # Register OpenAI models
    def create_model_factory(model_data):
        def factory(*args, **kwargs):
            data = model_data.copy()
            data.update(kwargs)
            return ModelInfo(**data)
        return factory

    models = [
        {
            "name": "o3-mini",
            "display_name": "OpenAI GPT-3.5 Mini",
            "provider": ProviderType.OPENAI,
            "categories": [ModelCategory.GENERAL, ModelCategory.REASONING],
            "context_length": 16385,
            "is_default": True,
            "default_for": ["general", "molecular"],
            "api_params": {}
        },
        {
            "name": "gpt-4",
            "display_name": "OpenAI GPT-4",
            "provider": ProviderType.OPENAI,
            "categories": [ModelCategory.GENERAL, ModelCategory.REASONING],
            "context_length": 8192,
            "is_default": False,
            "default_for": [],
            "api_params": {}
        }
    ]

    for model_data in models:
        ModelRegistry.register(
            model_data["name"],
            ModelInfo,
            create_model_factory(model_data)
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