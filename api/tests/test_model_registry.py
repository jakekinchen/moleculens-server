"""
Simple test script for the Model Registry implementation
"""

from typing import Dict, List, Optional, Type, TypeVar, Callable, Any
from pydantic import BaseModel, Field
from enum import Enum

# Define a type variable for models
T = TypeVar('T')

class ModelRegistry:
    """
    Central registry for managing model classes and their instantiation.
    Implements the Registry design pattern to decouple model references.
    """
    _registry = {}
    
    @classmethod
    def register(cls, name: str, model_cls, factory: Optional[Callable] = None):
        """
        Register a model class with an optional factory function.
        
        Args:
            name: Unique identifier for the model
            model_cls: The model class to register
            factory: Optional factory function to instantiate the model
        """
        cls._registry[name] = {
            "class": model_cls,
            "factory": factory or (lambda: model_cls())
        }
        return model_cls  # Return for decorator use
    
    @classmethod
    def get_model(cls, name: str):
        """
        Get a model class by name.
        
        Args:
            name: Name of the registered model
            
        Returns:
            The model class
            
        Raises:
            ValueError: If the model is not registered
        """
        model_data = cls._registry.get(name)
        if model_data is None:
            raise ValueError(f"Model '{name}' not registered.")
        return model_data["class"]
    
    @classmethod
    def create_instance(cls, model_name: str, *args, **kwargs):
        """
        Create an instance of a registered model.
        
        Args:
            model_name: Name of the registered model
            *args, **kwargs: Arguments to pass to the factory function
            
        Returns:
            An instance of the model
            
        Raises:
            ValueError: If the model is not registered
        """
        model_data = cls._registry.get(model_name)
        if model_data is None:
            raise ValueError(f"Model '{model_name}' not registered.")
        
        factory = model_data["factory"]
        return factory(*args, **kwargs)
    
    @classmethod
    def list_models(cls):
        """List all registered model names"""
        return list(cls._registry.keys())

# Test functionality
def test_model_registry():
    """Test the basic functionality of the model registry"""
    
    # Clear registry for testing
    ModelRegistry._registry = {}
    
    # Define a simple test model
    class TestModel(BaseModel):
        name: str
        value: int
    
    # Register model with a factory that accepts arguments
    ModelRegistry.register(
        "test-model", 
        TestModel,
        lambda *args, **kwargs: TestModel(name=kwargs.get("name", "default"), value=kwargs.get("value", 42))
    )
    
    # Get model class
    model_cls = ModelRegistry.get_model("test-model")
    assert model_cls == TestModel
    
    # Create instance
    instance = ModelRegistry.create_instance("test-model")
    assert isinstance(instance, TestModel)
    assert instance.name == "default"
    assert instance.value == 42
    
    # Create with custom parameters (passing kwargs instead of using named arguments)
    custom = ModelRegistry.create_instance("test-model", **{"name": "custom", "value": 99})
    assert custom.name == "custom"
    assert custom.value == 99
    
    print("ModelRegistry basic tests passed!")

# Define model category for demonstration
class ModelCategory(str, Enum):
    GENERAL = "general"
    CODE = "code"
    VISION = "vision"

# Simple model info class
class ModelInfo(BaseModel):
    name: str
    display_name: str
    categories: List[ModelCategory]

# Register sample models
def register_test_models():
    ModelRegistry.register(
        "model1",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "model1"),
            display_name=kwargs.get("display_name", "Model One"),
            categories=kwargs.get("categories", [ModelCategory.GENERAL])
        )
    )
    
    ModelRegistry.register(
        "model2",
        ModelInfo,
        lambda *args, **kwargs: ModelInfo(
            name=kwargs.get("name", "model2"),
            display_name=kwargs.get("display_name", "Model Two"),
            categories=kwargs.get("categories", [ModelCategory.CODE, ModelCategory.GENERAL])
        )
    )

# Test model registration
def test_model_info():
    register_test_models()
    
    models = ModelRegistry.list_models()
    assert "model1" in models
    assert "model2" in models
    
    model1 = ModelRegistry.create_instance("model1")
    assert model1.name == "model1"
    assert ModelCategory.GENERAL in model1.categories
    
    model2 = ModelRegistry.create_instance("model2")
    assert ModelCategory.CODE in model2.categories
    
    print("Model info tests passed!")

if __name__ == "__main__":
    test_model_registry()
    test_model_info()
    print("All tests passed!")