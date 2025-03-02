"""
Test module for the Model Registry implementation
"""

from api.agent_management.models import ModelRegistry, BaseModel
from api.agent_management.model_config import (
    ModelInfo, ModelCategory, register_models, 
    get_default_model, get_models_by_category
)

def test_model_registry():
    """Test the basic functionality of the model registry"""
    
    # Clear registry for testing
    ModelRegistry._registry = {}
    
    # Define a simple test model
    class TestModel(BaseModel):
        name: str
        value: int
    
    # Register model
    ModelRegistry.register(
        "test-model", 
        TestModel,
        lambda: TestModel(name="default", value=42)
    )
    
    # Get model class
    model_cls = ModelRegistry.get_model("test-model")
    assert model_cls == TestModel
    
    # Create instance
    instance = ModelRegistry.create_instance("test-model")
    assert isinstance(instance, TestModel)
    assert instance.name == "default"
    assert instance.value == 42
    
    # Create with custom parameters
    custom = ModelRegistry.create_instance("test-model", name="custom", value=99)
    assert custom.name == "custom"
    assert custom.value == 99
    
    print("ModelRegistry basic tests passed!")

def test_model_config():
    """Test the model configuration module"""
    
    # Initialize the registry
    register_models()
    
    # Verify models were registered
    models = ModelRegistry.list_models()
    assert len(models) > 0
    assert "o3-mini" in models
    assert "claude-3-7-sonnet-latest" in models
    
    # Get default model
    default = get_default_model()
    assert default in models
    
    # Get models by category
    code_models = get_models_by_category(ModelCategory.CODE)
    assert len(code_models) > 0
    # At least claude should be in code models
    assert "claude-3-7-sonnet-latest" in code_models or "gpt-4.5-preview" in code_models
    
    # Verify model info
    for model_name in models:
        model_info = ModelRegistry.create_instance(model_name)
        assert model_info.name == model_name
        assert model_info.display_name
        assert model_info.provider
        assert model_info.context_length > 0
        assert len(model_info.categories) > 0
    
    print("Model config tests passed!")

if __name__ == "__main__":
    test_model_registry()
    test_model_config()
    print("All tests passed!")