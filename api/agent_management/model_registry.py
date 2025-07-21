"""Model registry for managing model classes and factory functions."""

from typing import Any, Callable, Dict, List, Tuple, Type

from pydantic import BaseModel


class ModelRegistry:
    """Registry for managing model classes and factory functions."""

    _registry: Dict[str, Tuple[Type[BaseModel], Callable[..., Any]]] = {}

    @classmethod
    def register(
        cls,
        model_name: str,
        model_cls: Type[BaseModel],
        factory_func: Callable[..., Any],
    ) -> None:
        """Register a model class and its factory function.

        Args:
            model_name: Name to register the model under
            model_cls: The model class
            factory_func: Factory function to create model instances
        """
        cls._registry[model_name] = (model_cls, factory_func)

    @classmethod
    def get_model(cls, model_name: str) -> Type[BaseModel]:
        """Get the registered model class.

        Args:
            model_name: Name of the registered model

        Returns:
            The model class

        Raises:
            ValueError: If the model is not registered
        """
        if model_name not in cls._registry:
            raise ValueError(f"Model '{model_name}' not registered.")
        return cls._registry[model_name][0]

    @classmethod
    def create_instance(cls, model_name: str, *args: Any, **kwargs: Any) -> Any:
        """Create an instance of a registered model.

        Args:
            model_name: Name of the registered model
            *args: Arguments to pass to the factory function
            **kwargs: Keyword arguments to pass to the factory function

        Returns:
            An instance of the model

        Raises:
            ValueError: If the model is not registered
        """
        if model_name not in cls._registry:
            raise ValueError(f"Model '{model_name}' not registered.")
        _, factory_func = cls._registry[model_name]
        return factory_func(*args, **kwargs)

    @classmethod
    def list_models(cls) -> List[str]:
        """Get a list of all registered model names."""
        return list(cls._registry.keys())
