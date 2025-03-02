from fastapi import Depends, Query, HTTPException, Request
from agent_management.llm_service import LLMService
from agent_management.model_config import get_llm_service, get_default_model, get_models_by_category, ModelRegistry, ModelCategory
from pydantic import BaseModel
import os
from typing import Optional


async def use_llm(
    request: Request,
    model_name: Optional[str] = Query(None, alias="model")
):
    """
    Dependency that initializes LLMService based on the model query parameter.
    Uses the ModelRegistry to validate and create the appropriate LLM service.
    
    Features:
    - Can specify model directly via query parameter
    - Can specify preferred_model_category in request body
    - Falls back to default model if nothing specified
    
    Args:
        request: The FastAPI request object to access the body
        model_name: The name of the model to use, or None to use the default
        
    Returns:
        An initialized LLMService instance
        
    Raises:
        HTTPException: If the model name is invalid
    """
    # First priority: Use model_name query parameter if provided
    if model_name is not None:
        try:
            # This will raise ValueError if model doesn't exist
            ModelRegistry.get_model(model_name)
            return get_llm_service(model_name)
        except ValueError:
            # Get a list of available models for the error message
            available_models = ModelRegistry.list_models()
            raise HTTPException(
                status_code=400, 
                detail=f"Invalid model name '{model_name}'. Allowed values: {', '.join(available_models)}"
            )
    
    # Second priority: Check if request body contains preferred_model_category
    try:
        # Only try to parse body if content type is JSON
        content_type = request.headers.get("content-type", "")
        if "application/json" in content_type:
            body = await request.json()
            
            # Check if preferred_model_category is in body
            if "preferred_model_category" in body:
                preferred_category = body["preferred_model_category"]
                
                # Try to convert the string to ModelCategory enum
                try:
                    category = ModelCategory(preferred_category)
                    # Get models for this category
                    category_models = get_models_by_category(category)
                    
                    if category_models:
                        # Use the first model in this category
                        return get_llm_service(category_models[0])
                except (ValueError, KeyError):
                    # Invalid category or no models in category, fall through to default
                    pass
    except Exception:
        # If any error happens parsing the body, just fall through to default
        pass
    
    # Third priority: Use the default model
    default_model = get_default_model()
    return get_llm_service(default_model)
