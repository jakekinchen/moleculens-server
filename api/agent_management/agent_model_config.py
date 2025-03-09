"""
Agent-Model Configuration Module.

This module defines which models should be used by different agents in the pipeline,
allowing for centralized control over model selection for different tasks.
"""

from enum import Enum
from typing import Dict, List, Optional
from pydantic import BaseModel

from agent_management.llm_service import ProviderType
from agent_management.model_config import ModelCategory

class AgentType(str, Enum):
    """Types of agents in the visualization pipeline"""
    DOMAIN_VALIDATOR = "domain_validator"
    SCRIPT = "script" 
    ORCHESTRATION = "orchestration"
    GEOMETRY = "geometry"
    ANIMATION = "animation"
    CAPTION = "caption"
    AGGREGATOR = "aggregator"
    PUBCHEM = "pubchem"

class AgentModelConfig(BaseModel):
    """Configuration for which model an agent should use"""
    agent_type: AgentType
    preferred_model: str
    fallback_models: List[str] = []
    required_categories: List[ModelCategory] = []
    provider_preference: Optional[ProviderType] = None
    description: str

# Default agent model configurations
DEFAULT_AGENT_MODELS = [
    AgentModelConfig(
        agent_type=AgentType.DOMAIN_VALIDATOR,
        preferred_model="o3-mini",  # Fast, cost-effective
        fallback_models=["gpt-4o", "claude-3-5-sonnet-latest"],
        required_categories=[ModelCategory.GENERAL],
        description="Validates if user prompts are scientific in nature"
    ),
    
    AgentModelConfig(
        agent_type=AgentType.SCRIPT,
        preferred_model="o3-mini",  # High quality required for script generation
        fallback_models=["claude-3-7-sonnet-latest", "gpt-4.5-preview"],
        required_categories=[ModelCategory.REASONING],
        description="Generates animation scripts with timecodes and captions"
    ),
    
    AgentModelConfig(
        agent_type=AgentType.ORCHESTRATION,
        preferred_model="gpt-4o",  # Planning requires strong reasoning
        fallback_models=["claude-3-7-sonnet-latest", "gpt-4.5-preview"],
        required_categories=[ModelCategory.REASONING],
        description="Plans which objects are needed for the animation"
    ),
    
    AgentModelConfig(
        agent_type=AgentType.GEOMETRY,
        preferred_model="claude-3-5-sonnet-latest",  # Code generation quality
        fallback_models=["gpt-4.5-preview", "qwen-2.5-coder-32b"],
        required_categories=[ModelCategory.CODE],
        description="Generates Three.js geometry code for objects"
    ),
    
    AgentModelConfig(
        agent_type=AgentType.ANIMATION,
        preferred_model="claude-3-7-sonnet-latest",  # Best for code generation
        fallback_models=["gpt-4.5-preview", "gpt-4o"],
        required_categories=[ModelCategory.CODE],
        description="Generates animation code that manipulates objects over time"
    ),
    
    AgentModelConfig(
        agent_type=AgentType.CAPTION,
        preferred_model="o3-mini",  # Simple task, cost-effective model
        fallback_models=["gpt-4o", "llama3-70b-8192"],
        required_categories=[ModelCategory.GENERAL],
        description="Generates captions for animation frames"
    ),
    
    AgentModelConfig(
        agent_type=AgentType.AGGREGATOR,
        preferred_model="o3-mini",  # Simple aggregation task
        fallback_models=["gpt-4o", "llama3-70b-8192"],
        required_categories=[ModelCategory.GENERAL],
        description="Combines outputs from various agents"
    ),
    
    AgentModelConfig(
        agent_type=AgentType.PUBCHEM,
        preferred_model="o3-mini",  # Changed to OpenAI model for reliability
        fallback_models=["gpt-4o", "claude-3-7-sonnet-latest", "gpt-4.5-preview"],
        required_categories=[ModelCategory.REASONING],
        description="Handles PubChem data retrieval and processing"
    )
]

# Create a lookup dictionary for faster access
AGENT_MODEL_MAP: Dict[AgentType, AgentModelConfig] = {
    config.agent_type: config for config in DEFAULT_AGENT_MODELS
}

def get_model_for_agent(agent_type: AgentType) -> str:
    """
    Get the preferred model name for a specific agent type.
    
    Args:
        agent_type: The type of agent
        
    Returns:
        The name of the preferred model for this agent
        
    Raises:
        ValueError: If the agent type is not configured
    """
    if agent_type not in AGENT_MODEL_MAP:
        raise ValueError(f"No model configuration for agent type: {agent_type}")
    
    return AGENT_MODEL_MAP[agent_type].preferred_model

def get_agent_config(agent_type: AgentType) -> AgentModelConfig:
    """
    Get the full model configuration for a specific agent type.
    
    Args:
        agent_type: The type of agent
        
    Returns:
        The full agent model configuration
        
    Raises:
        ValueError: If the agent type is not configured
    """
    if agent_type not in AGENT_MODEL_MAP:
        raise ValueError(f"No model configuration for agent type: {agent_type}")
    
    return AGENT_MODEL_MAP[agent_type]

def create_agent_llm_service(agent_type: AgentType, override_model: Optional[str] = None) -> 'LLMService':
    """
    Create an LLMService configured for a specific agent.
    
    Args:
        agent_type: The type of agent
        override_model: Optional model name to override the default
        
    Returns:
        An initialized LLMService instance configured for the agent
        
    Raises:
        ValueError: If the agent type is not configured or if model is not found
    """
    from agent_management.model_config import get_llm_service, ModelRegistry
    
    # Get the agent's model configuration
    config = get_agent_config(agent_type)
    
    # Use the override model if provided, otherwise use the preferred model
    model_name = override_model or config.preferred_model
    
    # Try to create an LLM service with the selected model
    try:
        return get_llm_service(model_name)
    except ValueError:
        # If the preferred model is not available, try fallback models
        if not override_model and config.fallback_models:
            for fallback_model in config.fallback_models:
                try:
                    return get_llm_service(fallback_model)
                except ValueError:
                    continue
        
        # If we still don't have a model, try to find one that matches the requirements
        if config.required_categories:
            from agent_management.model_config import get_models_by_category
            
            for category in config.required_categories:
                models = get_models_by_category(category)
                if models:
                    return get_llm_service(models[0])
        
        # If we get here, we couldn't find a suitable model
        raise ValueError(
            f"No suitable model found for agent {agent_type}. "
            f"Tried preferred model '{config.preferred_model}' and fallbacks."
        )