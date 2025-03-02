"""
Agent Factory Module - Creates agent instances with appropriate model configurations.
"""

from typing import Dict, Optional, Any

from agent_management.agent_model_config import AgentType, create_agent_llm_service
from agent_management.agents.domain_bool_agent import DomainValidator
from agent_management.agents.script_agent import ScriptAgent
from agent_management.agents.orchestration_agent import OrchestrationAgent
from agent_management.agents.geometry_agent import GeometryAgent
from agent_management.agents.animation_agent import AnimationAgent
from agent_management.agents.caption_agent import CaptionAgent
from agent_management.agents.aggregator_agent import AggregatorAgent

class AgentFactory:
    """Factory for creating agent instances with appropriate model configurations."""
    
    @staticmethod
    def create_domain_validator(override_model: Optional[str] = None) -> DomainValidator:
        """Create a domain validator agent"""
        llm_service = create_agent_llm_service(AgentType.DOMAIN_VALIDATOR, override_model)
        return DomainValidator(llm_service)
    
    @staticmethod
    def create_script_agent(override_model: Optional[str] = None) -> ScriptAgent:
        """Create a script generation agent"""
        llm_service = create_agent_llm_service(AgentType.SCRIPT, override_model)
        return ScriptAgent(llm_service)
    
    @staticmethod
    def create_orchestration_agent(override_model: Optional[str] = None) -> OrchestrationAgent:
        """Create an orchestration agent"""
        llm_service = create_agent_llm_service(AgentType.ORCHESTRATION, override_model)
        return OrchestrationAgent(llm_service)
    
    @staticmethod
    def create_geometry_agent(override_model: Optional[str] = None) -> GeometryAgent:
        """Create a geometry agent"""
        llm_service = create_agent_llm_service(AgentType.GEOMETRY, override_model)
        return GeometryAgent(llm_service)
    
    @staticmethod
    def create_animation_agent(override_model: Optional[str] = None) -> AnimationAgent:
        """Create an animation agent"""
        llm_service = create_agent_llm_service(AgentType.ANIMATION, override_model)
        return AnimationAgent(llm_service)
    
    @staticmethod
    def create_caption_agent(override_model: Optional[str] = None) -> CaptionAgent:
        """Create a caption agent"""
        llm_service = create_agent_llm_service(AgentType.CAPTION, override_model)
        return CaptionAgent(llm_service)
    
    @staticmethod
    def create_aggregator_agent(override_model: Optional[str] = None) -> AggregatorAgent:
        """Create an aggregator agent"""
        llm_service = create_agent_llm_service(AgentType.AGGREGATOR, override_model)
        return AggregatorAgent(llm_service)
    
    @staticmethod
    def create_all_agents(global_override_model: Optional[str] = None) -> Dict[AgentType, Any]:
        """
        Create all agents with appropriate model configurations.
        
        Args:
            global_override_model: Optional model name to override all agents
            
        Returns:
            Dictionary of agent instances by type
        """
        agents = {
            AgentType.DOMAIN_VALIDATOR: AgentFactory.create_domain_validator(global_override_model),
            AgentType.SCRIPT: AgentFactory.create_script_agent(global_override_model),
            AgentType.ORCHESTRATION: AgentFactory.create_orchestration_agent(global_override_model),
            AgentType.GEOMETRY: AgentFactory.create_geometry_agent(global_override_model),
            AgentType.ANIMATION: AgentFactory.create_animation_agent(global_override_model),
            AgentType.CAPTION: AgentFactory.create_caption_agent(global_override_model),
            AgentType.AGGREGATOR: AgentFactory.create_aggregator_agent(global_override_model)
        }
        return agents