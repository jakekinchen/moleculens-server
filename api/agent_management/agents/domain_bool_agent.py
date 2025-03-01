"""
Domain Validator Agent - Validates whether a prompt is scientific in nature.
"""

import os
from agent_management.models import BooleanResponse
from agent_management.llm_service import (
    LLMService,
    StructuredLLMRequest,
    LLMModelConfig,
    ProviderType
)

class DomainValidator:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service
    
    def is_scientific(self, prompt: str) -> BooleanResponse:
        """
        Determine if a prompt is scientific in nature.
        
        Args:
            prompt: The user prompt to validate
            
        Returns:
            BooleanResponse with is_true=True if prompt is scientific, False otherwise
        """
        # Create a request for scientific content validation
        request = StructuredLLMRequest(
            user_prompt=f"""Is the following query scientific or not? '{prompt}'
            
Respond with a JSON object that has this exact structure:
{{
  "is_true": true,
  "confidence": 0.9,
  "reasoning": "explanation here"
}}

Where:
- is_true: boolean - true if the prompt is scientific in nature, false otherwise
- confidence: number between 0 and 1
- reasoning: string explaining your reasoning""",
            system_prompt="""You are a scientific content validation AI. Determine if user prompts are scientific in nature.

Scientific prompts typically involve:
- Natural phenomena (physics, chemistry, biology, astronomy, etc.)
- Mathematical concepts
- Engineering or technical scenarios
- Scientific experiments or observations
- Data visualization of scientific information

Non-scientific prompts typically involve:
- Fictional characters or scenarios
- Creative art without scientific context
- General everyday objects without scientific framing
- Abstract concepts without scientific basis
- Violent, harmful, or inappropriate content

Always respond with JSON exactly matching the requested format.""",
            llm_config=self.llm_service.config,
            response_model=BooleanResponse,
        )
        
        # Get the structured response
        return self.llm_service.generate_structured(request)
