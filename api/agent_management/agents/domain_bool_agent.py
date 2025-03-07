"""
Domain Validator Agent - Validates whether a prompt is scientific in nature.
"""

import os
from agent_management.models import BooleanResponse, MolecularStructure
from agent_management.llm_service import (
    LLMService,
    StructuredLLMRequest,
    LLMModelConfig,
    ProviderType
)

class DomainValidator:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service

    def is_macromolecule(self, prompt: str) -> BooleanResponse:
        """
        Determine if a prompt is about a macromolecule.
        """
        request = StructuredLLMRequest[BooleanResponse](
            user_prompt=f"""Does the following prompt contain something that can be built as a macromolecule? We are not talking about proteins in general, for instance, we would not consider the 20 standard amino acids as macromolecules, but we would consider DNA or RNA as macromolecules. '{prompt}'
            """,
            response_model=BooleanResponse
        )
        return self.llm_service.generate_structured(request)

    def molecular_structure(self, prompt: str) -> MolecularStructure:
        """
        Determine if a prompt can be built as molecular structure.
        """
        request = StructuredLLMRequest[MolecularStructure](
            user_prompt=f"""Only respond with a molecular structure in SMILES format. For instance: user prompt: 'Draw a molecule of aspirin' -> response: 'CC(=O)OC1=CC=CC=C1C(=O)O'
            The user prompt is: '{prompt}'
            """,
            response_model=MolecularStructure
        )
        return self.llm_service.generate_structured(request)
    
    def is_molecular(self, prompt: str) -> BooleanResponse:
        """
        Determine if a prompt can be built as molecular structure.
        
        Args:
            prompt: The user prompt to validate
            
        Returns:
            BooleanResponse with is_true=True if prompt is molecular, False otherwise
        """
        # Create a request for scientific content validation
        request = StructuredLLMRequest[BooleanResponse](
            user_prompt=f"""Does the following prompt contain something that can be built as a molecular structure? '{prompt}'
            
Respond with a JSON object that has this exact structure:
{{
  "is_true": true,
}}

Where:
- is_true: boolean - true if the prompt is molecular, false otherwise""",
            system_prompt="""You are a molecular structure validation AI. Determine if user prompts are molecular in nature or at least related to molecular structures.

Molecular prompts typically involve:
- Chemical compounds
- Biological molecules
- Physical structures
- Molecular structures

Non-molecular prompts typically involve:
- Fictional characters or scenarios
- Creative art without molecular context
- General everyday objects without molecular framing
- Abstract concepts without molecular basis
- Violent, harmful, or inappropriate content

Examples of molecular prompts:
- "I want to learn about the structure of water"
- "Teach me about ATP"
- "Display a protein structure"
- "Draw a molecule of sugar"
- "Draw a molecule of salt"
- "Show the structure of a DNA molecule"
- "Display a protein structure"
- "Draw a molecule of cocaine"
- "Draw an oil molecule"

Examples of non-molecular prompts:
- "Draw a car"
- "Show a house"
- "Display a person"


Always respond with JSON exactly matching the requested format.""",
            response_model=BooleanResponse,
        )
        
        # Get the structured response
        return self.llm_service.generate_structured(request)
