"""
Test script for analyzing the propane visualization generation process using Claude 3.7 Sonnet
"""

import os
import json
import sys
from pathlib import Path
from typing import Dict, Any, TypedDict, Union
from dotenv import load_dotenv
from agent_management.llm_service import LLMResponse

# Add the api directory to the Python path
api_dir = Path(__file__).parent.parent
sys.path.append(str(api_dir))

import pytest
from agent_management.llm_service import LLMService, LLMRequest, LLMModelConfig, ProviderType
from agent_management.agents.geometry_agent import GeometryAgent
from routers.prompt.routes import GeometryRequest
from agent_management.utils.code_extraction import extract_code_block

# Load environment variables
load_dotenv()

# Constants
CLAUDE_MODEL = "claude-3-5-sonnet-20241022"
PROPANE_PROMPT = "Teach me about propane's 3D structure"

class MoleculeStructure(TypedDict):
    atoms: Dict[str, int]
    bonds: Dict[str, int]

class CodeQuality(TypedDict):
    has_error_handling: bool
    has_proper_materials: bool
    has_proper_geometries: bool
    has_proper_positioning: bool
    has_proper_bonds: bool

class Analysis(TypedDict):
    has_required_components: bool
    molecule_structure: MoleculeStructure
    code_quality: CodeQuality

EXPECTED_ELEMENTS: MoleculeStructure = {
    'atoms': {
        'carbon': 3,
        'hydrogen': 8
    },
    'bonds': {
        'carbon_carbon': 2,
        'carbon_hydrogen': 8
    }
}

def analyze_threejs_code(code: str) -> Analysis:
    """
    Analyze the generated Three.js code for correctness and completeness
    """
    analysis: Analysis = {
        'has_required_components': False,
        'molecule_structure': {
            'atoms': {'carbon': 0, 'hydrogen': 0},
            'bonds': {'carbon_carbon': 0, 'carbon_hydrogen': 0}
        },
        'code_quality': {
            'has_error_handling': False,
            'has_proper_materials': False,
            'has_proper_geometries': False,
            'has_proper_positioning': False,
            'has_proper_bonds': False
        }
    }
    
    # Check for required components
    required_components = [
        'scene', 'THREE.MeshPhongMaterial',
        'THREE.SphereGeometry', 'THREE.Group',
        'createBond'
    ]
    analysis['has_required_components'] = all(comp in code for comp in required_components)
    
    # Count carbon atoms (looking for carbonGeometry instantiations)
    analysis['molecule_structure']['atoms']['carbon'] = code.count('new THREE.Mesh(carbonGeometry')
    
    # Count hydrogen atoms
    analysis['molecule_structure']['atoms']['hydrogen'] = code.count('new THREE.Mesh(hydrogenGeometry')
    
    # Count C-C bonds (looking for createBond between carbons)
    analysis['molecule_structure']['bonds']['carbon_carbon'] = code.count('c1c2Bond') + code.count('c2c3Bond')
    
    # Count C-H bonds (total createBond calls minus C-C bonds)
    analysis['molecule_structure']['bonds']['carbon_hydrogen'] = code.count('createBond') - analysis['molecule_structure']['bonds']['carbon_carbon']
    
    # Analyze code quality
    code_quality: CodeQuality = {
        'has_error_handling': 'throw new Error' in code,
        'has_proper_materials': all(x in code for x in ['carbonMaterial', 'hydrogenMaterial', 'bondMaterial']),
        'has_proper_geometries': all(x in code for x in ['carbonGeometry', 'hydrogenGeometry']),
        'has_proper_positioning': 'position.set' in code,
        'has_proper_bonds': 'createBond' in code
    }
    analysis['code_quality'] = code_quality
    
    return analysis

def test_propane_visualization_generation():
    """Test the complete flow of generating a propane visualization"""
    
    # 1. Initialize the LLM service with Claude 3.7 Sonnet
    llm_service = LLMService(
        config=LLMModelConfig(
            provider=ProviderType.ANTHROPIC,
            model_name=CLAUDE_MODEL,
            api_key=os.getenv("ANTHROPIC_API_KEY")
        )
    )
    
    # 2. Create the geometry agent
    geometry_agent = GeometryAgent(llm_service=llm_service)
    
    # 3. Create the request
    request = GeometryRequest(
        prompt=PROPANE_PROMPT,
        model=CLAUDE_MODEL
    )
    
    try:
        # 4. Generate the geometry
        result = geometry_agent.get_geometry_snippet(request.prompt)
        
        # Handle both string and LLMResponse types
        if isinstance(result, LLMResponse):
            code = result.content
        else:
            code = result
            
        # Extract code if it's in a code block
        code = extract_code_block(code, "javascript")
        
        # 5. Analyze the generated code
        analysis = analyze_threejs_code(code)
        
        print("\nPropane Visualization Analysis:")
        print("Generated Code:")
        print(code)
        print("\nAnalysis Results:")
        print(json.dumps(analysis, indent=2))
        
        return True
        
    except Exception as e:
        print(f"Error during propane visualization test: {str(e)}")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    import traceback
    test_propane_visualization_generation() 