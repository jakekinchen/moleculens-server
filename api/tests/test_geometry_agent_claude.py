"""
Test suite for the GeometryAgent with Claude 3.5 Sonnet.
Tests various molecular structures and geometry generation scenarios.
"""

import os
import pytest
from dotenv import load_dotenv
from agent_management.agent_factory import AgentFactory
from agent_management.agents.geometry_agent import GeometryAgent
from agent_management.llm_service import LLMService
from agent_management.utils.code_extraction import extract_code_block

def validate_threejs_code(code: str) -> list:
    """Validate Three.js code for common issues"""
    issues = []
    
    # Check for required components
    if 'scene.add' not in code:
        issues.append("Missing scene.add() call")
    
    # Check for proper variable declarations
    if not all(line.strip().startswith(('const ', 'let ', 'var ', 'function ', '//', '')) 
               for line in code.split('\n') if line.strip()):
        issues.append("Some variables may not be properly declared")
    
    # Check for proper semicolon usage in complete statements
    lines = code.split('\n')
    statement_buffer = []
    
    for i, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith('//'):
            continue
            
        statement_buffer.append(line)
        
        # Check if this is a complete statement
        if line.endswith(';') or line.endswith('{') or line.endswith('}'):
            statement_buffer = []
            continue
            
        # Check if this is the end of a multi-line statement
        if statement_buffer and line.endswith(')'):
            full_statement = ' '.join(statement_buffer)
            if not full_statement.endswith(';'):
                issues.append(f"Missing semicolon at end of multi-line statement near line {i}")
            statement_buffer = []
    
    return issues

@pytest.fixture(scope="module")
def geometry_agent():
    """Create a geometry agent with Claude 3.5 Sonnet for testing"""
    load_dotenv()
    return AgentFactory.create_geometry_agent("claude-3-sonnet-20240229")

def test_water_molecule(geometry_agent):
    """Test generating a water molecule with proper bond angles"""
    prompt = "Create a water molecule (H2O) with appropriate bond angles"
    code = geometry_agent.get_geometry_snippet(prompt)
    
    # Validate the code
    issues = validate_threejs_code(code)
    assert len(issues) == 0, f"Found issues in generated code: {issues}"
    
    # Check for required components
    assert "THREE.MeshPhongMaterial" in code, "Missing material definitions"
    assert "THREE.SphereGeometry" in code, "Missing sphere geometry"
    assert "THREE.Group" in code, "Missing group creation"
    assert "scene.add" in code, "Missing scene.add() call"
    
    # Check for proper bond angle
    assert "104.5" in code or "104.45" in code, "Incorrect or missing H2O bond angle"

def test_methane_molecule(geometry_agent):
    """Test generating a methane molecule with tetrahedral structure"""
    prompt = "Create a methane molecule (CH4) with tetrahedral structure"
    code = geometry_agent.get_geometry_snippet(prompt)
    
    # Validate the code
    issues = validate_threejs_code(code)
    assert len(issues) == 0, f"Found issues in generated code: {issues}"
    
    # Check for required components
    assert "THREE.MeshPhongMaterial" in code, "Missing material definitions"
    assert "THREE.SphereGeometry" in code, "Missing sphere geometry"
    assert "THREE.Group" in code, "Missing group creation"
    
    # Check for proper structure hints
    assert "109.5" in code or "109.47" in code, "Missing tetrahedral angle"

def test_error_handling(geometry_agent):
    """Test error handling with invalid prompts"""
    prompt = ""  # Empty prompt
    code = geometry_agent.get_geometry_snippet(prompt)
    
    # Should still generate valid Three.js code
    issues = validate_threejs_code(code)
    assert len(issues) == 0, f"Found issues in generated code: {issues}"

def test_complex_molecule(geometry_agent):
    """Test generating a more complex molecule (ethanol)"""
    prompt = "Create an ethanol molecule (C2H5OH) with proper bond angles"
    code = geometry_agent.get_geometry_snippet(prompt)
    
    # Validate the code
    issues = validate_threejs_code(code)
    assert len(issues) == 0, f"Found issues in generated code: {issues}"
    
    # Check for required components
    assert "THREE.MeshPhongMaterial" in code, "Missing material definitions"
    assert "THREE.SphereGeometry" in code, "Missing sphere geometry"
    assert "createBond" in code, "Missing bond creation function"
    
    # Check for proper structure
    assert code.count("const") >= 8, "Not enough atoms defined for ethanol"

def test_code_style(geometry_agent):
    """Test code style and formatting"""
    prompt = "Create a simple sphere"
    code = geometry_agent.get_geometry_snippet(prompt)
    
    # Check for consistent spacing and indentation
    lines = code.split('\n')
    indented_lines = [line for line in lines if line.strip() and line.startswith('    ')]
    assert len(indented_lines) > 0, "No properly indented lines found"
    
    # Check for comments
    assert code.count('//') >= 2, "Insufficient code documentation"

if __name__ == "__main__":
    pytest.main([__file__]) 