"""
Test cases for the prompt router endpoints
"""
import pytest
from fastapi.testclient import TestClient
from api.main import app
from unittest.mock import patch, MagicMock
from agent_management.agents.geometry_agent import GeometryAgent
from agent_management.llm_service import LLMService

client = TestClient(app)

def test_submit_prompt():
    """Test the basic prompt submission endpoint"""
    response = client.post(
        "/prompt/",
        json={"prompt": "test prompt"}
    )
    assert response.status_code == 200
    assert response.json() == {
        "prompt": "test prompt",
        "working": True
    }

@patch('agent_management.agents.geometry_agent.GeometryAgent.get_geometry_snippet')
def test_generate_geometry_success(mock_get_geometry):
    """Test successful geometry generation"""
    # Mock the geometry agent response
    mock_get_geometry.return_value = "// Test geometry code"
    
    response = client.post(
        "/prompt/generate-geometry/",
        json={"prompt": "create a simple sphere"}
    )
    
    assert response.status_code == 200
    assert "result" in response.json()
    assert isinstance(response.json()["result"], str)
    assert response.json()["result"] == "// Test geometry code"
    
    # Verify the mock was called correctly
    mock_get_geometry.assert_called_once_with("create a simple sphere")

@patch('agent_management.agents.geometry_agent.GeometryAgent.get_geometry_snippet')
def test_generate_geometry_empty_prompt(mock_get_geometry):
    """Test geometry generation with empty prompt"""
    mock_get_geometry.return_value = "// Empty geometry code"
    
    response = client.post(
        "/prompt/generate-geometry/",
        json={"prompt": ""}
    )
    assert response.status_code == 200
    assert response.json()["result"] == "// Empty geometry code"

def test_generate_geometry_invalid_request():
    """Test geometry generation with invalid request format"""
    response = client.post(
        "/prompt/generate-geometry/",
        json={"invalid_field": "test"}
    )
    assert response.status_code == 422  # Validation error

@patch('agent_management.agents.geometry_agent.GeometryAgent.get_geometry_snippet')
def test_generate_geometry_error_handling(mock_get_geometry):
    """Test error handling in geometry generation"""
    # Mock the geometry agent to raise an exception
    mock_get_geometry.side_effect = Exception("Test error")
    
    response = client.post(
        "/prompt/generate-geometry/",
        json={"prompt": "create a sphere"}
    )
    
    assert response.status_code == 500
    assert "detail" in response.json()
    assert response.json()["detail"] == "Test error" 