"""
Test cases for the prompt router endpoints
"""

import pytest
from fastapi.testclient import TestClient
from api.main import app
import os

client = TestClient(app)


def test_submit_prompt():
    """Test the basic prompt submission endpoint"""
    response = client.post("/prompt/", json={"prompt": "test prompt"})
    assert response.status_code == 200
    assert response.json() == {"prompt": "test prompt", "working": True}


def test_generate_geometry_success():
    """Test successful geometry generation"""
    response = client.post(
        "/prompt/generate-geometry/", json={"prompt": "create a simple sphere"}
    )

    assert response.status_code == 200
    assert "result" in response.json()
    assert isinstance(response.json()["result"], str)
    assert "THREE" in response.json()["result"]  # Should contain Three.js code


def test_generate_geometry_empty_prompt():
    """Test geometry generation with empty prompt"""
    response = client.post("/prompt/generate-geometry/", json={"prompt": ""})
    assert response.status_code == 200
    assert "result" in response.json()
    assert isinstance(response.json()["result"], str)


def test_generate_geometry_invalid_request():
    """Test geometry generation with invalid request format"""
    response = client.post("/prompt/generate-geometry/", json={"invalid_field": "test"})
    assert response.status_code == 422  # Validation error


def test_generate_geometry_no_api_key():
    """Test error handling when API key is not set"""
    # Temporarily unset the API key
    api_key = os.environ.get("OPENAI_API_KEY")
    if api_key:
        del os.environ["OPENAI_API_KEY"]

    try:
        response = client.post(
            "/prompt/generate-geometry/", json={"prompt": "create a sphere"}
        )

        assert response.status_code == 500
        assert (
            "OPENAI_API_KEY environment variable is not set"
            in response.json()["detail"]
        )
    finally:
        # Restore the API key
        if api_key:
            os.environ["OPENAI_API_KEY"] = api_key


def test_fetch_molecule_layout(monkeypatch):
    """Test the multi-molecule layout endpoint."""

    dummy_result = [
        {
            "query": "water",
            "box": {"x": 0, "y": 0, "width": 10, "height": 10},
            "atoms": [],
            "bonds": [],
            "name": "water",
            "cid": 1,
            "formula": "H2O",
        }
    ]

    def fake_layout(self, queries):
        return dummy_result

    monkeypatch.setattr(
        "agent_management.agents.pubchem_agent.PubChemAgent.get_molecules_2d_layout",
        fake_layout,
    )

    response = client.post(
        "/prompt/fetch-molecule-layout/",
        json={
            "molecules": [
                {"query": "water", "box": {"x": 0, "y": 0, "width": 10, "height": 10}}
            ]
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["molecules"] == dummy_result
