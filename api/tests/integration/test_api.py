"""Integration tests for API endpoints."""

from typing import Generator

import pytest
from fastapi.testclient import TestClient


@pytest.mark.integration
class TestPromptEndpoints:
    """Test suite for prompt-related endpoints."""

    def test_submit_prompt(self, test_client: TestClient, test_env_vars: None) -> None:
        """Test submitting a prompt."""
        response = test_client.post(
            "/prompt/", json={"prompt": "Show water molecule", "model": "gpt-3.5-turbo"}
        )
        assert response.status_code == 200
        data = response.json()
        assert "job_id" in data

    def test_get_job_status(self, test_client: TestClient, test_env_vars: None) -> None:
        """Test getting job status."""
        # First create a job
        create_response = test_client.post(
            "/prompt/", json={"prompt": "Show water molecule", "model": "gpt-3.5-turbo"}
        )
        job_id = create_response.json()["job_id"]

        # Then check its status
        status_response = test_client.get(f"/prompt/process/{job_id}")
        assert status_response.status_code == 200
        data = status_response.json()
        assert "status" in data
        assert data["status"] in ["pending", "processing", "completed", "error"]


@pytest.mark.integration
class TestPubChemEndpoints:
    """Test suite for PubChem-related endpoints."""

    def test_fetch_molecule_data(
        self, test_client: TestClient, test_env_vars: None
    ) -> None:
        """Test fetching molecule data."""
        response = test_client.post(
            "/prompt/fetch-molecule-data/", json={"query": "water"}
        )
        assert response.status_code == 200
        data = response.json()
        assert "molecule_data" in data


@pytest.mark.integration
class TestHealthCheck:
    """Test suite for health check endpoints."""

    def test_health_check(self, test_client: TestClient) -> None:
        """Test health check endpoint."""
        response = test_client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
