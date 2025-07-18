"""Test configuration and shared fixtures."""

import os
import sys
import types
from typing import Any, Dict, Generator

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

# Create a minimal FastAPI app for unit tests to avoid importing the full server
app = FastAPI()

# Set test environment
os.environ["ENVIRONMENT"] = "test"

# Stub the heavy routers package so importing api modules doesn't pull full FastAPI routes
sys.modules.setdefault("api.routers", types.ModuleType("api.routers"))
sys.modules.setdefault("routers", sys.modules["api.routers"])


@pytest.fixture
def test_client() -> Generator[TestClient, None, None]:
    """Create a test client for the FastAPI app."""
    with TestClient(app) as client:
        yield client


@pytest.fixture
def mock_openai_response() -> Dict[str, Any]:
    """Mock response from OpenAI API."""
    return {
        "choices": [
            {
                "message": {"content": "Test response", "role": "assistant"},
                "finish_reason": "stop",
                "index": 0,
            }
        ],
        "model": "gpt-3.5-turbo",
        "object": "chat.completion",
        "usage": {"completion_tokens": 10, "prompt_tokens": 20, "total_tokens": 30},
    }


@pytest.fixture
def mock_pubchem_response() -> Dict[str, Any]:
    """Mock response from PubChem API."""
    return {
        "PC_Compounds": [
            {
                "id": {"id": {"cid": 962}},
                "props": [{"urn": {"label": "IUPAC Name"}, "value": {"sval": "water"}}],
            }
        ]
    }


@pytest.fixture
def test_env_vars(monkeypatch: pytest.MonkeyPatch) -> None:
    """Set up test environment variables."""
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setenv("ENVIRONMENT", "test")
    monkeypatch.setenv("REDIS_HOST", "localhost")
    monkeypatch.setenv("REDIS_PORT", "6379")
