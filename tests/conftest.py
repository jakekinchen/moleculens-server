"""Pytest configuration and fixtures."""

import os
from pathlib import Path

import pytest


@pytest.fixture
def golden_dir() -> Path:
    """Path to golden test data directory."""
    return Path(__file__).parent / "golden"


@pytest.fixture
def water_sdf(golden_dir: Path) -> str:
    """Load water.sdf content."""
    return (golden_dir / "water.sdf").read_text()


@pytest.fixture
def artifacts_dir() -> Path:
    """Path to artifacts output directory."""
    path = Path(__file__).parent.parent / "artifacts" / "golden_run"
    path.mkdir(parents=True, exist_ok=True)
    return path


@pytest.fixture
def api_url() -> str:
    """API URL for integration tests."""
    return os.environ.get("API_URL", "http://localhost:8001")
