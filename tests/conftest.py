"""Test configuration and fixtures."""

import os
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

# Add the api directory to the Python path
api_dir = Path(__file__).parent.parent / "api"
sys.path.insert(0, str(api_dir))

# Set environment variables for testing
os.environ["ENVIRONMENT"] = "development"
os.environ["PYMOL_QUIET"] = "1"
os.environ["PYMOL_HEADLESS"] = "1"

# Mock PyMOL for testing if not available
try:
    import pymol
except ImportError:
    import types

    # Create a mock PyMOL module
    pymol = types.ModuleType("pymol")

    def _noop(*args, **kwargs):
        return None

    class _CmdStub:
        def reinitialize(self, *args, **kwargs):
            pass

        def do(self, *args, **kwargs):
            pass

        def png(self, *args, **kwargs):
            pass

        def save(self, *args, **kwargs):
            pass

        def get_view(self, *args, **kwargs):
            return [1.0] * 18

        def centerofmass(self, *args, **kwargs):
            return [0.0, 0.0, 0.0]

        def get_extent(self, *args, **kwargs):
            return [[-10.0, -10.0, -10.0], [10.0, 10.0, 10.0]]

        def set(self, *args, **kwargs):
            pass

        def ray(self, *args, **kwargs):
            pass

    pymol.finish_launching = _noop
    pymol.cmd = _CmdStub()
    sys.modules["pymol"] = pymol


@pytest.fixture
def client():
    """Create a test client for the FastAPI app."""
    from main import app

    return TestClient(app)


@pytest.fixture
def sample_yaml_spec():
    """Sample YAML specification for testing."""
    return """
meta:
  title: "Water Decomposition"
  version: "1.0.0"

canvas:
  w: 800
  h: 400

cells:
  - id: "diagram_1"
    type: "DIAGRAM"
    bbox:
      x: 0
      y: 0
      w: 800
      h: 400
    nodes:
      - id: "h2o"
        label: "H₂O"
      - id: "h2"
        label: "H₂"
      - id: "o2"
        label: "O₂"
    edges:
      - src: "h2o"
        dst: "h2"
        label: "decomposition"
      - src: "h2o"
        dst: "o2"
        label: "decomposition"
"""


@pytest.fixture
def sample_sdf_data():
    """Sample SDF data for testing."""
    return """
  Mrv2014 01010100002D

  3  2  0  0  0  0            999 V2000
   -0.4125    0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""
