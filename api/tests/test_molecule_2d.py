import json
from types import SimpleNamespace
import requests
import pytest
from typing import Dict, Any

from agent_management.agents.pubchem_agent import PubChemAgent
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType


def make_agent():
    config = LLMModelConfig(provider=ProviderType.OPENAI, model_name="gpt-3.5-turbo")
    llm_service = LLMService(config)
    return PubChemAgent(llm_service)


@pytest.fixture()
def sample_sdf():
    # Simple water molecule SDF
    return """
  -ISIS-  05202115492D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9428    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2357    0.8981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
"""


def test_get_molecule_2d_info_mock(monkeypatch, sample_sdf):
    """Test the 2D molecule info endpoint with mocked data."""
    agent = make_agent()

    class FakeCompound(SimpleNamespace):
        pass

    fake_compound = FakeCompound(
        cid=962,  # Real CID for water
        iupac_name="Water",
        molecular_formula="H2O",
    )

    def fake_interpret(query: str):
        return "water"

    def fake_search(query: str):
        return [fake_compound]

    class FakeResponse:
        status_code = 200
        text = sample_sdf

    def fake_get(url, timeout=30):
        return FakeResponse()

    monkeypatch.setattr(agent, "interpret_user_query", fake_interpret)
    monkeypatch.setattr(agent, "_search_with_fallbacks", fake_search)
    monkeypatch.setattr(requests, "get", fake_get)

    result = agent.get_molecule_2d_info("water")
    assert result["cid"] == 962
    assert len(result["atoms"]) == 3  # 1 O + 2 H
    assert len(result["bonds"]) == 2  # 2 O-H bonds
    assert result["formula"] == "H2O"
    
    # Validate atom types
    elements = [atom["element"] for atom in result["atoms"]]
    assert elements.count("O") == 1
    assert elements.count("H") == 2
    
    # Validate bond orders
    assert all(bond["order"] == 1.0 for bond in result["bonds"])  # All single bonds in water


def test_get_molecule_2d_info_real():
    """Test the 2D molecule info endpoint with real PubChem data for water."""
    agent = make_agent()
    result = agent.get_molecule_2d_info("water")
    
    # Basic structure checks
    assert isinstance(result, dict)
    assert all(key in result for key in ["atoms", "bonds", "name", "cid", "formula"])
    
    # Water molecule specific checks
    assert result["formula"] == "H2O"
    assert len(result["atoms"]) == 3  # 1 O + 2 H
    assert len(result["bonds"]) == 2  # 2 O-H bonds
    
    # Check atom structure
    for atom in result["atoms"]:
        assert all(key in atom for key in ["element", "x", "y"])
        assert atom["element"] in ["H", "O"]
        assert isinstance(atom["x"], (int, float))
        assert isinstance(atom["y"], (int, float))
    
    # Check bond structure
    for bond in result["bonds"]:
        assert all(key in bond for key in ["start", "end", "order"])
        assert isinstance(bond["start"], int)
        assert isinstance(bond["end"], int)
        assert bond["order"] == 1.0  # Single bonds in water


def test_get_molecule_2d_info_complex():
    """Test the 2D molecule info endpoint with a more complex molecule (benzene)."""
    agent = make_agent()
    result = agent.get_molecule_2d_info("benzene")
    
    # Basic structure checks
    assert isinstance(result, dict)
    assert all(key in result for key in ["atoms", "bonds", "name", "cid", "formula"])
    
    # Benzene specific checks
    assert result["formula"] == "C6H6"
    assert len(result["atoms"]) == 12  # 6 C + 6 H
    
    # Count carbon and hydrogen atoms
    elements = [atom["element"] for atom in result["atoms"]]
    assert elements.count("C") == 6
    assert elements.count("H") == 6
    
    # Check bonds
    carbon_bonds = []
    ch_bonds = []
    
    for bond in result["bonds"]:
        # Get the elements involved in this bond
        start_element = result["atoms"][bond["start"]]["element"]
        end_element = result["atoms"][bond["end"]]["element"]
        
        if start_element == "C" and end_element == "C":
            carbon_bonds.append(bond)
        elif (start_element == "C" and end_element == "H") or (start_element == "H" and end_element == "C"):
            ch_bonds.append(bond)
    
    # Verify we have 6 C-C bonds forming the ring and 6 C-H bonds
    assert len(carbon_bonds) == 6, "Should have 6 C-C bonds"
    assert len(ch_bonds) == 6, "Should have 6 C-H bonds"
    
    # All C-H bonds should be single bonds
    assert all(bond["order"] == 1.0 for bond in ch_bonds), "All C-H bonds should be single bonds"
    
    # Verify the ring structure by checking connectivity
    carbon_indices = [i for i, atom in enumerate(result["atoms"]) if atom["element"] == "C"]
    ring_bonds = set()
    for bond in carbon_bonds:
        if bond["start"] in carbon_indices and bond["end"] in carbon_indices:
            ring_bonds.add((min(bond["start"], bond["end"]), max(bond["start"], bond["end"])))
    
    assert len(ring_bonds) == 6, "Should have 6 unique bonds in the ring"


def test_get_molecule_2d_info_invalid():
    """Test the 2D molecule info endpoint with invalid input."""
    agent = make_agent()
    
    with pytest.raises(ValueError) as exc_info:
        agent.get_molecule_2d_info("not_a_real_molecule_name_123456789")
    assert "No compounds found" in str(exc_info.value)


def test_get_molecule_2d_info_large():
    """Test the 2D molecule info endpoint with a larger molecule (caffeine)."""
    agent = make_agent()
    result = agent.get_molecule_2d_info("caffeine")
    
    # Basic structure checks
    assert isinstance(result, dict)
    assert all(key in result for key in ["atoms", "bonds", "name", "cid", "formula"])
    
    # Caffeine specific checks
    assert result["formula"] == "C8H10N4O2"
    assert len(result["atoms"]) == 24  # 8 C + 10 H + 4 N + 2 O
    
    # Check coordinate system
    x_coords = [atom["x"] for atom in result["atoms"]]
    y_coords = [atom["y"] for atom in result["atoms"]]
    
    # Verify the molecule is properly laid out in 2D space
    assert max(x_coords) - min(x_coords) > 0  # Has width
    assert max(y_coords) - min(y_coords) > 0  # Has height


def validate_molecule_response(response: Dict[str, Any]) -> bool:
    """Helper function to validate the structure of a molecule response."""
    # Check required keys
    required_keys = ["atoms", "bonds", "name", "cid", "formula"]
    if not all(key in response for key in required_keys):
        return False
    
    # Validate atoms
    for atom in response["atoms"]:
        if not all(key in atom for key in ["element", "x", "y"]):
            return False
        if not isinstance(atom["element"], str):
            return False
        if not isinstance(atom["x"], (int, float)):
            return False
        if not isinstance(atom["y"], (int, float)):
            return False
    
    # Validate bonds
    for bond in response["bonds"]:
        if not all(key in bond for key in ["start", "end", "order"]):
            return False
        if not isinstance(bond["start"], int):
            return False
        if not isinstance(bond["end"], int):
            return False
        if not isinstance(bond["order"], (int, float)):
            return False
        
        # Validate bond indices
        if bond["start"] >= len(response["atoms"]) or bond["end"] >= len(response["atoms"]):
            return False
    
    return True


def test_molecule_response_validation():
    """Test the response validation for various molecules."""
    agent = make_agent()
    molecules = ["water", "methane", "ethanol", "benzene"]
    
    for molecule in molecules:
        result = agent.get_molecule_2d_info(molecule)
        assert validate_molecule_response(result), f"Invalid response structure for {molecule}"
