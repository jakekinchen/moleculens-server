import pytest
from fastapi.testclient import TestClient
from unittest.mock import patch, MagicMock
import base64
import re
import json # Import json for stringifying the mock response
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import pubchempy as pcp  # type: ignore
import os # For API key check and file path

from api.main import app
from api.routers.prompt.routes import DiagramPromptRequest, DiagramResponse, DiagramPlan, MoleculePlacement
from agent_management.agents.pubchem_agent import PubChemAgent

client = TestClient(app)

# Sample test data for a real water decomposition reaction
WATER_DECOMP_PROMPT = "show a water decomposition reaction for a group of middle schoolers"

# Define expected canvas sizes for the mock, matching a typical request
MOCK_CANVAS_WIDTH = 960
MOCK_CANVAS_HEIGHT = 640

WATER_DECOMP_LLM_RESPONSE_DICT = {
    "plan": "Create a diagram showing water (H2O) decomposing into hydrogen (H2) and oxygen (O2) molecules, with labels and an arrow indicating the reaction direction.",
    "molecule_list": [
        {
            "molecule": "H2O", "x": 320.0, "y": 360.0, "width": 200.0, "height": 200.0, 
            "label": "2H₂O", "label_position": "below"
        },
        {
            "molecule": "H2", "x": 640.0, "y": 260.0, "width": 200.0, "height": 200.0, 
            "label": "2H₂", "label_position": "above"
        },
        {
            "molecule": "O2", "x": 640.0, "y": 460.0, "width": 200.0, "height": 200.0, 
            "label": "O₂", "label_position": "below"
        }
    ],
    "arrows": [
        {
            "start": [420.0, 360.0], "end": [540.0, 360.0], 
            "style": "straight", "text": "electricity"
        }
    ],
    # Although the LLM won't return these, our endpoint adds them, so the *final* DiagramPlan object
    # that LLMService.generate_structured is mocked to return should represent this augmented state.
    "canvas_width": MOCK_CANVAS_WIDTH, 
    "canvas_height": MOCK_CANVAS_HEIGHT
}

# This string is what the LLM provider would *return* before pydantic validation by LLMService
WATER_DECOMP_LLM_RAW_JSON_STR = json.dumps(WATER_DECOMP_LLM_RESPONSE_DICT)
# This is the object LLMService.generate_structured should return after successful parsing & validation
EXPECTED_DIAGRAM_PLAN_OBJECT = DiagramPlan(**WATER_DECOMP_LLM_RESPONSE_DICT)

def verify_svg_content(svg_str: str) -> bool:
    """Verify SVG string contains expected elements."""
    stripped_svg = svg_str.strip()
    if not stripped_svg.startswith('<svg'): # Allow only <svg> start after stripping
        print(f"DEBUG: SVG does not start with <svg>. Starts with: {stripped_svg[:100]}")
        return False
    # Check for common SVG elements that should be present in a molecule drawing
    # A simple line drawing might not have <path>, but should have <line> or <circle> and <text>
    # Let's make this check more flexible for now or specific to what diagram_renderer produces.
    # For now, let's assume text for labels and line/circle for drawings are minimal.
    required_elements = ['<text'] # Labels should always be there
    # Check for drawing elements like <line>, <circle>, <path>
    has_drawing_element = any(elem_type in stripped_svg for elem_type in ['<line', '<circle', '<path'])
    
    if not has_drawing_element:
        print(f"DEBUG: SVG does not contain common drawing elements (line, circle, path).")
        return False
        
    return all(elem in stripped_svg for elem in required_elements)

def verify_molecule_structure(molecule_data):
    """Verify molecule data has valid structure."""
    assert "atoms" in molecule_data
    assert "bonds" in molecule_data
    assert len(molecule_data["atoms"]) > 0
    assert len(molecule_data["bonds"]) > 0
    
    # Verify atom data
    for atom in molecule_data["atoms"]:
        assert "element" in atom
        assert "x" in atom
        assert "y" in atom
        assert isinstance(atom["x"], (int, float))
        assert isinstance(atom["y"], (int, float))
    
    # Verify bond data
    for bond in molecule_data["bonds"]:
        assert "start" in bond
        assert "end" in bond
        assert "order" in bond
        assert bond["start"] < len(molecule_data["atoms"])
        assert bond["end"] < len(molecule_data["atoms"])

@pytest.fixture
def mock_openai_provider_structured():
    # Path to the actual provider class used by LLMService
    # Assuming OpenAIProvider is at agent_management.providers.openai_provider.OpenAIProvider
    with patch("agent_management.providers.openai_provider.OpenAIProvider.generate_structured") as mock_method:
        mock_method.return_value = EXPECTED_DIAGRAM_PLAN_OBJECT 
        yield mock_method

@pytest.fixture
def mock_openai_provider_error():
    with patch("agent_management.providers.openai_provider.OpenAIProvider.generate_structured") as mock_method:
        mock_method.side_effect = Exception("LLM provider simulated error")
        yield mock_method

def test_real_pubchem_fetch():
    """Test actual PubChem data fetching for water."""
    agent = PubChemAgent(None)  # LLM service not needed for this test
    
    # Test water molecule
    water_data = agent.get_molecule_2d_info("water")
    verify_molecule_structure(water_data)
    assert water_data["formula"] == "H2O"
    
    # Test hydrogen molecule
    h2_data = agent.get_molecule_2d_info("H2")
    verify_molecule_structure(h2_data)
    assert h2_data["formula"] == "H2"
    
    # Test oxygen molecule
    o2_data = agent.get_molecule_2d_info("O2")
    verify_molecule_structure(o2_data)
    assert o2_data["formula"] == "O2"

def test_real_molecule_layout():
    """Test actual molecule layout generation with RDKit."""
    agent = PubChemAgent(None)
    
    # Request layout for water decomposition molecules
    layout_requests = [
        {"query": "H2O", "box": {"x": 320, "y": 360, "width": 200, "height": 200}},
        {"query": "H2", "box": {"x": 640, "y": 260, "width": 200, "height": 200}},
        {"query": "O2", "box": {"x": 640, "y": 460, "width": 200, "height": 200}}
    ]
    
    layout_data = agent.get_molecules_2d_layout(layout_requests)
    assert len(layout_data) == 3
    
    # Verify each molecule's data
    for mol_data in layout_data:
        verify_molecule_structure(mol_data)
        assert "box" in mol_data
        assert all(k in mol_data["box"] for k in ["x", "y", "width", "height"])

def test_generate_molecule_diagram_endpoint(mock_openai_provider_structured):
    """Test the diagram generation endpoint with real PubChem/RDKit but mocked LLM."""
    request_data = {
        "prompt": WATER_DECOMP_PROMPT,
        "canvas_width": MOCK_CANVAS_WIDTH, # Use defined mock canvas size
        "canvas_height": MOCK_CANVAS_HEIGHT # Use defined mock canvas size
    }

    response = client.post("/prompt/generate-molecule-diagram/", json=request_data)
    assert response.status_code == 200

    data = response.json()
    assert isinstance(data, dict)
    assert "diagram_image" in data
    assert "diagram_plan" in data
    assert data["status"] == "completed"

    # Verify diagram image is valid SVG
    assert verify_svg_content(data["diagram_image"])

    # Verify plan structure
    plan = data["diagram_plan"]
    assert isinstance(plan, dict)
    assert "plan" in plan
    assert "molecule_list" in plan
    assert "arrows" in plan
    assert len(plan["molecule_list"]) == 3  # H2O, H2, O2
    assert len(plan["arrows"]) == 1  # One reaction arrow
    assert plan.get("canvas_width") == MOCK_CANVAS_WIDTH
    assert plan.get("canvas_height") == MOCK_CANVAS_HEIGHT

def test_real_molecule_rendering():
    """Test actual molecule rendering with RDKit."""
    # Create a simple water molecule
    mol = Chem.MolFromSmiles("O")  # Water
    Chem.AddHs(mol)  # Add hydrogens
    
    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)
    
    # Draw molecule
    img = Draw.MolToImage(mol)
    assert img.size[0] > 0 and img.size[1] > 0

def test_complex_reaction(mock_openai_provider_structured):
    """Test a more complex reaction diagram."""
    COMPLEX_CANVAS_WIDTH = 1200
    COMPLEX_CANVAS_HEIGHT = 800

    # Mock a specific DiagramPlan for this complex reaction scenario
    complex_diagram_plan_mock = DiagramPlan(
        plan="Plan for methane combustion", 
        molecule_list=[
            MoleculePlacement(molecule="CH4", x=100, y=100, label="CH₄"),
            MoleculePlacement(molecule="O2", x=300, y=100, label="2O₂"),
            MoleculePlacement(molecule="CO2", x=500, y=100, label="CO₂"),
            MoleculePlacement(molecule="H2O", x=700, y=100, label="2H₂O")
        ],
        arrows=[],
        canvas_width=COMPLEX_CANVAS_WIDTH, # Include canvas dimensions in the mock
        canvas_height=COMPLEX_CANVAS_HEIGHT
    )
    mock_openai_provider_structured.return_value = complex_diagram_plan_mock

    request_data = {
        "prompt": "Show the reaction of methane combustion (CH4 + O2 -> CO2 + H2O)",
        "canvas_width": COMPLEX_CANVAS_WIDTH,
        "canvas_height": COMPLEX_CANVAS_HEIGHT
    }
    response = client.post("/prompt/generate-molecule-diagram/", json=request_data)
    assert response.status_code == 200
    data = response.json()
    assert verify_svg_content(data["diagram_image"])
    
    plan = data["diagram_plan"]
    molecules_in_plan = [m["molecule"] for m in plan.get("molecule_list", [])]
    expected_molecules = ["CH4", "O2", "CO2", "H2O"]
    assert all(mol in molecules_in_plan for mol in expected_molecules)
    assert plan.get("plan") == "Plan for methane combustion"
    assert plan.get("canvas_width") == COMPLEX_CANVAS_WIDTH
    assert plan.get("canvas_height") == COMPLEX_CANVAS_HEIGHT

def test_invalid_prompt():
    """Test endpoint behavior with empty prompt."""
    request_data = {
        "prompt": "",
        "canvas_width": 960,
        "canvas_height": 640
    }

    response = client.post("/prompt/generate-molecule-diagram/", json=request_data)
    assert response.status_code == 422  # Validation error

def test_llm_error_handling(mock_openai_provider_error):
    """Test handling of LLM service errors."""
    mock_openai_provider_error.side_effect = Exception("LLM provider simulated error")
    
    request_data = {
        "prompt": WATER_DECOMP_PROMPT,
        "canvas_width": 960,
        "canvas_height": 640
    }

    response = client.post("/prompt/generate-molecule-diagram/", json=request_data)
    assert response.status_code == 500
    # Check if the *specific* simulated error message is part of the detail
    # or if the generic "Unexpected server error" is acceptable.
    # For now, let's check for the generic prefix + part of our message.
    json_response = response.json()
    assert "detail" in json_response
    assert "Unexpected server error: LLM provider simulated error" in json_response["detail"]

def test_custom_canvas_size(mock_openai_provider_structured):
    """Test endpoint with non-default canvas dimensions."""
    CUSTOM_CANVAS_WIDTH = 1200
    CUSTOM_CANVAS_HEIGHT = 800

    # Ensure the mock returns a DiagramPlan that includes these custom dimensions
    # This is important because the endpoint logic now relies on the LLM-returned plan object (even if mocked)
    # to eventually form the basis of the response plan.
    custom_plan_mock = DiagramPlan(
        plan=EXPECTED_DIAGRAM_PLAN_OBJECT.plan,
        molecule_list=EXPECTED_DIAGRAM_PLAN_OBJECT.molecule_list,
        arrows=EXPECTED_DIAGRAM_PLAN_OBJECT.arrows,
        canvas_width=CUSTOM_CANVAS_WIDTH,
        canvas_height=CUSTOM_CANVAS_HEIGHT
    )
    mock_openai_provider_structured.return_value = custom_plan_mock

    request_data = {
        "prompt": WATER_DECOMP_PROMPT,
        "canvas_width": CUSTOM_CANVAS_WIDTH,
        "canvas_height": CUSTOM_CANVAS_HEIGHT
    }
    response = client.post("/prompt/generate-molecule-diagram/", json=request_data)
    assert response.status_code == 200
    data = response.json()
    assert "diagram_image" in data
    
    plan = data["diagram_plan"]
    assert plan.get("canvas_width") == CUSTOM_CANVAS_WIDTH
    assert plan.get("canvas_height") == CUSTOM_CANVAS_HEIGHT
    
    svg_str = data["diagram_image"]
    width_match = re.search(r'width="(\d+)"', svg_str)
    height_match = re.search(r'height="(\d+)"', svg_str)
    assert width_match and int(width_match.group(1)) == CUSTOM_CANVAS_WIDTH
    assert height_match and int(height_match.group(1)) == CUSTOM_CANVAS_HEIGHT

@pytest.mark.skipif(not os.environ.get("OPENAI_API_KEY"), reason="OPENAI_API_KEY not set, skipping real API call test")
def test_real_openai_diagram_and_download():
    """Tests diagram generation with a real OpenAI call and downloads the SVG."""
    request_data = {
        "prompt": "Draw a simple water molecule H2O", # A very simple prompt for the real call
        "canvas_width": 300,
        "canvas_height": 200,
        "model": "gpt-4o" # Explicitly suggest a model known to work well with structured output if possible
                           # Or ensure your default LLMService config points to a capable OpenAI model
    }

    print("\nAttempting real OpenAI API call for diagram generation...")
    response = client.post("/prompt/generate-molecule-diagram/", json=request_data)

    print(f"Response Status Code: {response.status_code}")
    try:
        data = response.json()
        print(f"Response JSON: {json.dumps(data, indent=2)}")
    except json.JSONDecodeError:
        print(f"Response Text (not JSON): {response.text}")
        data = None

    assert response.status_code == 200, f"API call failed with status {response.status_code} and text: {response.text}"
    
    assert data is not None, "Response was not valid JSON."
    assert isinstance(data, dict)
    assert data.get("status") == "completed", f"Diagram generation status was not 'completed'. Full response: {data}"
    
    diagram_image_svg = data.get("diagram_image")
    assert diagram_image_svg, "No diagram_image found in response."
    assert isinstance(diagram_image_svg, str), "diagram_image is not a string."
    assert "<svg" in diagram_image_svg.strip().lower(), "diagram_image does not look like an SVG."

    # Save the SVG to a file in the root directory
    # The workspace root is the parent of the 'api' directory
    # Assuming the test is run from the workspace root.
    file_path = "test_diagram_output.svg"
    try:
        with open(file_path, "w") as f:
            f.write(diagram_image_svg)
        print(f"Diagram SVG saved to: {os.path.abspath(file_path)}")
    except Exception as e:
        pytest.fail(f"Failed to write SVG to file: {e}")

    # Basic check that file was created and is not empty
    assert os.path.exists(file_path), f"SVG file was not created at {file_path}"
    assert os.path.getsize(file_path) > 0, f"SVG file {file_path} is empty"
