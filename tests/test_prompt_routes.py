"""Integration tests for prompt routes."""


class TestPromptRoutes:
    """Test cases for /prompt routes."""

    def test_generate_from_pubchem_success(self, client):
        """Test successful molecule name generation from PubChem."""
        request_data = {"prompt": "Show me caffeine and aspirin molecules"}

        response = client.post("/prompt/generate-from-pubchem/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "molecule_names" in data
        assert isinstance(data["molecule_names"], list)
        assert len(data["molecule_names"]) > 0

    def test_generate_from_pubchem_with_model(self, client):
        """Test molecule generation with specific model."""
        request_data = {"prompt": "Generate molecules for diabetes treatment", "model": "gpt-4"}

        response = client.post("/prompt/generate-from-pubchem/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "molecule_names" in data
        assert isinstance(data["molecule_names"], list)

    def test_generate_from_pubchem_non_molecular(self, client):
        """Test rejection of non-molecular prompts."""
        request_data = {"prompt": "Tell me about the weather today"}

        response = client.post("/prompt/generate-from-pubchem/", json=request_data)

        # Should return 400 for non-molecular prompt
        assert response.status_code == 400
        assert "Non-molecular prompt rejected" in response.json()["detail"]

    def test_validate_scientific_molecular(self, client):
        """Test validation of molecular/scientific prompt."""
        request_data = {"prompt": "Show me the structure of glucose C6H12O6"}

        response = client.post("/prompt/validate-scientific/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["is_molecular"] is True

    def test_validate_scientific_non_molecular(self, client):
        """Test validation of non-scientific prompt."""
        request_data = {"prompt": "What's the weather like today?"}

        response = client.post("/prompt/validate-scientific/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["is_molecular"] is False

    def test_validate_scientific_chemistry_terms(self, client):
        """Test validation with chemistry-related terms."""
        request_data = {"prompt": "Explain protein folding and enzyme catalysis"}

        response = client.post("/prompt/validate-scientific/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["is_molecular"] is True

    def test_generate_molecule_diagram_success(self, client):
        """Test successful molecule diagram generation."""
        request_data = {
            "prompt": "Show water H2O decomposing into hydrogen H2 and oxygen O2",
            "canvas_width": 800,
            "canvas_height": 400,
        }

        response = client.post("/prompt/generate-molecule-diagram/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "diagram_image" in data
        assert "diagram_plan" in data
        assert data["status"] == "completed"
        assert data["diagram_plan"]["canvas_width"] == 800
        assert data["diagram_plan"]["canvas_height"] == 400
        assert len(data["diagram_plan"]["molecule_list"]) > 0

    def test_generate_molecule_diagram_empty_prompt(self, client):
        """Test diagram generation with empty prompt."""
        request_data = {"prompt": "", "canvas_width": 800, "canvas_height": 400}

        response = client.post("/prompt/generate-molecule-diagram/", json=request_data)

        assert response.status_code == 422

    def test_generate_molecule_diagram_custom_dimensions(self, client):
        """Test diagram generation with custom canvas dimensions."""
        request_data = {
            "prompt": "ATP hydrolysis to ADP and phosphate",
            "canvas_width": 1200,
            "canvas_height": 600,
            "model": "gpt-4",
        }

        response = client.post("/prompt/generate-molecule-diagram/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["diagram_plan"]["canvas_width"] == 1200
        assert data["diagram_plan"]["canvas_height"] == 600
        assert data["status"] == "completed"

    def test_fetch_molecule_data_success(self, client):
        """Test fetching molecule data."""
        request_data = {"query": "caffeine"}

        response = client.post("/prompt/fetch-molecule-data/", json=request_data)

        assert response.status_code == 200
        # Response structure depends on PubChemClient implementation

    def test_fetch_molecule_data_chemical_formula(self, client):
        """Test fetching molecule data with chemical formula."""
        request_data = {
            "query": "C8H10N4O2"  # Caffeine formula
        }

        response = client.post("/prompt/fetch-molecule-data/", json=request_data)

        assert response.status_code == 200

    def test_fetch_molecule_2d_success(self, client):
        """Test fetching 2D molecule data."""
        request_data = {"query": "aspirin"}

        response = client.post("/prompt/fetch-molecule-2d/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["query"] == "aspirin"
        assert "status" in data

    def test_fetch_molecule_layout_success(self, client):
        """Test fetching molecule layout data."""
        request_data = {
            "molecules": [
                {"query": "glucose", "box": {"x": 100.0, "y": 100.0, "width": 200.0, "height": 150.0}},
                {"query": "fructose", "box": {"x": 400.0, "y": 100.0, "width": 200.0, "height": 150.0}},
            ]
        }

        response = client.post("/prompt/fetch-molecule-layout/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "molecules" in data
        assert len(data["molecules"]) == 2

    def test_fetch_molecule_layout_empty_list(self, client):
        """Test fetching molecule layout with empty molecule list."""
        request_data = {"molecules": []}

        response = client.post("/prompt/fetch-molecule-layout/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "molecules" in data
        assert len(data["molecules"]) == 0

    def test_sdf_to_pdb_conversion_success(self, client, sample_sdf_data):
        """Test successful SDF to PDB conversion."""
        request_data = {"sdf": sample_sdf_data}

        response = client.post("/prompt/sdf-to-pdb/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "pdb_data" in data
        assert len(data["pdb_data"]) > 0

    def test_sdf_to_pdb_invalid_sdf(self, client):
        """Test SDF to PDB conversion with invalid SDF data."""
        request_data = {"sdf": "invalid sdf data"}

        response = client.post("/prompt/sdf-to-pdb/", json=request_data)

        assert response.status_code == 400

    def test_sdf_to_pdb_empty_sdf(self, client):
        """Test SDF to PDB conversion with empty SDF data."""
        request_data = {"sdf": ""}

        response = client.post("/prompt/sdf-to-pdb/", json=request_data)

        assert response.status_code == 400

    def test_generate_molecule_diagram_complex_reaction(self, client):
        """Test diagram generation for complex biochemical reaction."""
        request_data = {
            "prompt": "Show the Calvin cycle with RuBP, CO2, and glucose formation",
            "canvas_width": 1400,
            "canvas_height": 800,
        }

        response = client.post("/prompt/generate-molecule-diagram/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "completed"
        assert len(data["diagram_plan"]["molecule_list"]) >= 3  # At least RuBP, CO2, glucose
        assert len(data["diagram_plan"]["arrows"]) > 0  # Should have connecting arrows
