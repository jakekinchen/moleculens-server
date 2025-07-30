"""Integration tests for graphic routes."""


class TestGraphicRoutes:
    """Test cases for /graphic routes."""

    def test_plan_graphic_success(self, client):
        """Test successful graphic plan generation."""
        request_data = {
            "brief": "Show water H2O splitting into hydrogen H2 and oxygen O2",
            "width": 800,
            "height": 400,
            "model_name": "o3-mini",
        }

        response = client.post("/graphic/plan", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "yaml_spec" in data
        assert data["status"] == "completed"
        assert data["error"] is None
        assert len(data["yaml_spec"]) > 0

    def test_plan_graphic_empty_brief(self, client):
        """Test graphic plan with empty brief."""
        request_data = {"brief": "", "width": 800, "height": 400}

        response = client.post("/graphic/plan", json=request_data)

        assert response.status_code == 422

    def test_plan_graphic_with_context(self, client):
        """Test graphic plan with custom context and theme."""
        request_data = {
            "brief": "Calvin cycle photosynthesis process",
            "context": "Educational biochemistry diagram",
            "theme": "Modern scientific illustration",
            "width": 1200,
            "height": 800,
            "sections": "Input molecules, Enzyme reactions, Output products",
            "notes": "Use clear arrows and labels",
        }

        response = client.post("/graphic/plan", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "completed"
        assert "yaml_spec" in data

    def test_validate_graphic_valid_spec(self, client, sample_yaml_spec):
        """Test validation of a valid YAML spec."""
        request_data = {"yaml_spec": sample_yaml_spec}

        response = client.post("/graphic/validate", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert len(data["errors"]) == 0
        assert data["status"] == "completed"

    def test_validate_graphic_invalid_yaml(self, client):
        """Test validation of invalid YAML."""
        request_data = {"yaml_spec": "invalid: yaml: content: ["}

        response = client.post("/graphic/validate", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False
        assert len(data["errors"]) > 0
        assert "YAML syntax error" in data["errors"][0]

    def test_validate_graphic_empty_spec(self, client):
        """Test validation with empty spec."""
        request_data = {"yaml_spec": ""}

        response = client.post("/graphic/validate", json=request_data)

        assert response.status_code == 422

    def test_render_graphic_svg(self, client, sample_yaml_spec):
        """Test rendering graphic as SVG."""
        request_data = {"yaml_spec": sample_yaml_spec, "output_format": "svg"}

        response = client.post("/graphic/render", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "svg_content" in data
        assert data["status"] == "completed"
        assert data["png_base64"] is None
        assert len(data["svg_content"]) > 0

    def test_render_graphic_png(self, client, sample_yaml_spec):
        """Test rendering graphic as PNG."""
        request_data = {"yaml_spec": sample_yaml_spec, "output_format": "png"}

        response = client.post("/graphic/render", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "svg_content" in data
        assert data["png_base64"] is not None
        assert data["status"] == "completed"

    def test_render_graphic_both_formats(self, client, sample_yaml_spec):
        """Test rendering graphic in both SVG and PNG formats."""
        request_data = {"yaml_spec": sample_yaml_spec, "output_format": "both"}

        response = client.post("/graphic/render", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "svg_content" in data
        assert data["png_base64"] is not None
        assert data["status"] == "completed"
        assert len(data["svg_content"]) > 0

    def test_render_graphic_invalid_yaml(self, client):
        """Test rendering with invalid YAML."""
        request_data = {"yaml_spec": "invalid: yaml: ["}

        response = client.post("/graphic/render", json=request_data)

        assert response.status_code == 422

    def test_make_graphic_success(self, client):
        """Test full pipeline graphic creation."""
        request_data = {
            "brief": "Show glucose C6H12O6 converting to pyruvate",
            "width": 960,
            "height": 640,
            "output_format": "both",
        }

        response = client.post("/graphic/make", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "yaml_spec" in data
        assert "svg_content" in data
        assert data["status"] == "completed"
        assert len(data["yaml_spec"]) > 0
        assert len(data["svg_content"]) > 0

    def test_make_graphic_svg_only(self, client):
        """Test full pipeline with SVG output only."""
        request_data = {"brief": "ATP hydrolysis reaction", "width": 800, "height": 400, "output_format": "svg"}

        response = client.post("/graphic/make", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["png_base64"] is None
        assert "svg_content" in data
        assert data["status"] == "completed"

    def test_make_graphic_empty_brief(self, client):
        """Test full pipeline with empty brief."""
        request_data = {"brief": "", "width": 800, "height": 400}

        response = client.post("/graphic/make", json=request_data)

        assert response.status_code == 422

    def test_get_job_status(self, client):
        """Test job status endpoint."""
        job_id = "test-job-123"

        response = client.get(f"/graphic/job/{job_id}")

        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] == job_id
        assert data["status"] == "completed"
        assert "message" in data

    def test_make_graphic_custom_dimensions(self, client):
        """Test graphic creation with custom dimensions."""
        request_data = {
            "brief": "Photosynthesis light and dark reactions",
            "width": 1400,
            "height": 900,
            "model_name": "o3-mini",
            "output_format": "svg",
        }

        response = client.post("/graphic/make", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "completed"
        assert "yaml_spec" in data
        assert "svg_content" in data
