"""Full integration tests combining multiple routes."""


class TestFullIntegration:
    """Test cases for full workflow integration across multiple routes."""

    def test_full_molecular_workflow(self, client):
        """Test complete workflow: validate -> generate -> render."""
        # Step 1: Validate scientific prompt
        validate_request = {"prompt": "Show the structure of caffeine and its interaction with adenosine receptors"}

        validate_response = client.post("/prompt/validate-scientific/", json=validate_request)
        assert validate_response.status_code == 200
        assert validate_response.json()["is_molecular"] is True

        # Step 2: Generate molecule names
        generate_request = {"prompt": "Show caffeine and adenosine molecules"}

        generate_response = client.post("/prompt/generate-from-pubchem/", json=generate_request)
        assert generate_response.status_code == 200
        molecule_names = generate_response.json()["molecule_names"]
        assert len(molecule_names) > 0

        # Step 3: Create molecular diagram
        diagram_request = {
            "prompt": "Show caffeine binding to adenosine receptor",
            "canvas_width": 1000,
            "canvas_height": 600,
        }

        diagram_response = client.post("/prompt/generate-molecule-diagram/", json=diagram_request)
        assert diagram_response.status_code == 200
        assert diagram_response.json()["status"] == "completed"

        # Step 4: Render 3D structure
        render_request = {
            "description": "Show caffeine molecule in 3D with stick representation",
            "format": "image",
            "resolution": [800, 600],
        }

        render_response = client.post("/render", json=render_request)
        assert render_response.status_code == 200

    def test_protein_structure_workflow(self, client):
        """Test protein structure analysis workflow."""
        # Step 1: Fetch protein structure from RCSB
        structure_request = {"identifier": "1ubq", "format": "pdb"}

        structure_response = client.post("/rcsb/fetch-structure/", json=structure_request)
        assert structure_response.status_code == 200
        pdb_data = structure_response.json()["data"]
        assert len(pdb_data) > 0

        # Step 2: Get metadata for the structure
        metadata_response = client.get("/rcsb/entry/1ubq")
        assert metadata_response.status_code == 200
        metadata = metadata_response.json()["metadata"]
        assert isinstance(metadata, dict)

        # Step 3: Get sequence annotations
        annotations_response = client.get("/rcsb/annotations/1ubq")
        assert annotations_response.status_code == 200
        annotations = annotations_response.json()["annotations"]
        assert isinstance(annotations, dict)

        # Step 4: Render the structure
        render_request = {
            "description": "Show ubiquitin structure with cartoon representation and surface",
            "format": "image",
            "ray_trace": True,
            "resolution": [1024, 768],
        }

        render_response = client.post("/render", json=render_request)
        assert render_response.status_code == 200

    def test_diagram_creation_workflow(self, client):
        """Test complete diagram creation workflow."""
        # Step 1: Create diagram using graphic/make endpoint
        make_request = {
            "brief": "Show glycolysis pathway from glucose to pyruvate",
            "width": 1400,
            "height": 800,
            "output_format": "both",
        }

        make_response = client.post("/graphic/make", json=make_request)
        assert make_response.status_code == 200
        make_data = make_response.json()
        assert make_data["status"] == "completed"
        assert "yaml_spec" in make_data
        assert "svg_content" in make_data

        # Step 2: Validate the generated YAML spec
        validate_request = {"yaml_spec": make_data["yaml_spec"]}

        validate_response = client.post("/graphic/validate", json=validate_request)
        assert validate_response.status_code == 200
        assert validate_response.json()["valid"] is True

        # Step 3: Re-render with different format
        render_request = {"yaml_spec": make_data["yaml_spec"], "output_format": "png"}

        render_response = client.post("/graphic/render", json=render_request)
        assert render_response.status_code == 200
        assert render_response.json()["png_base64"] is not None

    def test_molecular_comparison_workflow(self, client):
        """Test workflow for comparing multiple molecules."""
        # Step 1: Generate molecules for comparison
        generate_request = {"prompt": "Compare glucose, fructose, and galactose structures"}

        generate_response = client.post("/prompt/generate-from-pubchem/", json=generate_request)
        assert generate_response.status_code == 200
        molecules = generate_response.json()["molecule_names"]

        # Step 2: Fetch data for each molecule
        for molecule in molecules[:3]:  # Limit to first 3 molecules
            fetch_request = {"query": molecule}

            fetch_response = client.post("/prompt/fetch-molecule-data/", json=fetch_request)
            assert fetch_response.status_code == 200

        # Step 3: Create comparative diagram
        diagram_request = {
            "prompt": f"Show structural comparison of {', '.join(molecules[:3])}",
            "canvas_width": 1200,
            "canvas_height": 600,
        }

        diagram_response = client.post("/prompt/generate-molecule-diagram/", json=diagram_request)
        assert diagram_response.status_code == 200

        # Step 4: Render 3D comparison
        render_request = {
            "description": f"Show {', '.join(molecules[:3])} molecules side by side in 3D",
            "format": "image",
            "resolution": [1600, 800],
        }

        render_response = client.post("/render", json=render_request)
        assert render_response.status_code == 200

    def test_biochemical_pathway_workflow(self, client):
        """Test creation of biochemical pathway diagrams."""
        # Step 1: Plan the pathway diagram
        plan_request = {
            "brief": "Create a detailed diagram of the citric acid cycle showing all intermediates and enzymes",
            "context": "Biochemistry education",
            "theme": "Clear metabolic pathway visualization",
            "width": 1600,
            "height": 1000,
            "sections": "Input molecules, Cycle intermediates, Enzymes, Products",
            "notes": "Include ATP, NADH, and FADH2 production points",
        }

        plan_response = client.post("/graphic/plan", json=plan_request)
        assert plan_response.status_code == 200
        yaml_spec = plan_response.json()["yaml_spec"]

        # Step 2: Validate the complex pathway spec
        validate_request = {"yaml_spec": yaml_spec}

        validate_response = client.post("/graphic/validate", json=validate_request)
        assert validate_response.status_code == 200

        # Step 3: Render the pathway
        render_request = {"yaml_spec": yaml_spec, "output_format": "both"}

        render_response = client.post("/graphic/render", json=render_request)
        assert render_response.status_code == 200

        # Step 4: Generate 3D structures for key molecules
        key_molecules = ["acetyl-CoA", "citrate", "alpha-ketoglutarate"]
        for molecule in key_molecules:
            render_3d_request = {
                "description": f"Show {molecule} molecule structure in 3D",
                "format": "image",
                "resolution": [600, 600],
            }

            render_3d_response = client.post("/render", json=render_3d_request)
            assert render_3d_response.status_code == 200

    def test_protein_analysis_workflow(self, client):
        """Test comprehensive protein analysis workflow."""
        protein_id = "1crn"  # Crambin

        # Step 1: Fetch protein structure
        structure_request = {"identifier": protein_id, "format": "pdb"}

        structure_response = client.post("/rcsb/fetch-structure/", json=structure_request)
        assert structure_response.status_code == 200

        # Step 2: Get comprehensive metadata
        metadata_response = client.get(f"/rcsb/entry/{protein_id}")
        assert metadata_response.status_code == 200

        # Step 3: Get feature annotations
        features_response = client.get(f"/rcsb/feature-annotations/{protein_id}")
        assert features_response.status_code == 200

        # Step 4: Upload and analyze custom structure
        custom_pdb = """HEADER    TEST PROTEIN                    01-JAN-00   TEST
ATOM      1  N   ALA A   1      20.154  16.967  14.421  1.00 20.00           N
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C
END"""

        upload_request = {"data": custom_pdb, "filename": "test_protein.pdb"}

        upload_response = client.post("/rcsb/upload-structure/", json=upload_request)
        assert upload_response.status_code == 200

        # Step 5: Render multiple views of the protein
        render_views = [
            "Show crambin protein with cartoon representation",
            "Show crambin protein with surface representation",
            "Show crambin protein with stick representation highlighting disulfide bonds",
        ]

        for view_description in render_views:
            render_request = {
                "description": view_description,
                "format": "image",
                "resolution": [800, 600],
                "ray_trace": True,
            }

            render_response = client.post("/render", json=render_request)
            assert render_response.status_code == 200

    def test_error_handling_workflow(self, client):
        """Test error handling across different routes."""
        # Test invalid inputs across routes

        # Invalid graphic plan
        invalid_plan_response = client.post("/graphic/plan", json={"brief": ""})
        assert invalid_plan_response.status_code == 422

        # Invalid YAML validation
        invalid_yaml_response = client.post("/graphic/validate", json={"yaml_spec": "invalid: ["})
        assert invalid_yaml_response.status_code == 200
        assert invalid_yaml_response.json()["valid"] is False

        # Invalid structure fetch
        invalid_structure_response = client.post(
            "/rcsb/fetch-structure/", json={"identifier": "invalid_id_12345", "format": "pdb"}
        )
        assert invalid_structure_response.status_code == 500

        # Invalid SDF conversion
        invalid_sdf_response = client.post("/prompt/sdf-to-pdb/", json={"sdf": "invalid"})
        assert invalid_sdf_response.status_code == 400

    def test_performance_workflow(self, client):
        """Test performance with multiple concurrent-like requests."""
        # Test multiple small requests
        requests = [
            {"description": "Show water molecule", "format": "image", "resolution": [200, 200]},
            {"description": "Show methane molecule", "format": "image", "resolution": [200, 200]},
            {"description": "Show ammonia molecule", "format": "image", "resolution": [200, 200]},
        ]

        for request_data in requests:
            response = client.post("/render", json=request_data)
            assert response.status_code == 200

        # Test caching by repeating first request
        response = client.post("/render", json=requests[0])
        assert response.status_code == 200
