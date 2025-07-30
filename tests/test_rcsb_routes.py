"""Integration tests for RCSB routes."""


class TestRCSBRoutes:
    """Test cases for /rcsb routes."""

    def test_fetch_structure_pdb_format(self, client):
        """Test fetching structure in PDB format."""
        request_data = {"identifier": "1ubq", "format": "pdb"}

        response = client.post("/rcsb/fetch-structure/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "data" in data
        assert len(data["data"]) > 0

    def test_fetch_structure_cif_format(self, client):
        """Test fetching structure in CIF format."""
        request_data = {"identifier": "1ubq", "format": "cif"}

        response = client.post("/rcsb/fetch-structure/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "data" in data
        assert len(data["data"]) > 0

    def test_fetch_structure_invalid_identifier(self, client):
        """Test fetching structure with invalid identifier."""
        request_data = {"identifier": "invalid_id_12345", "format": "pdb"}

        response = client.post("/rcsb/fetch-structure/", json=request_data)

        assert response.status_code == 500  # Should handle gracefully

    def test_fetch_model_alphafold(self, client):
        """Test fetching AlphaFold model."""
        request_data = {
            "uniprot_id": "P69905",  # Human hemoglobin
            "format": "pdb",
        }

        response = client.post("/rcsb/fetch-model/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "data" in data

    def test_fetch_model_cif_format(self, client):
        """Test fetching AlphaFold model in CIF format."""
        request_data = {"uniprot_id": "P69905", "format": "cif"}

        response = client.post("/rcsb/fetch-model/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "data" in data

    def test_entry_metadata(self, client):
        """Test fetching entry metadata."""
        identifier = "1ubq"

        response = client.get(f"/rcsb/entry/{identifier}")

        assert response.status_code == 200
        data = response.json()
        assert "metadata" in data
        assert isinstance(data["metadata"], dict)

    def test_entry_metadata_invalid_id(self, client):
        """Test fetching metadata for invalid entry."""
        identifier = "invalid_entry_id"

        response = client.get(f"/rcsb/entry/{identifier}")

        assert response.status_code == 500

    def test_sequence_annotations(self, client):
        """Test fetching sequence annotations."""
        identifier = "1ubq"

        response = client.get(f"/rcsb/annotations/{identifier}")

        assert response.status_code == 200
        data = response.json()
        assert "annotations" in data
        assert isinstance(data["annotations"], dict)

    def test_sequence_annotations_invalid_id(self, client):
        """Test fetching annotations for invalid identifier."""
        identifier = "invalid_id"

        response = client.get(f"/rcsb/annotations/{identifier}")

        assert response.status_code == 500

    def test_computed_model(self, client):
        """Test fetching computed model via GraphQL."""
        request_data = {"identifier": "1ubq", "model_id": "1"}

        response = client.post("/rcsb/computed-model/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "metadata" in data
        assert isinstance(data["metadata"], dict)

    def test_computed_model_invalid_data(self, client):
        """Test computed model with invalid data."""
        request_data = {"identifier": "invalid", "model_id": "999"}

        response = client.post("/rcsb/computed-model/", json=request_data)

        assert response.status_code == 500

    def test_fetch_esm_model(self, client):
        """Test fetching ESM model."""
        request_data = {"uniprot_id": "P69905", "format": "pdb"}

        response = client.post("/rcsb/fetch-esm-model/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "data" in data

    def test_fetch_esm_model_cif(self, client):
        """Test fetching ESM model in CIF format."""
        request_data = {"uniprot_id": "P69905", "format": "cif"}

        response = client.post("/rcsb/fetch-esm-model/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "data" in data

    def test_upload_structure_pdb(self, client):
        """Test uploading a PDB structure."""
        sample_pdb_data = """HEADER    TRANSFERASE/DNA                 20-JUL-95   1ABC
ATOM      1  N   ALA A   1      20.154  16.967  14.421  1.00 20.00           N
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C
END"""

        request_data = {"data": sample_pdb_data, "filename": "test_structure.pdb"}

        response = client.post("/rcsb/upload-structure/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "upload_id" in data
        assert len(data["upload_id"]) > 0

    def test_upload_structure_default_filename(self, client):
        """Test uploading structure with default filename."""
        sample_data = "ATOM      1  N   ALA A   1      20.154  16.967  14.421  1.00 20.00           N"

        request_data = {"data": sample_data}

        response = client.post("/rcsb/upload-structure/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "upload_id" in data

    def test_pairwise_alignment(self, client):
        """Test pairwise structure alignment."""
        request_data = {"identifier1": "1ubq", "identifier2": "1crn"}

        response = client.post("/rcsb/align/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "metadata" in data
        assert isinstance(data["metadata"], dict)

    def test_pairwise_alignment_same_structure(self, client):
        """Test alignment of structure with itself."""
        request_data = {"identifier1": "1ubq", "identifier2": "1ubq"}

        response = client.post("/rcsb/align/", json=request_data)

        assert response.status_code == 200
        data = response.json()
        assert "metadata" in data

    def test_group_entries(self, client):
        """Test fetching group entries."""
        group_id = "G_1002001"

        response = client.get(f"/rcsb/group/{group_id}")

        assert response.status_code == 200
        data = response.json()
        assert "metadata" in data
        assert isinstance(data["metadata"], dict)

    def test_group_entries_invalid_id(self, client):
        """Test fetching invalid group entries."""
        group_id = "invalid_group_id"

        response = client.get(f"/rcsb/group/{group_id}")

        assert response.status_code == 500

    def test_feature_annotations(self, client):
        """Test fetching feature annotations."""
        identifier = "1ubq"

        response = client.get(f"/rcsb/feature-annotations/{identifier}")

        assert response.status_code == 200
        data = response.json()
        assert "annotations" in data
        assert isinstance(data["annotations"], dict)

    def test_feature_annotations_invalid_id(self, client):
        """Test fetching feature annotations for invalid identifier."""
        identifier = "invalid_feature_id"

        response = client.get(f"/rcsb/feature-annotations/{identifier}")

        assert response.status_code == 500

    def test_fetch_structure_multiple_formats(self, client):
        """Test fetching the same structure in different formats."""
        identifier = "2itx"

        # Test PDB format
        pdb_response = client.post("/rcsb/fetch-structure/", json={"identifier": identifier, "format": "pdb"})

        # Test CIF format
        cif_response = client.post("/rcsb/fetch-structure/", json={"identifier": identifier, "format": "cif"})

        assert pdb_response.status_code == 200
        assert cif_response.status_code == 200

        pdb_data = pdb_response.json()
        cif_data = cif_response.json()

        assert "data" in pdb_data
        assert "data" in cif_data
        assert len(pdb_data["data"]) > 0
        assert len(cif_data["data"]) > 0
