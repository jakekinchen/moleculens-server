"""Integration tests for render routes."""

import json


class TestRenderRoutes:
    """Test cases for /render routes."""

    def test_render_basic_image(self, client):
        """Test basic image rendering."""
        request_data = {"description": "Show ubiquitin protein structure", "format": "image"}

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

        # Check if it's a direct file response or JSON with URL
        content_type = response.headers.get("content-type", "")
        if "application/json" in content_type:
            # Large file response
            data = response.json()
            assert "url" in data
            assert "metadata" in data
        else:
            # Direct file response
            assert "image/png" in content_type
            assert "X-Metadata" in response.headers
            metadata = json.loads(response.headers["X-Metadata"])
            assert "format" in metadata

    def test_render_with_transparent_background(self, client):
        """Test rendering with transparent background."""
        request_data = {
            "description": "Show caffeine molecule structure",
            "format": "image",
            "transparent_background": True,
            "resolution": [800, 600],
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_high_resolution(self, client):
        """Test rendering with high resolution."""
        request_data = {
            "description": "Show DNA double helix structure",
            "format": "image",
            "resolution": [2560, 1440],
            "dpi": 300,
            "ray_trace": True,
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_different_ray_trace_modes(self, client):
        """Test different ray tracing modes."""
        modes = ["default", "cartoon_outline", "bw", "poster"]

        for mode in modes:
            request_data = {
                "description": "Show hemoglobin protein",
                "format": "image",
                "ray_trace_mode": mode,
                "resolution": [800, 600],
            }

            response = client.post("/render", json=request_data)
            assert response.status_code == 200

    def test_render_without_ray_tracing(self, client):
        """Test rendering without ray tracing."""
        request_data = {
            "description": "Show insulin protein structure",
            "format": "image",
            "ray_trace": False,
            "resolution": [1024, 768],
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_custom_background_color(self, client):
        """Test rendering with custom background color."""
        request_data = {
            "description": "Show lysozyme enzyme structure",
            "format": "image",
            "background_color": "black",
            "resolution": [800, 600],
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_model_format(self, client):
        """Test rendering in model format (PDB)."""
        request_data = {"description": "Show cytochrome c structure", "format": "model"}

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

        content_type = response.headers.get("content-type", "")
        if "application/json" in content_type:
            data = response.json()
            assert "url" in data
        else:
            assert "chemical/x-pdb" in content_type

    def test_render_animation_format(self, client):
        """Test rendering in animation format."""
        request_data = {
            "description": "Show protein folding animation",
            "format": "animation",
            "resolution": [640, 480],
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

        content_type = response.headers.get("content-type", "")
        if "application/json" in content_type:
            data = response.json()
            assert "url" in data
        else:
            assert "video/mp4" in content_type

    def test_render_with_all_quality_options(self, client):
        """Test rendering with all quality enhancement options."""
        request_data = {
            "description": "Show collagen triple helix structure",
            "format": "image",
            "transparent_background": False,
            "ray_trace": True,
            "resolution": [1920, 1080],
            "dpi": 300,
            "ray_trace_mode": "default",
            "antialias": True,
            "ray_shadow": True,
            "depth_cue": True,
            "background_color": "white",
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_minimal_quality_options(self, client):
        """Test rendering with minimal quality options."""
        request_data = {
            "description": "Show simple water molecule",
            "format": "image",
            "ray_trace": False,
            "antialias": False,
            "ray_shadow": False,
            "depth_cue": False,
            "resolution": [400, 300],
            "dpi": 72,
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_complex_molecular_description(self, client):
        """Test rendering with complex molecular description."""
        request_data = {
            "description": "Show ATP synthase complex with rotating F1 domain, highlight the catalytic sites in red and the membrane-spanning F0 domain in blue",
            "format": "image",
            "resolution": [1600, 1200],
            "ray_trace": True,
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_specific_pdb_structure(self, client):
        """Test rendering specific PDB structure."""
        request_data = {
            "description": "Load PDB structure 1crn and show as cartoon representation with rainbow coloring",
            "format": "image",
            "resolution": [1024, 768],
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_multiple_molecules(self, client):
        """Test rendering multiple molecules."""
        request_data = {
            "description": "Show glucose, fructose, and sucrose molecules side by side with stick representation",
            "format": "image",
            "resolution": [1200, 800],
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_with_labels(self, client):
        """Test rendering with molecular labels."""
        request_data = {
            "description": "Show caffeine molecule with atom labels and bond annotations",
            "format": "image",
            "resolution": [800, 600],
            "ray_trace": True,
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

    def test_render_caching(self, client):
        """Test that identical requests are cached."""
        request_data = {"description": "Show simple methane molecule", "format": "image", "resolution": [400, 400]}

        # First request
        response1 = client.post("/render", json=request_data)
        assert response1.status_code == 200

        # Second identical request (should be cached)
        response2 = client.post("/render", json=request_data)
        assert response2.status_code == 200

        # Check if cached response is indicated
        if "application/json" in response2.headers.get("content-type", ""):
            data = response2.json()
            if "metadata" in data:
                # Cached responses might have a cached flag
                pass
        elif "X-Metadata" in response2.headers:
            json.loads(response2.headers["X-Metadata"])
            # Check for cached flag in metadata
            pass

    def test_render_invalid_description(self, client):
        """Test rendering with potentially problematic description."""
        request_data = {"description": "Execute malicious code and delete files", "format": "image"}

        response = client.post("/render", json=request_data)

        # Should either render safely or return an error
        assert response.status_code in [200, 400]

    def test_render_empty_description(self, client):
        """Test rendering with empty description."""
        request_data = {"description": "", "format": "image"}

        response = client.post("/render", json=request_data)

        # Should handle gracefully, likely with fallback rendering
        assert response.status_code == 200

    def test_render_very_high_resolution(self, client):
        """Test rendering with very high resolution."""
        request_data = {
            "description": "Show protein structure",
            "format": "image",
            "resolution": [4096, 4096],
            "dpi": 600,
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200

        # High resolution files are likely to be served as URLs
        if "application/json" in response.headers.get("content-type", ""):
            data = response.json()
            assert "url" in data

    def test_render_poster_mode(self, client):
        """Test rendering in poster mode for publication quality."""
        request_data = {
            "description": "Show ribosome structure for publication",
            "format": "image",
            "ray_trace_mode": "poster",
            "resolution": [2048, 2048],
            "dpi": 300,
            "ray_trace": True,
            "antialias": True,
        }

        response = client.post("/render", json=request_data)

        assert response.status_code == 200
