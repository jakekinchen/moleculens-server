"""Golden example integration test.

This test runs against a live API + worker stack and validates the complete
orbital computation pipeline using the water molecule test case.
"""

import base64
import gzip
import json
import time
from pathlib import Path
from typing import Any

import httpx
import numpy as np
import pytest


def decode_mesh_array(b64_data: str, dtype: str = "float32") -> np.ndarray:
    """Decode a compressed mesh array."""
    compressed = base64.b64decode(b64_data)
    decompressed = gzip.decompress(compressed)
    np_dtype = np.float32 if dtype == "float32" else np.uint32
    return np.frombuffer(decompressed, dtype=np_dtype)


class TestGoldenExample:
    """Golden example integration tests."""

    @pytest.fixture
    def expected_metrics(self, golden_dir: Path) -> dict[str, Any]:
        """Load expected metrics from golden file."""
        metrics_file = golden_dir / "expected_metrics.json"
        return json.loads(metrics_file.read_text())

    def test_health_check(self, api_url: str) -> None:
        """Test that the API is healthy."""
        response = httpx.get(f"{api_url}/health", timeout=10)
        assert response.status_code == 200

        data = response.json()
        assert data["status"] in ("healthy", "degraded")
        assert "version" in data

    def test_compute_water_orbitals(
        self,
        api_url: str,
        water_sdf: str,
        expected_metrics: dict[str, Any],
        artifacts_dir: Path,
    ) -> None:
        """Test computing water molecule orbitals end-to-end."""
        # Submit computation request
        request_data = {
            "sdfContent": water_sdf,
            "method": expected_metrics["method"],
            "basis": expected_metrics["basis"],
            "orbitals": expected_metrics["expectedOrbitals"]
            + (["density"] if expected_metrics["expectedDensity"] else []),
            "gridSpacing": expected_metrics["gridSpacing"],
            "isovalue": expected_metrics["isovalue"],
        }

        response = httpx.post(
            f"{api_url}/v1/orbitals/compute",
            json=request_data,
            timeout=30,
        )
        assert response.status_code == 200, f"Submit failed: {response.text}"

        job_data = response.json()
        job_id = job_data["jobId"]

        assert "jobId" in job_data
        assert "cacheKey" in job_data
        assert "status" in job_data
        assert job_data["status"] in ("pending", "queued", "running", "done")

        # Poll for completion
        max_wait = 300  # 5 minutes max
        poll_interval = 2
        elapsed = 0

        while elapsed < max_wait:
            response = httpx.get(
                f"{api_url}/v1/orbitals/jobs/{job_id}",
                timeout=10,
            )
            assert response.status_code == 200

            job_data = response.json()
            status = job_data["status"]

            if status == "done":
                break
            elif status == "error":
                pytest.fail(f"Job failed: {job_data.get('errorMessage', 'Unknown error')}")

            time.sleep(poll_interval)
            elapsed += poll_interval

        assert status == "done", f"Job did not complete in {max_wait}s (status: {status})"

        # Validate result structure
        result = job_data.get("result")
        assert result is not None, "No result in completed job"

        # Save result to artifacts
        result_file = artifacts_dir / "result.json"
        result_file.write_text(json.dumps(job_data, indent=2))

        # Validate orbitals
        orbitals = result.get("orbitals", {})
        tolerances = expected_metrics["tolerances"]

        for orbital_name in expected_metrics["expectedOrbitals"]:
            assert orbital_name in orbitals, f"Missing orbital: {orbital_name}"
            orbital = orbitals[orbital_name]

            # Check for positive lobe (should exist for most orbitals)
            if orbital.get("positive"):
                self._validate_mesh(
                    orbital["positive"],
                    tolerances,
                    f"{orbital_name}_positive",
                    artifacts_dir,
                )

            # Check for negative lobe
            if orbital.get("negative"):
                self._validate_mesh(
                    orbital["negative"],
                    tolerances,
                    f"{orbital_name}_negative",
                    artifacts_dir,
                )

            # Check energy
            if orbital_name == "homo":
                energy = orbital.get("energyEv")
                if energy is not None:
                    assert (
                        tolerances["energyHomo"]["min"] <= energy <= tolerances["energyHomo"]["max"]
                    ), f"HOMO energy {energy} eV out of range"

            if orbital_name == "lumo":
                energy = orbital.get("energyEv")
                if energy is not None:
                    assert (
                        tolerances["energyLumo"]["min"] <= energy <= tolerances["energyLumo"]["max"]
                    ), f"LUMO energy {energy} eV out of range"

        # Validate density if expected
        if expected_metrics["expectedDensity"]:
            density = result.get("density")
            if density:
                self._validate_mesh(density, tolerances, "density", artifacts_dir)

        # Validate metadata
        meta = result.get("meta", {})
        assert meta.get("method") == expected_metrics["method"]
        assert expected_metrics["basis"].replace("*", "star") in meta.get("basis", "").lower()
        assert meta.get("computeTimeMs") is not None

        # Generate preview images
        self._generate_preview_images(orbitals, result.get("density"), artifacts_dir)

    def _validate_mesh(
        self,
        mesh: dict[str, Any],
        tolerances: dict[str, Any],
        name: str,
        artifacts_dir: Path,
    ) -> None:
        """Validate mesh data and check bounding box."""
        assert "vertices" in mesh
        assert "normals" in mesh
        assert "indices" in mesh
        assert "vertexCount" in mesh
        assert "triangleCount" in mesh

        # Check counts are within tolerance
        vertex_count = mesh["vertexCount"]
        triangle_count = mesh["triangleCount"]

        assert (
            tolerances["vertexCount"]["min"] <= vertex_count <= tolerances["vertexCount"]["max"]
        ), f"{name}: vertex count {vertex_count} out of range"
        assert (
            tolerances["triangleCount"]["min"]
            <= triangle_count
            <= tolerances["triangleCount"]["max"]
        ), f"{name}: triangle count {triangle_count} out of range"

        # Decode and validate vertices
        vertices = decode_mesh_array(mesh["vertices"], "float32")
        assert len(vertices) == vertex_count * 3, f"{name}: vertex array size mismatch"

        # Reshape to (n, 3)
        vertices = vertices.reshape(-1, 3)

        # Check bounding box
        bbox_min = vertices.min(axis=0)
        bbox_max = vertices.max(axis=0)
        bbox_size = bbox_max - bbox_min

        for dim, size in enumerate(bbox_size):
            assert (
                tolerances["boundingBox"]["minDimension"]
                <= size
                <= tolerances["boundingBox"]["maxDimension"]
            ), f"{name}: dimension {dim} size {size} out of range"

        # Check that vertices are finite
        assert np.all(np.isfinite(vertices)), f"{name}: non-finite vertices found"

        # Decode and validate normals
        normals = decode_mesh_array(mesh["normals"], "float32")
        assert len(normals) == vertex_count * 3, f"{name}: normal array size mismatch"

        normals = normals.reshape(-1, 3)
        assert np.all(np.isfinite(normals)), f"{name}: non-finite normals found"

        # Decode and validate indices
        indices = decode_mesh_array(mesh["indices"], "uint32")
        assert len(indices) == triangle_count * 3, f"{name}: index array size mismatch"

        # All indices should be valid
        assert np.all(indices < vertex_count), f"{name}: invalid vertex indices"

        # Save mesh stats
        stats = {
            "vertexCount": vertex_count,
            "triangleCount": triangle_count,
            "boundingBoxMin": bbox_min.tolist(),
            "boundingBoxMax": bbox_max.tolist(),
            "boundingBoxSize": bbox_size.tolist(),
        }
        stats_file = artifacts_dir / f"{name}_stats.json"
        stats_file.write_text(json.dumps(stats, indent=2))

    def _generate_preview_images(
        self,
        orbitals: dict[str, Any],
        density: dict[str, Any] | None,
        artifacts_dir: Path,
    ) -> None:
        """Generate simple 2D projection preview images."""
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            return  # Skip if matplotlib not available

        for orbital_name, orbital in orbitals.items():
            for lobe_name in ["positive", "negative"]:
                lobe = orbital.get(lobe_name)
                if lobe is None:
                    continue

                vertices = decode_mesh_array(lobe["vertices"], "float32").reshape(-1, 3)

                fig, axes = plt.subplots(1, 3, figsize=(12, 4))

                # XY projection
                axes[0].scatter(vertices[:, 0], vertices[:, 1], s=1, alpha=0.5)
                axes[0].set_xlabel("X (Å)")
                axes[0].set_ylabel("Y (Å)")
                axes[0].set_title(f"{orbital_name} {lobe_name} - XY")
                axes[0].set_aspect("equal")

                # XZ projection
                axes[1].scatter(vertices[:, 0], vertices[:, 2], s=1, alpha=0.5)
                axes[1].set_xlabel("X (Å)")
                axes[1].set_ylabel("Z (Å)")
                axes[1].set_title(f"{orbital_name} {lobe_name} - XZ")
                axes[1].set_aspect("equal")

                # YZ projection
                axes[2].scatter(vertices[:, 1], vertices[:, 2], s=1, alpha=0.5)
                axes[2].set_xlabel("Y (Å)")
                axes[2].set_ylabel("Z (Å)")
                axes[2].set_title(f"{orbital_name} {lobe_name} - YZ")
                axes[2].set_aspect("equal")

                plt.tight_layout()
                plt.savefig(artifacts_dir / f"preview_{orbital_name}_{lobe_name}.png", dpi=100)
                plt.close()

        # Density preview
        if density:
            vertices = decode_mesh_array(density["vertices"], "float32").reshape(-1, 3)

            fig, axes = plt.subplots(1, 3, figsize=(12, 4))

            axes[0].scatter(vertices[:, 0], vertices[:, 1], s=1, alpha=0.5, c="blue")
            axes[0].set_xlabel("X (Å)")
            axes[0].set_ylabel("Y (Å)")
            axes[0].set_title("Density - XY")
            axes[0].set_aspect("equal")

            axes[1].scatter(vertices[:, 0], vertices[:, 2], s=1, alpha=0.5, c="blue")
            axes[1].set_xlabel("X (Å)")
            axes[1].set_ylabel("Z (Å)")
            axes[1].set_title("Density - XZ")
            axes[1].set_aspect("equal")

            axes[2].scatter(vertices[:, 1], vertices[:, 2], s=1, alpha=0.5, c="blue")
            axes[2].set_xlabel("Y (Å)")
            axes[2].set_ylabel("Z (Å)")
            axes[2].set_title("Density - YZ")
            axes[2].set_aspect("equal")

            plt.tight_layout()
            plt.savefig(artifacts_dir / "preview_density.png", dpi=100)
            plt.close()

    def test_job_deduplication(self, api_url: str, water_sdf: str) -> None:
        """Test that duplicate requests return the same job."""
        request_data = {
            "sdfContent": water_sdf,
            "method": "scf",
            "basis": "sto-3g",
            "orbitals": ["homo"],
            "gridSpacing": 0.3,
            "isovalue": 0.05,
        }

        # Submit first request
        response1 = httpx.post(
            f"{api_url}/v1/orbitals/compute",
            json=request_data,
            timeout=30,
        )
        assert response1.status_code == 200
        job1 = response1.json()

        # Submit duplicate request
        response2 = httpx.post(
            f"{api_url}/v1/orbitals/compute",
            json=request_data,
            timeout=30,
        )
        assert response2.status_code == 200
        job2 = response2.json()

        # Should have same cache key
        assert job1["cacheKey"] == job2["cacheKey"]

    def test_invalid_sdf_returns_error(self, api_url: str) -> None:
        """Test that invalid SDF content returns an error."""
        request_data = {
            "sdfContent": "not valid sdf content",
            "method": "scf",
            "basis": "sto-3g",
        }

        response = httpx.post(
            f"{api_url}/v1/orbitals/compute",
            json=request_data,
            timeout=30,
        )

        # Should return 400 for invalid input
        assert response.status_code == 400

    def test_missing_job_returns_404(self, api_url: str) -> None:
        """Test that missing job ID returns 404."""
        response = httpx.get(
            f"{api_url}/v1/orbitals/jobs/nonexistent-job-id",
            timeout=10,
        )
        assert response.status_code == 404
