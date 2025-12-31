"""Unit tests for marching cubes mesh generation."""

import numpy as np
import pytest

from moleculens.compute.cube_parser import CubeData
from moleculens.compute.marching_cubes import (
    MeshData,
    decode_mesh_array,
    generate_density_mesh,
    generate_orbital_meshes,
)


class TestMarchingCubes:
    """Tests for marching cubes mesh generation."""

    @pytest.fixture
    def simple_cube_data(self) -> CubeData:
        """Create simple cube data with a sphere-like field."""
        n = 10
        origin = np.array([-2.0, -2.0, -2.0])
        axes = np.array(
            [
                [0.4, 0.0, 0.0],
                [0.0, 0.4, 0.0],
                [0.0, 0.0, 0.4],
            ]
        )

        # Create a spherical field: positive in center, negative outside
        x = np.linspace(-2, 2, n)
        y = np.linspace(-2, 2, n)
        z = np.linspace(-2, 2, n)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        data = 1.0 - np.sqrt(X**2 + Y**2 + Z**2)

        return CubeData(
            origin=origin,
            axes=axes,
            n_points=(n, n, n),
            data=data,
            atom_positions=np.array([[0.0, 0.0, 0.0]]),
            atom_numbers=[6],
        )

    @pytest.fixture
    def orbital_cube_data(self) -> CubeData:
        """Create orbital-like cube data with positive and negative lobes."""
        n = 15
        origin = np.array([-3.0, -3.0, -3.0])
        axes = np.array(
            [
                [0.4, 0.0, 0.0],
                [0.0, 0.4, 0.0],
                [0.0, 0.0, 0.4],
            ]
        )

        # Create p-orbital-like field: positive on +z, negative on -z
        x = np.linspace(-3, 3, n)
        y = np.linspace(-3, 3, n)
        z = np.linspace(-3, 3, n)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        r = np.sqrt(X**2 + Y**2 + Z**2) + 0.1  # Avoid division by zero
        data = Z * np.exp(-r)

        return CubeData(
            origin=origin,
            axes=axes,
            n_points=(n, n, n),
            data=data,
            atom_positions=np.array([[0.0, 0.0, 0.0]]),
            atom_numbers=[6],
        )

    def test_generate_orbital_meshes(self, orbital_cube_data: CubeData) -> None:
        """Test generating orbital meshes with positive and negative lobes."""
        meshes = generate_orbital_meshes(orbital_cube_data, isovalue=0.1)

        assert "positive" in meshes
        assert "negative" in meshes

        # Check positive lobe
        pos_mesh = meshes["positive"]
        assert pos_mesh is not None
        assert pos_mesh.vertex_count > 0
        assert pos_mesh.triangle_count > 0
        assert len(pos_mesh.vertices_b64) > 0
        assert len(pos_mesh.normals_b64) > 0
        assert len(pos_mesh.indices_b64) > 0

        # Check negative lobe
        neg_mesh = meshes["negative"]
        assert neg_mesh is not None
        assert neg_mesh.vertex_count > 0

    def test_generate_density_mesh(self, simple_cube_data: CubeData) -> None:
        """Test generating density mesh."""
        mesh = generate_density_mesh(simple_cube_data, isovalue=0.5)

        assert mesh is not None
        assert mesh.vertex_count > 0
        assert mesh.triangle_count > 0

    def test_decode_mesh_array_vertices(self, orbital_cube_data: CubeData) -> None:
        """Test decoding compressed vertex array."""
        meshes = generate_orbital_meshes(orbital_cube_data, isovalue=0.1)
        pos_mesh = meshes["positive"]
        assert pos_mesh is not None

        vertices = decode_mesh_array(pos_mesh.vertices_b64, "float32")

        # Should have 3 coordinates per vertex
        assert len(vertices) == pos_mesh.vertex_count * 3
        assert vertices.dtype == np.float32

    def test_decode_mesh_array_indices(self, orbital_cube_data: CubeData) -> None:
        """Test decoding compressed index array."""
        meshes = generate_orbital_meshes(orbital_cube_data, isovalue=0.1)
        pos_mesh = meshes["positive"]
        assert pos_mesh is not None

        indices = decode_mesh_array(pos_mesh.indices_b64, "uint32")

        # Should have 3 indices per triangle
        assert len(indices) == pos_mesh.triangle_count * 3
        assert indices.dtype == np.uint32

    def test_compression_deterministic(self, orbital_cube_data: CubeData) -> None:
        """Test that compression is deterministic (mtime=0)."""
        meshes1 = generate_orbital_meshes(orbital_cube_data, isovalue=0.1)
        meshes2 = generate_orbital_meshes(orbital_cube_data, isovalue=0.1)

        pos1 = meshes1["positive"]
        pos2 = meshes2["positive"]

        assert pos1 is not None and pos2 is not None
        assert pos1.vertices_b64 == pos2.vertices_b64
        assert pos1.normals_b64 == pos2.normals_b64
        assert pos1.indices_b64 == pos2.indices_b64

    def test_no_surface_returns_none(self) -> None:
        """Test that no isosurface returns None."""
        # Create cube data with all values below isovalue
        cube = CubeData(
            origin=np.array([0, 0, 0]),
            axes=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            n_points=(5, 5, 5),
            data=np.zeros((5, 5, 5)),  # All zeros
            atom_positions=np.zeros((0, 3)),
            atom_numbers=[],
        )

        mesh = generate_density_mesh(cube, isovalue=0.5)
        assert mesh is None


class TestMeshData:
    """Tests for MeshData dataclass."""

    def test_mesh_data_fields(self) -> None:
        """Test MeshData structure."""
        mesh = MeshData(
            vertices_b64="abc",
            normals_b64="def",
            indices_b64="ghi",
            vertex_count=10,
            triangle_count=5,
        )

        assert mesh.vertices_b64 == "abc"
        assert mesh.normals_b64 == "def"
        assert mesh.indices_b64 == "ghi"
        assert mesh.vertex_count == 10
        assert mesh.triangle_count == 5
