"""Marching cubes mesh generation for molecular orbitals.

Converts volumetric data to Three.js-compatible meshes.
"""

import base64
import gzip
from dataclasses import dataclass
from io import BytesIO

import numpy as np
from numpy.typing import NDArray
from skimage import measure

from moleculens.compute.cube_parser import BOHR_TO_ANGSTROM, CubeData


@dataclass
class MeshData:
    """Mesh data in Three.js-compatible format."""

    # Compressed base64-encoded typed arrays
    vertices_b64: str  # Float32Array, gzipped
    normals_b64: str  # Float32Array, gzipped
    indices_b64: str  # Uint32Array, gzipped

    # Metadata
    vertex_count: int
    triangle_count: int


def _compress_array(arr: NDArray[np.floating | np.integer]) -> str:
    """Compress array to gzipped base64 string.

    Uses deterministic gzip (mtime=0) for reproducible results.
    """
    # Ensure correct dtype for Three.js
    if arr.dtype in (np.float64, np.float32):
        arr = arr.astype(np.float32)
    elif arr.dtype in (np.int64, np.int32, np.uint64, np.uint32):
        arr = arr.astype(np.uint32)

    # Flatten and get bytes
    data = arr.flatten().tobytes()

    # Gzip with mtime=0 for determinism
    buffer = BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as gz:
        gz.write(data)

    return base64.b64encode(buffer.getvalue()).decode("ascii")


def _compute_normals(
    vertices: NDArray[np.float64], faces: NDArray[np.int64]
) -> NDArray[np.float64]:
    """Compute vertex normals by averaging face normals.

    Args:
        vertices: Vertex positions, shape (n_vertices, 3)
        faces: Triangle indices, shape (n_faces, 3)

    Returns:
        Vertex normals, shape (n_vertices, 3), normalized
    """
    # Initialize normals
    normals = np.zeros_like(vertices)

    # Get vertices for each face
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]

    # Compute face normals
    edge1 = v1 - v0
    edge2 = v2 - v0
    face_normals = np.cross(edge1, edge2)

    # Accumulate normals for each vertex
    for i in range(3):
        np.add.at(normals, faces[:, i], face_normals)

    # Normalize
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    lengths = np.where(lengths > 0, lengths, 1)  # Avoid division by zero
    normals = normals / lengths

    return normals


def _extract_isosurface(
    cube: CubeData, isovalue: float, invert: bool = False
) -> tuple[NDArray[np.float64], NDArray[np.int64]] | None:
    """Extract isosurface using marching cubes.

    Args:
        cube: Cube data with volumetric field
        isovalue: Isosurface value (positive)
        invert: If True, negate the field (for negative lobe)

    Returns:
        Tuple of (vertices, faces) in Angstrom, or None if no surface found
    """
    data = cube.data
    if invert:
        data = -data

    # Check if isosurface exists in data
    if data.min() >= isovalue or data.max() <= isovalue:
        return None

    try:
        # Run marching cubes
        # spacing is the step size along each axis
        spacing = (
            np.linalg.norm(cube.axes[0]),
            np.linalg.norm(cube.axes[1]),
            np.linalg.norm(cube.axes[2]),
        )

        verts, faces, _, _ = measure.marching_cubes(
            data,
            level=isovalue,
            spacing=spacing,
            allow_degenerate=False,
        )
    except (ValueError, RuntimeError):
        return None

    if len(verts) == 0:
        return None

    # Transform vertices to physical coordinates
    # marching_cubes gives vertices in grid spacing units from origin (0,0,0)
    # We need to add the cube origin and handle axis orientation
    verts_physical = np.zeros_like(verts)

    for i in range(len(verts)):
        # Get grid position (in units of spacing, but we used actual spacing)
        # So verts are already in Bohr from grid origin at (0,0,0)
        # We need to add cube.origin and account for non-orthogonal grids

        # For orthogonal grids (which Psi4 produces), this simplifies:
        grid_pos = verts[i] / spacing  # Back to grid indices

        # Transform to physical coordinates
        pos_bohr = (
            cube.origin
            + grid_pos[0] * cube.axes[0]
            + grid_pos[1] * cube.axes[1]
            + grid_pos[2] * cube.axes[2]
        )
        verts_physical[i] = pos_bohr

    # Convert Bohr to Angstrom
    verts_angstrom = verts_physical * BOHR_TO_ANGSTROM

    return verts_angstrom, faces


def generate_orbital_meshes(cube: CubeData, isovalue: float = 0.05) -> dict[str, MeshData | None]:
    """Generate positive and negative lobe meshes for an orbital.

    Args:
        cube: Cube data with orbital volumetric data
        isovalue: Isosurface value for mesh generation

    Returns:
        Dict with 'positive' and 'negative' MeshData (or None if no surface)
    """
    result = {}

    # Positive lobe
    pos_surface = _extract_isosurface(cube, isovalue, invert=False)
    if pos_surface is not None:
        verts, faces = pos_surface
        normals = _compute_normals(verts, faces)
        result["positive"] = MeshData(
            vertices_b64=_compress_array(verts),
            normals_b64=_compress_array(normals),
            indices_b64=_compress_array(faces.flatten()),
            vertex_count=len(verts),
            triangle_count=len(faces),
        )
    else:
        result["positive"] = None

    # Negative lobe (negate field, then extract positive isovalue)
    neg_surface = _extract_isosurface(cube, isovalue, invert=True)
    if neg_surface is not None:
        verts, faces = neg_surface
        normals = _compute_normals(verts, faces)
        result["negative"] = MeshData(
            vertices_b64=_compress_array(verts),
            normals_b64=_compress_array(normals),
            indices_b64=_compress_array(faces.flatten()),
            vertex_count=len(verts),
            triangle_count=len(faces),
        )
    else:
        result["negative"] = None

    return result


def generate_density_mesh(cube: CubeData, isovalue: float = 0.02) -> MeshData | None:
    """Generate mesh for electron density isosurface.

    Args:
        cube: Cube data with density volumetric data
        isovalue: Isosurface value for mesh generation

    Returns:
        MeshData or None if no surface found
    """
    surface = _extract_isosurface(cube, isovalue, invert=False)
    if surface is None:
        return None

    verts, faces = surface
    normals = _compute_normals(verts, faces)

    return MeshData(
        vertices_b64=_compress_array(verts),
        normals_b64=_compress_array(normals),
        indices_b64=_compress_array(faces.flatten()),
        vertex_count=len(verts),
        triangle_count=len(faces),
    )


def decode_mesh_array(b64_data: str, dtype: str = "float32") -> NDArray:
    """Decode a compressed mesh array (for testing/validation).

    Args:
        b64_data: Base64-encoded gzipped array
        dtype: Target dtype ('float32' or 'uint32')

    Returns:
        Decoded numpy array
    """
    compressed = base64.b64decode(b64_data)
    decompressed = gzip.decompress(compressed)

    np_dtype = np.float32 if dtype == "float32" else np.uint32
    return np.frombuffer(decompressed, dtype=np_dtype)
