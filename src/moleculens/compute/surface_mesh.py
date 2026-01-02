"""Molecular surface mesh generation using signed distance fields.

Generates a smooth molecular surface as a union of spheres (van der Waals radii)
with optional probe radius, then extracts the isosurface via marching cubes.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from skimage import measure

from moleculens.compute.marching_cubes import _compress_array, _compute_normals

# Van der Waals radii in Angstrom (common elements)
VDW_RADII: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
    # Metals and others
    "Li": 1.82,
    "Na": 2.27,
    "K": 2.75,
    "Mg": 1.73,
    "Ca": 2.31,
    "Fe": 2.04,
    "Zn": 1.39,
    "Cu": 1.40,
    "B": 1.92,
    "Si": 2.10,
    "Se": 1.90,
}

DEFAULT_VDW_RADIUS = 1.70  # Default for unknown elements


@dataclass
class MolecularSurfaceMesh:
    """Result of molecular surface mesh generation."""

    vertices_b64: str  # gzipped base64 Float32Array [x0,y0,z0,...]
    normals_b64: str  # gzipped base64 Float32Array
    indices_b64: str  # gzipped base64 Uint32Array
    vertex_count: int
    triangle_count: int


def get_vdw_radius(symbol: str) -> float:
    """Get van der Waals radius for an element."""
    return VDW_RADII.get(symbol.capitalize(), DEFAULT_VDW_RADIUS)


def generate_molecular_surface(
    atom_positions: NDArray[np.float64],
    atom_symbols: list[str],
    grid_spacing: float = 0.25,
    probe_radius: float = 1.4,
    padding: float = 3.0,
) -> tuple[MolecularSurfaceMesh, NDArray[np.float64]]:
    """Generate a molecular surface mesh using union of spheres SDF.

    Args:
        atom_positions: Atom coordinates in Angstrom, shape (n_atoms, 3)
        atom_symbols: Element symbols for each atom
        grid_spacing: Grid spacing in Angstrom
        probe_radius: Probe radius in Angstrom (typically 1.4 for water)
        padding: Padding around molecule in Angstrom

    Returns:
        Tuple of (MolecularSurfaceMesh, vertex_positions_array)
        vertex_positions_array is shape (n_vertices, 3) for potential calculation
    """
    n_atoms = len(atom_positions)
    if n_atoms == 0:
        raise ValueError("No atoms provided")

    # Get vdW radii for each atom
    radii = np.array([get_vdw_radius(s) + probe_radius for s in atom_symbols])

    # Compute bounding box
    min_coords = atom_positions.min(axis=0) - radii.max() - padding
    max_coords = atom_positions.max(axis=0) + radii.max() + padding

    # Create 3D grid
    nx = int(np.ceil((max_coords[0] - min_coords[0]) / grid_spacing)) + 1
    ny = int(np.ceil((max_coords[1] - min_coords[1]) / grid_spacing)) + 1
    nz = int(np.ceil((max_coords[2] - min_coords[2]) / grid_spacing)) + 1

    # Limit grid size for safety
    max_grid_size = 200
    if max(nx, ny, nz) > max_grid_size:
        scale = max_grid_size / max(nx, ny, nz)
        nx = int(nx * scale)
        ny = int(ny * scale)
        nz = int(nz * scale)
        grid_spacing = (max_coords[0] - min_coords[0]) / (nx - 1)

    # Create coordinate arrays
    x = np.linspace(min_coords[0], max_coords[0], nx)
    y = np.linspace(min_coords[1], max_coords[1], ny)
    z = np.linspace(min_coords[2], max_coords[2], nz)

    # Create meshgrid
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # Compute signed distance field (union of spheres)
    # SDF = min over all atoms of (distance to atom center - radius)
    sdf = np.full((nx, ny, nz), np.inf)

    for i in range(n_atoms):
        atom_pos = atom_positions[i]
        radius = radii[i]

        # Distance from each grid point to this atom
        dist = np.sqrt((X - atom_pos[0]) ** 2 + (Y - atom_pos[1]) ** 2 + (Z - atom_pos[2]) ** 2)

        # SDF for this sphere
        sphere_sdf = dist - radius

        # Union: take minimum
        sdf = np.minimum(sdf, sphere_sdf)

    # Extract isosurface at SDF = 0
    try:
        spacing = (grid_spacing, grid_spacing, grid_spacing)
        verts, faces, _, _ = measure.marching_cubes(
            sdf,
            level=0.0,
            spacing=spacing,
            allow_degenerate=False,
        )
    except (ValueError, RuntimeError) as e:
        raise ValueError(f"Marching cubes failed: {e}") from e

    if len(verts) == 0:
        raise ValueError("No surface generated - check atom positions")

    # Transform vertices to world coordinates
    # marching_cubes returns vertices in grid units * spacing from origin (0,0,0)
    # We need to add min_coords to get world coordinates
    verts_world = verts + min_coords

    # Compute normals
    normals = _compute_normals(verts_world, faces)

    # Compress arrays
    mesh = MolecularSurfaceMesh(
        vertices_b64=_compress_array(verts_world.flatten()),
        normals_b64=_compress_array(normals.flatten()),
        indices_b64=_compress_array(faces.flatten()),
        vertex_count=len(verts_world),
        triangle_count=len(faces),
    )

    return mesh, verts_world
