"""Gaussian cube file parser.

Parses .cube files produced by Psi4's cubeprop for orbital and density data.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# Conversion factor from Bohr to Angstrom
BOHR_TO_ANGSTROM = 0.529177249


@dataclass
class CubeData:
    """Parsed Gaussian cube file data."""

    # Grid origin in Bohr
    origin: NDArray[np.float64]  # Shape: (3,)

    # Grid axis vectors (each row is an axis, columns are step sizes in Bohr)
    # Shape: (3, 3) - axes[i] gives direction vector for axis i
    axes: NDArray[np.float64]

    # Number of points along each axis
    n_points: tuple[int, int, int]

    # Volumetric data (scalar field values)
    data: NDArray[np.float64]  # Shape: (nx, ny, nz)

    # Atom positions in Bohr (for reference)
    atom_positions: NDArray[np.float64]  # Shape: (n_atoms, 3)
    atom_numbers: list[int]

    @property
    def origin_angstrom(self) -> NDArray[np.float64]:
        """Grid origin in Angstrom."""
        return self.origin * BOHR_TO_ANGSTROM

    @property
    def axes_angstrom(self) -> NDArray[np.float64]:
        """Grid axis vectors in Angstrom."""
        return self.axes * BOHR_TO_ANGSTROM

    def grid_to_angstrom(self, i: int, j: int, k: int) -> NDArray[np.float64]:
        """Convert grid indices to Angstrom coordinates."""
        pos_bohr = self.origin + i * self.axes[0] + j * self.axes[1] + k * self.axes[2]
        return pos_bohr * BOHR_TO_ANGSTROM


class CubeParseError(Exception):
    """Raised when cube file parsing fails."""


def parse_cube_file(path: Path | str) -> CubeData:
    """Parse a Gaussian cube file.

    Args:
        path: Path to the .cube file.

    Returns:
        CubeData with parsed volumetric data.

    Raises:
        CubeParseError: If the file format is invalid.
    """
    path = Path(path)

    if not path.exists():
        raise CubeParseError(f"Cube file not found: {path}")

    with path.open("r") as f:
        lines = f.readlines()

    if len(lines) < 6:
        raise CubeParseError("Cube file too short")

    try:
        # Lines 1-2: Comments (ignore)
        # Line 3: Number of atoms and origin
        parts = lines[2].split()
        n_atoms = int(parts[0])
        # Negative n_atoms means orbital data (we handle both)
        n_atoms = abs(n_atoms)
        origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])

        # Lines 4-6: Number of voxels and axis vectors
        axes = np.zeros((3, 3))
        n_points = []

        for i in range(3):
            parts = lines[3 + i].split()
            n_points.append(int(parts[0]))
            axes[i] = [float(parts[1]), float(parts[2]), float(parts[3])]

        n_points = tuple(n_points)

        # Lines 7 to 6+n_atoms: Atom data
        atom_numbers = []
        atom_positions = []

        for i in range(n_atoms):
            parts = lines[6 + i].split()
            atom_numbers.append(int(parts[0]))
            # Parts: atomic_number, charge, x, y, z
            atom_positions.append([float(parts[2]), float(parts[3]), float(parts[4])])

        atom_positions = np.array(atom_positions) if atom_positions else np.zeros((0, 3))

        # Remaining lines: Volumetric data
        # Data is stored in row-major order: for ix, for iy, for iz
        data_start = 6 + n_atoms

        # Some cube files have an extra line for MO coefficients after atoms
        # Check if first data line looks like orbital info
        if data_start < len(lines):
            first_data = lines[data_start].split()
            # Orbital info line typically has fewer numbers or starts with integer
            if len(first_data) == 2 and first_data[0].isdigit():
                data_start += 1

        # Parse volumetric data
        all_values: list[float] = []
        for line in lines[data_start:]:
            values = line.split()
            all_values.extend(float(v) for v in values)

        expected_size = n_points[0] * n_points[1] * n_points[2]
        if len(all_values) < expected_size:
            raise CubeParseError(
                f"Insufficient data: expected {expected_size} values, got {len(all_values)}"
            )

        # Reshape to 3D array
        data = np.array(all_values[:expected_size]).reshape(n_points)

    except (ValueError, IndexError) as e:
        raise CubeParseError(f"Failed to parse cube file: {e}") from e

    return CubeData(
        origin=origin,
        axes=axes,
        n_points=n_points,  # type: ignore[arg-type]
        data=data,
        atom_positions=atom_positions,
        atom_numbers=atom_numbers,
    )
