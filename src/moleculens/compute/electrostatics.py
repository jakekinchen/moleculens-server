"""Electrostatics computation runner.

Computes molecular surface colored by electrostatic potential (ESP/MEP)
and dipole moment.
"""

import hashlib
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from moleculens.compute.charges import (
    calculate_charges_gasteiger,
    calculate_charges_xtb,
    calculate_coulombic_potential,
)
from moleculens.compute.marching_cubes import _compress_array
from moleculens.compute.sdf_parser import parse_sdf
from moleculens.compute.surface_mesh import generate_molecular_surface
from moleculens.core import get_logger

logger = get_logger(__name__)


@dataclass
class ElectrostaticsResult:
    """Result of electrostatics computation."""

    # Surface mesh with potential
    surface: dict[str, Any] = field(default_factory=dict)

    # Dipole information
    dipole: dict[str, Any] = field(default_factory=dict)

    # Metadata
    meta: dict[str, Any] = field(default_factory=dict)

    # Artifact files created
    artifact_files: list[str] = field(default_factory=list)

    # Timing
    compute_time_ms: float = 0.0


class ElectrostaticsComputationError(Exception):
    """Raised when electrostatics computation fails."""


def compute_sdf_hash(sdf_content: str) -> str:
    """Compute stable hash of SDF content."""
    return hashlib.sha256(sdf_content.strip().encode()).hexdigest()[:16]


def run_electrostatics_computation(
    sdf_content: str,
    charge: int = 0,
    multiplicity: int = 1,
    method: str = "xtb-gfn2",
    grid_spacing: float = 0.25,
    probe_radius: float = 1.4,
    padding: float = 3.0,
    softening_epsilon: float = 0.3,
    clamp_percentiles: tuple[float, float] = (5.0, 95.0),
    output_dir: Path | None = None,
) -> ElectrostaticsResult:
    """Run electrostatics computation for ESP surface and dipole.

    Args:
        sdf_content: SDF file content with 3D coordinates
        charge: Total molecular charge
        multiplicity: Spin multiplicity
        method: Charge method ("xtb-gfn2" or "gasteiger")
        grid_spacing: Surface grid spacing in Angstrom
        probe_radius: Probe radius for surface in Angstrom
        padding: Padding around molecule in Angstrom
        softening_epsilon: Softening parameter for potential calculation
        clamp_percentiles: Percentiles for potential clamping (p_low, p_high)
        output_dir: Directory to store output files

    Returns:
        ElectrostaticsResult with surface mesh, potential, and dipole

    Raises:
        ElectrostaticsComputationError: If computation fails
    """
    start_time = time.time()
    timings: dict[str, float] = {}

    # Parse SDF
    try:
        molecule = parse_sdf(sdf_content)
    except Exception as e:
        raise ElectrostaticsComputationError(f"Failed to parse SDF: {e}") from e

    if molecule.num_atoms == 0:
        raise ElectrostaticsComputationError("No atoms in molecule")

    # Validate 3D coordinates
    positions = np.array([[a.x, a.y, a.z] for a in molecule.atoms])
    if np.allclose(positions[:, 2], 0.0):
        # All z-coordinates are zero - likely 2D
        logger.warning("Molecule appears to be 2D (all z=0)")

    atom_symbols = [a.symbol for a in molecule.atoms]

    logger.info(
        "Starting electrostatics computation",
        method=method,
        num_atoms=molecule.num_atoms,
    )

    # Step 1: Calculate charges and dipole
    t0 = time.time()
    charge_result = None

    if method.lower().startswith("xtb"):
        charge_result = calculate_charges_xtb(molecule, charge, multiplicity)

    if charge_result is None:
        # Fall back to Gasteiger
        logger.info("Using Gasteiger charges (xTB not available or failed)")
        try:
            charge_result = calculate_charges_gasteiger(molecule, sdf_content)
        except ImportError as e:
            raise ElectrostaticsComputationError(
                "Neither xTB nor RDKit available for charge calculation"
            ) from e

    timings["charges"] = (time.time() - t0) * 1000

    logger.info(
        "Charges calculated",
        model=charge_result.charge_model,
        dipole_magnitude=f"{charge_result.dipole_magnitude:.2f} D",
    )

    # Step 2: Generate molecular surface
    t0 = time.time()
    try:
        surface_mesh, vertex_positions = generate_molecular_surface(
            atom_positions=positions,
            atom_symbols=atom_symbols,
            grid_spacing=grid_spacing,
            probe_radius=probe_radius,
            padding=padding,
        )
    except Exception as e:
        raise ElectrostaticsComputationError(f"Surface generation failed: {e}") from e

    timings["surface"] = (time.time() - t0) * 1000

    logger.info(
        "Surface generated",
        vertices=surface_mesh.vertex_count,
        triangles=surface_mesh.triangle_count,
    )

    # Step 3: Calculate potential at surface vertices
    t0 = time.time()
    potential = calculate_coulombic_potential(
        vertex_positions=vertex_positions,
        atom_positions=positions,
        charges=charge_result.charges,
        softening_epsilon=softening_epsilon,
    )
    timings["potential"] = (time.time() - t0) * 1000

    # Compute statistics for clamping
    pot_min = float(np.min(potential))
    pot_max = float(np.max(potential))
    pot_p05 = float(np.percentile(potential, clamp_percentiles[0]))
    pot_p95 = float(np.percentile(potential, clamp_percentiles[1]))

    logger.info(
        "Potential calculated",
        min=f"{pot_min:.2f}",
        max=f"{pot_max:.2f}",
        p05=f"{pot_p05:.2f}",
        p95=f"{pot_p95:.2f}",
    )

    # Build result
    total_time = (time.time() - start_time) * 1000
    timings["total"] = total_time

    result = ElectrostaticsResult()
    result.compute_time_ms = total_time

    result.surface = {
        "mesh": {
            "vertices": surface_mesh.vertices_b64,
            "normals": surface_mesh.normals_b64,
            "indices": surface_mesh.indices_b64,
            "potential": _compress_array(potential.astype(np.float32)),
            "vertexCount": surface_mesh.vertex_count,
            "triangleCount": surface_mesh.triangle_count,
        },
        "range": {
            "min": pot_min,
            "max": pot_max,
            "p05": pot_p05,
            "p95": pot_p95,
        },
        "units": {
            "coords": "angstrom",
            "potential": "kcal/mol/e",
        },
    }

    result.dipole = {
        "origin_angstrom": charge_result.dipole_origin.tolist(),
        "vector_debye": charge_result.dipole_vector.tolist(),
        "magnitude_debye": charge_result.dipole_magnitude,
        "convention": charge_result.convention,
    }

    result.meta = {
        "charge_model": charge_result.charge_model,
        "method": method,
        "molecular_charge": charge,
        "multiplicity": multiplicity,
        "surface_params": {
            "grid_spacing": grid_spacing,
            "probe_radius": probe_radius,
            "padding": padding,
        },
        "potential_params": {
            "softening_epsilon": softening_epsilon,
            "clamp_percentiles": list(clamp_percentiles),
        },
        "timings_ms": timings,
        "sdf_hash": compute_sdf_hash(sdf_content),
    }

    # Save to output directory if provided
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

        result_data = {
            "surface": result.surface,
            "dipole": result.dipole,
            "meta": result.meta,
        }

        result_file = output_dir / "electrostatics_result.json"
        with result_file.open("w") as f:
            json.dump(result_data, f, indent=2)
        result.artifact_files.append("electrostatics_result.json")

    logger.info(
        "Electrostatics computation complete",
        total_time_ms=total_time,
    )

    return result
