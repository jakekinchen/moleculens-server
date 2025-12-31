"""Psi4 computation runner for orbital and density calculations.

This module interfaces with Psi4 to compute molecular orbitals and electron density,
then generates meshes from the resulting cube files.
"""

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from moleculens.compute.cube_parser import parse_cube_file
from moleculens.compute.marching_cubes import (
    generate_density_mesh,
    generate_orbital_meshes,
)
from moleculens.compute.sdf_parser import parse_sdf
from moleculens.core import get_logger, settings

logger = get_logger(__name__)

# Conversion factor from Angstrom to Bohr
ANGSTROM_TO_BOHR = 1.8897259886


@dataclass
class ComputationResult:
    """Result of a Psi4 orbital computation."""

    orbitals: dict[str, dict[str, Any]] = field(default_factory=dict)
    density: dict[str, Any] | None = None
    meta: dict[str, Any] = field(default_factory=dict)
    artifact_files: list[str] = field(default_factory=list)
    compute_time_ms: float = 0.0


class Psi4ComputationError(Exception):
    """Raised when Psi4 computation fails."""


def run_psi4_computation(
    sdf_content: str,
    method: str,
    basis: str,
    grid_spacing: float,
    isovalue: float,
    orbitals: list[str],
    output_dir: Path,
    scratch_dir: Path | None = None,
) -> ComputationResult:
    """Run Psi4 computation for orbital/density data.

    Args:
        sdf_content: SDF file content with 3D coordinates
        method: Computational method (scf, hf, b3lyp, etc.)
        basis: Basis set name
        grid_spacing: Grid spacing in Angstrom
        isovalue: Isosurface value for mesh generation
        orbitals: List of orbital types to compute
        output_dir: Directory to store output files
        scratch_dir: Psi4 scratch directory

    Returns:
        ComputationResult with mesh data and metadata

    Raises:
        Psi4ComputationError: If computation fails
    """
    start_time = time.time()

    # Parse SDF
    try:
        molecule = parse_sdf(sdf_content)
    except Exception as e:
        raise Psi4ComputationError(f"Failed to parse SDF: {e}") from e

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set up scratch directory
    if scratch_dir is None:
        scratch_dir = settings.psi_scratch
    scratch_dir.mkdir(parents=True, exist_ok=True)

    # Import psi4 (only available in worker environment)
    try:
        import psi4
    except ImportError as e:
        raise Psi4ComputationError("Psi4 not available in this environment") from e

    # Configure Psi4
    psi4.set_memory(f"{settings.psi4_memory_mb} MB")
    psi4.set_num_threads(settings.psi4_num_threads)
    psi4.core.set_output_file(str(output_dir / "psi4_output.dat"), False)

    # Set scratch directory
    psi4_io = psi4.core.IOManager.shared_object()
    psi4_io.set_default_path(str(scratch_dir))

    # Create Psi4 geometry
    geom_str = molecule.to_psi4_geometry()
    try:
        mol = psi4.geometry(geom_str)
    except Exception as e:
        raise Psi4ComputationError(f"Failed to create Psi4 geometry: {e}") from e

    logger.info(
        "Starting Psi4 computation",
        method=method,
        basis=basis,
        num_atoms=molecule.num_atoms,
    )

    # Run energy calculation to get wavefunction
    try:
        # Handle different method specifications
        method_spec = method.lower()
        if method_spec in ("scf", "hf"):
            energy, wfn = psi4.energy("scf", molecule=mol, basis=basis, return_wfn=True)
        else:
            energy, wfn = psi4.energy(method_spec, molecule=mol, basis=basis, return_wfn=True)
    except Exception as e:
        raise Psi4ComputationError(f"Psi4 energy calculation failed: {e}") from e

    logger.info("Energy calculation complete", energy_hartree=energy)

    # Determine orbital indices
    # Closed shell assumption: HOMO = nalpha, LUMO = nalpha + 1
    nalpha = wfn.nalpha()
    homo_idx = nalpha  # 1-indexed in Psi4
    lumo_idx = nalpha + 1

    # Get orbital energies
    epsilon = wfn.epsilon_a()
    homo_energy = float(epsilon.get(homo_idx - 1))  # 0-indexed for epsilon
    lumo_energy = float(epsilon.get(lumo_idx - 1)) if lumo_idx <= epsilon.dim() else None

    # Convert to eV (1 Hartree = 27.211386 eV)
    HARTREE_TO_EV = 27.211386
    homo_energy_ev = homo_energy * HARTREE_TO_EV
    lumo_energy_ev = lumo_energy * HARTREE_TO_EV if lumo_energy else None

    # Configure cubeprop
    # Grid spacing conversion: Angstrom to Bohr
    grid_spacing_bohr = grid_spacing * ANGSTROM_TO_BOHR

    # Build list of orbitals to compute
    orbital_indices = []
    requested_orbitals = set(orbitals)

    if "homo" in requested_orbitals or "homo-1" in requested_orbitals:
        orbital_indices.append(homo_idx)
    if "lumo" in requested_orbitals or "lumo+1" in requested_orbitals:
        orbital_indices.append(lumo_idx)
    if "homo-1" in requested_orbitals and homo_idx > 1:
        orbital_indices.append(homo_idx - 1)
    if "lumo+1" in requested_orbitals and lumo_idx < epsilon.dim():
        orbital_indices.append(lumo_idx + 1)

    # Remove duplicates and sort
    orbital_indices = sorted(set(orbital_indices))

    # Set cubeprop options
    psi4.set_options(
        {
            "CUBEPROP_TASKS": ["ORBITALS"]
            + (["DENSITY"] if "density" in requested_orbitals else []),
            "CUBEPROP_ORBITALS": orbital_indices,
            "CUBIC_GRID_SPACING": [grid_spacing_bohr] * 3,
            "CUBIC_GRID_OVERAGE": [4.0, 4.0, 4.0],  # Bohr overage
            "CUBEPROP_FILEPATH": str(output_dir),
        }
    )

    # Run cubeprop
    try:
        psi4.cubeprop(wfn)
    except Exception as e:
        raise Psi4ComputationError(f"Cubeprop failed: {e}") from e

    # Process cube files and generate meshes
    result = ComputationResult()
    result.meta = {
        "method": method,
        "basis": basis,
        "gridSpacingAngstrom": grid_spacing,
        "psi4NumThreads": settings.psi4_num_threads,
        "totalEnergy": energy,
    }

    # Map orbital indices to names
    idx_to_name = {}
    if homo_idx in orbital_indices:
        idx_to_name[homo_idx] = "homo"
    if lumo_idx in orbital_indices:
        idx_to_name[lumo_idx] = "lumo"
    if (homo_idx - 1) in orbital_indices:
        idx_to_name[homo_idx - 1] = "homo-1"
    if (lumo_idx + 1) in orbital_indices:
        idx_to_name[lumo_idx + 1] = "lumo+1"

    # Process each orbital cube file
    for idx in orbital_indices:
        orb_name = idx_to_name.get(idx, f"mo_{idx}")

        # Psi4 names cube files as Psi_<mol>_<idx>_<type>.cube
        # Find the cube file for this orbital
        cube_pattern = f"Psi_*_{idx}_*.cube"
        cube_files = list(output_dir.glob(cube_pattern))

        if not cube_files:
            # Try alternative naming
            cube_files = list(output_dir.glob(f"*{idx}*.cube"))

        if not cube_files:
            logger.warning("Cube file not found for orbital", orbital_idx=idx)
            continue

        cube_file = cube_files[0]
        logger.info("Processing orbital cube file", path=str(cube_file))

        try:
            cube_data = parse_cube_file(cube_file)
            meshes = generate_orbital_meshes(cube_data, isovalue)

            orbital_data: dict[str, Any] = {"isovalue": isovalue}

            # Get orbital energy
            if orb_name == "homo":
                orbital_data["energyEv"] = homo_energy_ev
            elif orb_name == "lumo" and lumo_energy_ev is not None:
                orbital_data["energyEv"] = lumo_energy_ev

            # Add mesh data
            if meshes.get("positive"):
                orbital_data["positive"] = {
                    "vertices": meshes["positive"].vertices_b64,
                    "normals": meshes["positive"].normals_b64,
                    "indices": meshes["positive"].indices_b64,
                    "vertexCount": meshes["positive"].vertex_count,
                    "triangleCount": meshes["positive"].triangle_count,
                }
            if meshes.get("negative"):
                orbital_data["negative"] = {
                    "vertices": meshes["negative"].vertices_b64,
                    "normals": meshes["negative"].normals_b64,
                    "indices": meshes["negative"].indices_b64,
                    "vertexCount": meshes["negative"].vertex_count,
                    "triangleCount": meshes["negative"].triangle_count,
                }

            result.orbitals[orb_name] = orbital_data
            result.artifact_files.append(cube_file.name)

        except Exception as e:
            logger.error("Failed to process orbital", orbital=orb_name, error=str(e))

    # Process density cube file
    if "density" in requested_orbitals:
        density_files = list(output_dir.glob("*Dt.cube")) + list(output_dir.glob("*density*.cube"))
        if density_files:
            try:
                cube_data = parse_cube_file(density_files[0])
                density_mesh = generate_density_mesh(cube_data, isovalue)

                if density_mesh:
                    result.density = {
                        "vertices": density_mesh.vertices_b64,
                        "normals": density_mesh.normals_b64,
                        "indices": density_mesh.indices_b64,
                        "vertexCount": density_mesh.vertex_count,
                        "triangleCount": density_mesh.triangle_count,
                        "isovalue": isovalue,
                    }
                    result.artifact_files.append(density_files[0].name)
            except Exception as e:
                logger.error("Failed to process density", error=str(e))

    # Calculate computation time
    result.compute_time_ms = (time.time() - start_time) * 1000
    result.meta["computeTimeMs"] = result.compute_time_ms

    # Save result.json
    result_data = {
        "orbitals": result.orbitals,
        "density": result.density,
        "meta": result.meta,
    }
    with (output_dir / "result.json").open("w") as f:
        json.dump(result_data, f, indent=2)
    result.artifact_files.append("result.json")

    logger.info(
        "Computation complete",
        compute_time_ms=result.compute_time_ms,
        num_orbitals=len(result.orbitals),
        has_density=result.density is not None,
    )

    # Clean up Psi4 scratch files
    psi4.core.clean()

    return result
