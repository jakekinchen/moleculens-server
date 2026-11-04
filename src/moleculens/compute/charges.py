"""Atomic charge and dipole moment calculation.

Supports:
- xTB (GFN2-xTB) via xtb-python or subprocess
- RDKit Gasteiger charges when explicitly requested
"""

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from moleculens.compute.sdf_parser import ELEMENTS, SDFMolecule
from moleculens.core import get_logger

logger = get_logger(__name__)

# Conversion: Debye = e·Å × 4.80320471257
DEBYE_PER_E_ANGSTROM = 4.80320471257


@dataclass
class ChargeResult:
    """Result of charge calculation."""

    charges: NDArray[np.float64]  # Atomic partial charges (e)
    dipole_vector: NDArray[np.float64]  # Dipole moment vector (Debye)
    dipole_magnitude: float  # Dipole magnitude (Debye)
    dipole_origin: NDArray[np.float64]  # Center of mass / geometry
    charge_model: str  # e.g., "xtb-gfn2", "rdkit-gasteiger"
    convention: str  # "physics" (negative to positive) or "chemistry" (opposite)


def calculate_charges_xtb(
    molecule: SDFMolecule,
    charge: int = 0,
    multiplicity: int = 1,
) -> ChargeResult | None:
    """Calculate charges using xTB (GFN2-xTB).

    Args:
        molecule: Parsed SDF molecule
        charge: Total molecular charge
        multiplicity: Spin multiplicity

    Returns:
        ChargeResult or None if xTB not available/fails
    """
    try:
        # Try xtb-python first
        return _calculate_xtb_python(molecule, charge, multiplicity)
    except ImportError:
        pass

    # Fall back to subprocess
    try:
        return _calculate_xtb_subprocess(molecule, charge, multiplicity)
    except (FileNotFoundError, subprocess.SubprocessError) as e:
        logger.warning("xTB calculation failed", error=str(e))
        return None


def _calculate_xtb_python(
    molecule: SDFMolecule,
    charge: int,
    multiplicity: int,
) -> ChargeResult:
    """Calculate charges using xtb-python bindings."""
    from xtb.interface import Calculator, Param
    from xtb.libxtb import VERBOSITY_MUTED

    # Build atomic numbers and positions
    numbers = np.array([_element_to_number(a.symbol) for a in molecule.atoms])
    positions = np.array([[a.x, a.y, a.z] for a in molecule.atoms])

    # Convert positions from Angstrom to Bohr
    BOHR_PER_ANGSTROM = 1.8897259886
    positions_bohr = positions * BOHR_PER_ANGSTROM

    # Set up calculator
    calc = Calculator(Param.GFN2xTB, numbers, positions_bohr, charge=charge, uhf=multiplicity - 1)
    calc.set_verbosity(VERBOSITY_MUTED)

    # Calculate
    res = calc.singlepoint()

    # Get charges (Mulliken by default in xtb-python)
    charges = res.get_charges()

    # Get dipole moment (in e·Bohr, convert to Debye)
    dipole_bohr = res.get_dipole()
    dipole_angstrom = dipole_bohr / BOHR_PER_ANGSTROM
    dipole_debye = dipole_angstrom * DEBYE_PER_E_ANGSTROM

    # Translate the engine dipole from coordinate zero to the declared origin.
    origin = positions.mean(axis=0)
    dipole_debye -= charge * origin * DEBYE_PER_E_ANGSTROM

    return ChargeResult(
        charges=charges,
        dipole_vector=dipole_debye,
        dipole_magnitude=float(np.linalg.norm(dipole_debye)),
        dipole_origin=origin,
        charge_model="xtb-gfn2-mulliken",
        convention="physics",  # Electric dipole vector points negative to positive
    )


def _calculate_xtb_subprocess(
    molecule: SDFMolecule,
    charge: int,
    multiplicity: int,
) -> ChargeResult:
    """Calculate charges using xTB command line tool."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        # Write XYZ file
        xyz_file = tmppath / "molecule.xyz"
        with xyz_file.open("w") as f:
            f.write(f"{molecule.num_atoms}\n")
            f.write(f"charge={charge}\n")
            for atom in molecule.atoms:
                f.write(f"{atom.symbol} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}\n")

        # Run xTB
        uhf = multiplicity - 1
        cmd = [
            "xtb",
            str(xyz_file),
            "--gfn2",
            "--chrg",
            str(charge),
            "--uhf",
            str(uhf),
            "--json",
        ]

        result = subprocess.run(
            cmd,
            cwd=tmppath,
            capture_output=True,
            text=True,
            timeout=120,
        )

        if result.returncode != 0:
            raise subprocess.SubprocessError(f"xTB failed: {result.stderr}")

        # Parse xtbout.json
        import json

        json_file = tmppath / "xtbout.json"
        if not json_file.exists():
            raise FileNotFoundError("xTB did not produce xtbout.json")

        with json_file.open() as f:
            data = json.load(f)

        # xTB 6.7 JSON reports the vector in atomic units, not Debye.
        # Unknown/missing schemas are errors; never synthesize a zero dipole.
        charges = np.asarray(data["partial charges"], dtype=float)
        dipole_au = np.asarray(data["dipole / a.u."], dtype=float)
        if charges.shape != (molecule.num_atoms,) or dipole_au.shape != (3,):
            raise ValueError("Unsupported xTB result shape")
        if not np.isfinite(charges).all() or not np.isfinite(dipole_au).all():
            raise ValueError("Non-finite xTB result")
        positions = np.array([[a.x, a.y, a.z] for a in molecule.atoms])
        origin = positions.mean(axis=0)
        dipole_debye = dipole_au * (DEBYE_PER_E_ANGSTROM / 1.8897259886)
        dipole_debye -= charge * origin * DEBYE_PER_E_ANGSTROM
        dipole_mag = float(np.linalg.norm(dipole_debye))

        return ChargeResult(
            charges=charges,
            dipole_vector=dipole_debye,
            dipole_magnitude=dipole_mag,
            dipole_origin=origin,
            charge_model="xtb-gfn2",
            convention="physics",
        )


def calculate_charges_gasteiger(
    molecule: SDFMolecule, sdf_content: str | None = None
) -> ChargeResult:
    """Calculate Gasteiger charges using RDKit.

    This is an explicitly requested approximation, never an implicit xTB fallback.

    Args:
        molecule: Parsed SDF molecule
        sdf_content: Original SDF content for RDKit parsing (optional)

    Returns:
        ChargeResult with Gasteiger charges and computed dipole
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ImportError("RDKit required for Gasteiger charges") from e

    positions = np.array([[a.x, a.y, a.z] for a in molecule.atoms])

    # Try parsing the SDF directly if content is provided
    mol = None
    if sdf_content:
        try:
            mol = Chem.MolFromMolBlock(sdf_content, sanitize=True, removeHs=False)
        except Exception:
            mol = None

    if mol is None or mol.GetNumAtoms() != molecule.num_atoms:
        raise ValueError("Gasteiger requires a valid molecular graph with matching atom order")
    AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=True)
    # Condense implicit-H charges onto their parent positions. This convention is
    # explicit in charge_model; total formal charge remains conserved.
    charges = np.array(
        [
            a.GetDoubleProp("_GasteigerCharge") + a.GetDoubleProp("_GasteigerHCharge")
            for a in mol.GetAtoms()
        ]
    )
    if not np.isfinite(charges).all():
        raise ValueError("Gasteiger produced non-finite charges")
    if not np.isclose(charges.sum(), Chem.GetFormalCharge(mol), atol=1e-5):
        raise ValueError("Gasteiger charges do not conserve formal charge")
    origin = positions.mean(axis=0)
    dipole_debye = np.sum(charges[:, None] * (positions - origin), axis=0) * DEBYE_PER_E_ANGSTROM

    return ChargeResult(
        charges=charges,
        dipole_vector=dipole_debye,
        dipole_magnitude=float(np.linalg.norm(dipole_debye)),
        dipole_origin=origin,
        charge_model="rdkit-gasteiger-condensed-h",
        convention="physics",
    )


def calculate_coulombic_potential(
    vertex_positions: NDArray[np.float64],
    atom_positions: NDArray[np.float64],
    charges: NDArray[np.float64],
    softening_epsilon: float = 0.3,
) -> NDArray[np.float64]:
    """Calculate Coulombic potential at surface vertices.

    V(r) = 332.06371 * Σ (qᵢ / sqrt(|r-rᵢ|² + ε²))

    where 332.06371 is the conversion factor for kcal/mol/e when
    distances are in Angstrom and charges in elementary charge units.

    Args:
        vertex_positions: Surface vertex positions, shape (n_verts, 3)
        atom_positions: Atom positions, shape (n_atoms, 3)
        charges: Atomic partial charges, shape (n_atoms,)
        softening_epsilon: Softening parameter to avoid singularities (Å)

    Returns:
        Potential values at each vertex in kcal/mol/e, shape (n_verts,)
    """
    # Coulomb constant in kcal/mol/e * Å
    COULOMB_CONST = 332.06371

    n_verts = len(vertex_positions)
    n_atoms = len(atom_positions)

    # Compute pairwise distances with softening
    # Shape: (n_verts, n_atoms)
    potential = np.zeros(n_verts)

    eps_sq = softening_epsilon**2

    for i in range(n_atoms):
        # Distance from each vertex to this atom
        diff = vertex_positions - atom_positions[i]
        dist_sq = np.sum(diff**2, axis=1) + eps_sq
        dist = np.sqrt(dist_sq)

        # Add contribution from this atom
        potential += COULOMB_CONST * charges[i] / dist

    return potential


def _element_to_number(symbol: str) -> int:
    """Convert element symbol to atomic number."""
    try:
        return ELEMENTS.index(symbol) + 1
    except ValueError as exc:
        raise ValueError(f"Unknown element: {symbol}") from exc
