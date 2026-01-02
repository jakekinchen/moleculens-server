"""Atomic charge and dipole moment calculation.

Supports:
- xTB (GFN2-xTB) via xtb-python or subprocess
- RDKit Gasteiger charges as fallback
"""

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from moleculens.compute.sdf_parser import SDFMolecule
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
    convention: str  # "physics" (+ to -) or "chemistry" (- to +)


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
    calc = Calculator(Param.GFN2xTB, numbers, positions_bohr, charge=charge)
    calc.set_verbosity(VERBOSITY_MUTED)

    # Calculate
    res = calc.singlepoint()

    # Get charges (Mulliken by default in xtb-python)
    charges = res.get_charges()

    # Get dipole moment (in e·Bohr, convert to Debye)
    dipole_bohr = res.get_dipole()
    dipole_angstrom = dipole_bohr / BOHR_PER_ANGSTROM
    dipole_debye = dipole_angstrom * DEBYE_PER_E_ANGSTROM

    # Compute center of geometry as dipole origin
    origin = positions.mean(axis=0)

    return ChargeResult(
        charges=charges,
        dipole_vector=dipole_debye,
        dipole_magnitude=float(np.linalg.norm(dipole_debye)),
        dipole_origin=origin,
        charge_model="xtb-gfn2-mulliken",
        convention="physics",  # xTB uses physics convention (+ to -)
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

        # Extract charges
        charges = np.array(data.get("partial charges", []))

        # Extract dipole (in Debye)
        dipole = data.get("dipole moment", {})
        dipole_debye = np.array(
            [
                dipole.get("x", 0.0),
                dipole.get("y", 0.0),
                dipole.get("z", 0.0),
            ]
        )
        dipole_mag = dipole.get("total", float(np.linalg.norm(dipole_debye)))

        # Compute center of geometry
        positions = np.array([[a.x, a.y, a.z] for a in molecule.atoms])
        origin = positions.mean(axis=0)

        return ChargeResult(
            charges=charges,
            dipole_vector=dipole_debye,
            dipole_magnitude=dipole_mag,
            dipole_origin=origin,
            charge_model="xtb-gfn2",
            convention="physics",
        )


def calculate_charges_gasteiger(molecule: SDFMolecule) -> ChargeResult:
    """Calculate Gasteiger charges using RDKit.

    This is the fallback when xTB is not available.

    Args:
        molecule: Parsed SDF molecule

    Returns:
        ChargeResult with Gasteiger charges and computed dipole
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as e:
        raise ImportError("RDKit required for Gasteiger charges") from e

    # Build RDKit molecule from atoms (simple approach)
    positions = np.array([[a.x, a.y, a.z] for a in molecule.atoms])

    # Create editable molecule
    mol = Chem.RWMol()
    conf = Chem.Conformer(molecule.num_atoms)

    for i, atom in enumerate(molecule.atoms):
        rdatom = Chem.Atom(_element_to_number(atom.symbol))
        mol.AddAtom(rdatom)
        conf.SetAtomPosition(i, (atom.x, atom.y, atom.z))

    mol.AddConformer(conf, assignId=True)

    # Compute Gasteiger charges
    AllChem.ComputeGasteigerCharges(mol)

    charges = np.array(
        [
            float(mol.GetAtomWithIdx(i).GetDoubleProp("_GasteigerCharge"))
            for i in range(mol.GetNumAtoms())
        ]
    )

    # Handle NaN charges (can happen for some atoms)
    charges = np.nan_to_num(charges, nan=0.0)

    # Compute dipole moment from charges and positions
    # μ = Σ qᵢ rᵢ (in e·Å, then convert to Debye)
    dipole_e_angstrom = np.sum(charges[:, np.newaxis] * positions, axis=0)
    dipole_debye = dipole_e_angstrom * DEBYE_PER_E_ANGSTROM

    # Center of geometry as origin
    origin = positions.mean(axis=0)

    return ChargeResult(
        charges=charges,
        dipole_vector=dipole_debye,
        dipole_magnitude=float(np.linalg.norm(dipole_debye)),
        dipole_origin=origin,
        charge_model="rdkit-gasteiger",
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
    elements = {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Ga": 31,
        "Ge": 32,
        "As": 33,
        "Se": 34,
        "Br": 35,
        "Kr": 36,
        "Rb": 37,
        "Sr": 38,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "In": 49,
        "Sn": 50,
        "Sb": 51,
        "Te": 52,
        "I": 53,
        "Xe": 54,
    }
    return elements.get(symbol.capitalize(), 6)  # Default to carbon
