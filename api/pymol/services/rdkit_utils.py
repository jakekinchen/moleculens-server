"""RDKit utilities for molecular structure conversion and processing."""

import logging
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)


def sdf_to_pdb_block(sdf_data: str) -> str:
    """Convert SDF data (string) to a single PDB block using RDKit in-memory.

    Args:
        sdf_data: SDF format molecular structure data

    Returns:
        PDB format string, or empty string if conversion fails
    """
    try:
        mol = Chem.MolFromMolBlock(sdf_data, sanitize=True, removeHs=False)
        if mol is None:
            logger.warning("Failed to parse SDF data")
            return ""

        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        AllChem.MMFFOptimizeMolecule(mol)

        pdb_data = Chem.MolToPDBBlock(mol)
        return pdb_data if pdb_data else ""

    except Exception as e:
        logger.error(f"Error converting SDF to PDB: {str(e)}")
        return ""


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """Convert SMILES string to RDKit Mol object.

    Args:
        smiles: SMILES string representation

    Returns:
        RDKit Mol object or None if conversion fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Failed to parse SMILES: {smiles}")
        return mol
    except Exception as e:
        logger.error(f"Error parsing SMILES {smiles}: {str(e)}")
        return None


def mol_to_smarts(mol: Chem.Mol) -> Optional[str]:
    """Convert RDKit Mol object to SMARTS pattern.

    Args:
        mol: RDKit Mol object

    Returns:
        SMARTS pattern string or None if conversion fails
    """
    try:
        return Chem.MolToSmarts(mol)
    except Exception as e:
        logger.error(f"Error converting mol to SMARTS: {str(e)}")
        return None


def generate_2d_coords(mol: Chem.Mol) -> None:
    """Generate 2D coordinates for a molecule in-place.

    Args:
        mol: RDKit Mol object to modify
    """
    try:
        AllChem.Compute2DCoords(mol)
    except Exception as e:
        logger.error(f"Error generating 2D coordinates: {str(e)}")
