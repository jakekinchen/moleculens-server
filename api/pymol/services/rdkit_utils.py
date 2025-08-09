"""RDKit utilities for molecular structure conversion and processing."""

import logging
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor

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
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())  # type: ignore[attr-defined]

        AllChem.MMFFOptimizeMolecule(mol)  # type: ignore[attr-defined]

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
        rdDepictor.Compute2DCoords(mol)
    except Exception as e:
        logger.error(f"Error generating 2D coordinates: {str(e)}")


def _compile_smarts(smarts: str) -> Optional[Chem.Mol]:
    """Safely compile a SMARTS pattern to an RDKit Mol.

    Returns None if compilation fails.
    """
    try:
        pattern = Chem.MolFromSmarts(smarts)
        return pattern
    except Exception as exc:  # pragma: no cover - defensive
        logger.error(f"Failed to compile SMARTS '{smarts}': {exc}")
        return None


def get_functional_group_patterns() -> Dict[str, Chem.Mol]:
    """Return common functional group SMARTS patterns as RDKit Mol objects.

    Notes:
        Patterns are intentionally simple/representative to balance recall and precision.
    """
    patterns: Dict[str, str] = {
        # Oxygen-containing
        "alcohol": "[CX4;!$(C=*)][OX2H]",
        "phenol": "[cX3H][OX2H]",
        "aldehyde": "[CX3H1](=O)[#6]",
        "ketone": "[#6][CX3](=O)[#6]",
        "carboxylic_acid": "[CX3](=O)[OX2H1]",
        "ester": "[CX3](=O)[OX2][#6]",
        "ether": "[OD2]([#6])[#6]",
        "carbonate": "[OX2][CX3](=O)[OX2][#6]",
        # Nitrogen-containing
        "amine_primary": "[NX3;H2][#6]",
        "amine_secondary": "[NX3;H1]([#6])[#6]",
        "amine_tertiary": "[NX3;H0]([#6])([#6])[#6]",
        "amide": "[NX3][CX3](=O)[#6]",
        "imine": "[CX2]=[NX3]",
        "urea": "[NX3][CX3](=O)[NX3]",
        # Halogens and related
        "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
        "alkyl_halide": "[F,Cl,Br,I][CX4]",
        # Sulfur-containing
        "thiol": "[SX2H]",
        "disulfide": "[SX2]-S-[SX2]",
        "sulfoxide": "[SX3](=O)[#6]",
        "sulfone": "[SX4](=O)(=O)[#6]",
        # Phosphorus-containing
        "phosphate": "[PX4](=O)([OX2])[OX2]",
        # Unsaturations / rings
        "alkene": "[CX3]=[CX3]",
        "alkyne": "[CX2]#[CX2]",
        "aromatic_ring": "a1aaaaa1",
        # Other
        "anhydride": "[CX3](=O)O[CX3](=O)",
        "nitrile": "[CX2]#N",
        "nitro": "[N+](=O)[O-]",
        "ether_aromatic": "[OD2]([#6])[a]",
    }

    compiled: Dict[str, Chem.Mol] = {}
    for name, smarts in patterns.items():
        patt = _compile_smarts(smarts)
        if patt is not None:
            compiled[name] = patt
    return compiled


def identify_functional_groups(
    mol: Chem.Mol, include_matches: bool = False
) -> Tuple[Dict[str, int], List[Tuple[str, Tuple[int, ...]]]]:
    """Identify functional groups present in a molecule via SMARTS matching.

    Args:
        mol: RDKit Mol
        include_matches: When True, return per-match atom indices.

    Returns:
        A tuple of (group_counts, match_list). `match_list` contains tuples of
        (group_name, atom_indices) and is empty when include_matches is False.
    """
    if mol is None:
        return {}, []

    patterns = get_functional_group_patterns()
    group_counts: Dict[str, int] = {}
    match_list: List[Tuple[str, Tuple[int, ...]]] = []

    for name, patt in patterns.items():
        try:
            matches = mol.GetSubstructMatches(patt, uniquify=True)
        except Exception as exc:  # pragma: no cover - defensive
            logger.error(f"Substructure search failed for {name}: {exc}")
            continue

        count = len(matches)
        if count > 0:
            group_counts[name] = count
            if include_matches:
                for m in matches:
                    match_list.append((name, m))

    return group_counts, match_list
