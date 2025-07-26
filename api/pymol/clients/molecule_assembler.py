"""Molecule assembler for converting CIDs to structured molecule data."""

import logging
from typing import Any

from rdkit import Chem

from ..services import (
    PubChemSearchService,
    generate_2d_coords,
    mol_to_smarts,
    sdf_to_pdb_block,
    smiles_to_mol,
)

logger = logging.getLogger(__name__)


class MoleculeAssembler:
    """Assembles comprehensive molecule data from PubChem CIDs."""

    def __init__(self, search_service: PubChemSearchService):
        """Initialize the molecule assembler.

        Args:
            search_service: PubChem search service instance
        """
        self.search_service = search_service

    def assemble_molecule_data(self, cid: int) -> dict[str, Any]:
        """Assemble comprehensive molecule data from a PubChem CID.

        Args:
            cid: PubChem Compound ID

        Returns:
            Dictionary containing comprehensive molecule data

        Raises:
            ValueError: If molecule data cannot be assembled
        """
        try:
            # Get compound details
            compound = self.search_service.get_compound_details(cid)
            if not compound:
                raise ValueError(f"Could not retrieve compound details for CID {cid}")

            # Fetch SDF data (try 3D first, fallback to 2D)
            sdf_data = self.search_service.fetch_sdf(cid, "3d")
            if not sdf_data:
                sdf_data = self.search_service.fetch_sdf(cid, "2d")
                if not sdf_data:
                    raise ValueError(f"Could not retrieve SDF data for CID {cid}")

            # Convert to PDB
            pdb_data = sdf_to_pdb_block(sdf_data)
            if not pdb_data:
                logger.warning(f"Failed to convert SDF to PDB for CID {cid}")

            # Process SMILES if available
            smarts_pattern = None
            if hasattr(compound, "isomeric_smiles") and compound.isomeric_smiles:
                mol = smiles_to_mol(compound.isomeric_smiles)
                if mol:
                    smarts_pattern = mol_to_smarts(mol)

            # Assemble comprehensive data
            molecule_data = {
                "name": getattr(compound, "iupac_name", str(cid)),
                "cid": cid,
                "pdb_data": pdb_data,
                "sdf": sdf_data,
                "smiles": getattr(compound, "isomeric_smiles", None),
                "canonical_smiles": getattr(compound, "canonical_smiles", None),
                "smarts_pattern": smarts_pattern,
                "iupac_name": getattr(compound, "iupac_name", None),
                "molecular_formula": getattr(compound, "molecular_formula", None),
                "molecular_weight": getattr(compound, "molecular_weight", None),
                "elements": getattr(compound, "elements", None),
                "charge": getattr(compound, "charge", None),
                "synonyms": getattr(compound, "synonyms", None),
            }

            # Add atomic details if available
            if hasattr(compound, "atoms") and compound.atoms:
                molecule_data["atoms"] = [
                    {
                        "number": getattr(atom, "number", None),
                        "element": getattr(atom, "element", None),
                        "x": getattr(atom, "x", None),
                        "y": getattr(atom, "y", None),
                        "z": getattr(atom, "z", None),
                        "charge": getattr(atom, "charge", None),
                    }
                    for atom in compound.atoms
                ]

            # Add bond details if available
            if hasattr(compound, "bonds") and compound.bonds:
                molecule_data["bonds"] = [
                    {
                        "aid1": getattr(bond, "aid1", None),
                        "aid2": getattr(bond, "aid2", None),
                        "order": getattr(bond, "order", None),
                        "style": getattr(bond, "style", None),
                    }
                    for bond in compound.bonds
                ]

            return molecule_data

        except Exception as e:
            logger.error(f"Error assembling molecule data for CID {cid}: {str(e)}")
            raise ValueError(f"Failed to assemble molecule data for CID {cid}: {str(e)}") from e

    def assemble_2d_molecule_data(self, cid: int) -> dict[str, Any]:
        """Assemble 2D molecule data for diagram rendering.

        Args:
            cid: PubChem Compound ID

        Returns:
            Dictionary containing 2D molecular structure data

        Raises:
            ValueError: If 2D data cannot be assembled
        """
        try:
            # Get compound details
            compound = self.search_service.get_compound_details(cid)
            if not compound:
                raise ValueError(f"Could not retrieve compound details for CID {cid}")

            # Fetch 2D SDF data
            sdf_data = self.search_service.fetch_sdf(cid, "2d")
            if not sdf_data:
                raise ValueError(f"Could not retrieve 2D SDF data for CID {cid}")

            # Parse with RDKit to get 2D coordinates
            mol = Chem.MolFromMolBlock(sdf_data, sanitize=True, removeHs=False)
            if mol is None:
                raise ValueError("Unable to parse SDF data with RDKit")

            if mol.GetNumConformers() == 0:
                generate_2d_coords(mol)

            conf = mol.GetConformer()
            atoms = []
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                atoms.append(
                    {
                        "element": atom.GetSymbol(),
                        "x": pos.x,
                        "y": pos.y,
                    }
                )

            bonds = []
            for bond in mol.GetBonds():
                bonds.append(
                    {
                        "start": bond.GetBeginAtomIdx(),
                        "end": bond.GetEndAtomIdx(),
                        "order": bond.GetBondTypeAsDouble(),
                    }
                )

            return {
                "atoms": atoms,
                "bonds": bonds,
                "name": getattr(compound, "iupac_name", str(cid)),
                "cid": cid,
                "formula": getattr(compound, "molecular_formula", ""),
            }

        except Exception as e:
            logger.error(f"Error assembling 2D molecule data for CID {cid}: {str(e)}")
            raise ValueError(f"Failed to assemble 2D molecule data for CID {cid}: {str(e)}") from e

    def assemble_molecules_2d_layout(self, queries: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Assemble 2D data for multiple molecules with layout information.

        Args:
            queries: List of dictionaries with 'query' and 'box' keys

        Returns:
            List of dictionaries with molecule data and box information
        """
        layout = []
        for item in queries:
            query = item.get("query")
            box = item.get("box")
            if query is None:
                continue

            try:
                # Search for the molecule first
                compounds = self.search_service.search_with_fallbacks(query)
                if not compounds:
                    logger.warning(f"No compounds found for query: {query}")
                    continue

                cid = compounds[0].cid
                molecule_data = self.assemble_2d_molecule_data(cid)
                layout.append(
                    {
                        **molecule_data,
                        "box": box,
                        "query": query,
                    }
                )
            except Exception as e:
                logger.error(f"Error processing query '{query}': {str(e)}")
                continue

        return layout
