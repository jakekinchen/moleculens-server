"""
PubChem Agent - Handles interaction with PubChem database and molecular structure retrieval.
"""

import datetime
import json
import logging
import os
import re
import tempfile
import urllib.parse
from typing import Any, Dict, List, NamedTuple, Optional

import pubchempy as pcp
import requests
from agent_management.agents.pubchem_agent_helper import validate_and_convert_script
from agent_management.agents.script_agent import ScriptAgent
from agent_management.debug_utils import DEBUG_PUBCHEM, write_debug_file
from agent_management.llm_service import (
    LLMModelConfig,
    LLMRequest,
    LLMService,
    ProviderType,
    StructuredLLMRequest,
)
from agent_management.models import PubChemCompound, PubChemSearchResult
from agent_management.molecule_visualizer import MoleculeVisualizer
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Fragments


def _sdf_to_pdb_block(sdf_data: str) -> str:
    """
    Convert SDF data (string) to a single PDB block using RDKit in-memory.

    Returns an empty string if conversion fails.
    """
    mol = Chem.MolFromMolBlock(sdf_data, sanitize=True, removeHs=False)
    if mol is None:
        return ""

    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    AllChem.MMFFOptimizeMolecule(mol)

    pdb_data = Chem.MolToPDBBlock(mol)
    return pdb_data if pdb_data else ""


class MoleculePackage(NamedTuple):
    """Container for molecule visualization package data"""

    pdb_data: str
    html: str
    title: str


class PubChemAgent:
    """Agent for retrieving and processing molecular data from PubChem"""

    FORMULA_EXTRACTION_PROMPT = """Given a chemical compound name, provide its molecular formula.
Only respond with the formula, nothing else. If you're not certain, respond with 'unknown'.

Example inputs and outputs:
Input: "water"
Output: H2O

Input: "glucose"
Output: C6H12O6

Input: "random chemical"
Output: unknown

Input: "{query}"
Output:"""

    def __init__(
        self,
        llm_service: LLMService,
        use_element_labels: bool = True,
        convert_back_to_indices: bool = False,
        script_model: Optional[str] = None,
    ):
        """
        Initialize the PubChem agent

        Args:
            llm_service: LLM service for generating interpretations and scripts
            use_element_labels: Whether to use element-based labels (C1, O1) instead of indices
            convert_back_to_indices: Whether to convert element-labels back to numeric indices
            script_model: Optional model override for script agent
        """
        self.llm_service = llm_service
        self.use_element_labels = use_element_labels
        self.convert_back_to_indices = convert_back_to_indices  # Convert back to numeric indices after script generation
        self.script_model = script_model  # Optional model override for script agent
        self.logger = logging.getLogger(__name__)

    def _normalize_query(self, query: str) -> List[str]:
        """
        Normalize a chemical query string to improve search results.

        Args:
            query: Raw query string

        Returns:
            List of normalized query strings to try
        """
        # Remove any extra whitespace
        query = query.strip()

        # List to store variations
        variations = []

        # Add the original query
        variations.append(query)

        # Handle special cases for chemical notation
        if "[" in query and "]" in query:
            # For cases like [2]Rotaxane or cryptand[2.2.2], try multiple formats
            no_space = re.sub(r"\[\s*(\d+)\s*\]", r"[\1]", query)
            with_space = re.sub(r"\[\s*(\d+)\s*\]", r"[ \1 ]", query)
            no_brackets = re.sub(r"\[\s*(\d+)\s*\]", r"\1-", query)
            variations.extend([no_space, with_space, no_brackets])

            # For cryptands, try alternative notations
            if "cryptand" in query.lower():
                cryptand_base = query.lower().replace("cryptand", "").strip("[]")
                variations.extend(
                    [
                        f"Cryptand-{cryptand_base}",
                        f"Kryptofix {cryptand_base}",
                        "4,7,13,16,21,24-Hexaoxa-1,10-diazabicyclo[8.8.8]hexacosane",  # Common cryptand
                    ]
                )

        # Handle complex names with multiple parts
        if " and " in query.lower() or " with " in query.lower():
            # Try the first part of compound names
            parts = re.split(r"\s+(?:and|with)\s+", query, flags=re.IGNORECASE)
            if parts:
                variations.append(parts[0].strip())

        # Handle specific chemical classes
        if "phosphazene" in query.lower():
            variations.extend(
                [
                    "Hexachlorocyclotriphosphazene",  # Common phosphazene
                    "Cyclophosphazene",
                    "Phosphonitrilic chloride trimer",
                ]
            )

        if "crown ether" in query.lower():
            variations.extend(["18-crown-6", "15-crown-5", "12-crown-4"])

        if "metal-carbonyl" in query.lower():
            variations.extend(["Fe(CO)5", "Ni(CO)4", "Cr(CO)6"])

        # Remove duplicates while preserving order
        seen = set()
        unique_variations = []
        for v in variations:
            if v not in seen:
                seen.add(v)
                unique_variations.append(v)

        return unique_variations

    def _search_pubchem_rest(
        self, query: str, max_results: int = 5
    ) -> List[Dict[str, Any]]:
        """
        Search PubChem using the REST API directly.

        Args:
            query: Search query
            max_results: Maximum number of results to return

        Returns:
            List of compound data dictionaries
        """
        self.logger.info(f"[DEBUG] Searching PubChem REST API for: {query}")

        # Encode the query for URL
        encoded_query = urllib.parse.quote(query)

        # First, search for compounds matching the query
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_query}/cids/JSON"

        try:
            search_response = requests.get(search_url, timeout=30)
            if search_response.status_code != 200:
                self.logger.warning(
                    f"PubChem REST search failed with status {search_response.status_code}"
                )
                return []

            data = search_response.json()
            if "IdentifierList" not in data:
                return []

            cids = data["IdentifierList"].get("CID", [])
            if not cids:
                return []

            # Limit the number of CIDs
            cids = cids[:max_results]

            # Get detailed information for each CID
            compounds = []
            for cid in cids:
                try:
                    # Get compound details using pubchempy for consistency
                    compound = pcp.Compound.from_cid(cid)
                    compounds.append(compound)
                except Exception as e:
                    self.logger.warning(
                        f"Failed to get details for CID {cid}: {str(e)}"
                    )
                    continue

            return compounds

        except Exception as e:
            self.logger.error(f"Error in PubChem REST search: {str(e)}")
            return []

    def _search_pubchem_direct(self, query: str) -> List[Dict[str, Any]]:
        """
        Search PubChem using direct compound search that matches the web interface.

        Args:
            query: Search query (e.g. 'cryptand[2.2.2]')

        Returns:
            List of compound data dictionaries
        """
        self.logger.info(f"[DEBUG] Direct PubChem search for: {query}")

        try:
            # First try exact name match
            search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(query)}/JSON"
            response = requests.get(search_url, timeout=30)

            if response.status_code == 200:
                data = response.json()
                if "PC_Compounds" in data:
                    cid = data["PC_Compounds"][0]["id"]["id"]["cid"]
                    self.logger.info(f"[DEBUG] Found exact match CID: {cid}")
                    compound = pcp.Compound.from_cid(int(cid))
                    return [compound]

            # If exact match fails, try the autocomplete API
            autocomplete_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{urllib.parse.quote(query)}/json"
            response = requests.get(autocomplete_url, timeout=30)

            if response.status_code == 200:
                data = response.json()
                if (
                    "dictionary_terms" in data
                    and "compound" in data["dictionary_terms"]
                ):
                    # Get the first suggestion's CID
                    for suggestion in data["dictionary_terms"]["compound"]:
                        if "cid" in suggestion:
                            cid = suggestion["cid"]
                            self.logger.info(f"[DEBUG] Found suggested CID: {cid}")
                            compound = pcp.Compound.from_cid(int(cid))
                            return [compound]

            return []

        except Exception as e:
            self.logger.error(f"Error in direct PubChem search: {str(e)}")
            return []

    def interpret_user_query(self, user_input: str) -> str:
        """
        Use LLM to interpret user input into a molecule name or identifier.
        Handles provider failures gracefully with fallback behavior.

        Args:
            user_input: The user's query about a molecule

        Returns:
            A molecule name or identifier, or the original input if interpretation fails

        Raises:
            Exception: If all interpretation attempts fail and the input cannot be used directly
        """
        try:
            request = LLMRequest(
                user_prompt=f"""Given this user input is related to molecular structures or trying to learn something about the microscopic structures of molecules or some topic related, give us the molecule that would best help the user learn what they are trying to learn about, return 'N/A' if you can't determine a valid molecule.

Important guidelines:
1. Prefer simple, well-known examples over complex ones
2. For complex structures, choose a representative simple example
3. Use common names over systematic names when possible
4. Avoid special characters unless absolutely necessary
5. For metal complexes, prefer simple coordination compounds
6. For organic structures, prefer parent compounds over derivatives

Examples:
- For "metal carbonyl complexes" → "Fe(CO)5" (not a complex derivative)
- For "crown ethers" → "18-crown-6" (the most common example)
- For "phosphazenes" → "hexachlorocyclophosphazene" (the parent compound)
- For "cryptands" → "cryptand-222" (standard notation)

User input: '{user_input}'
Only respond with the molecule name or 'N/A', no other text.""",
                system_prompt="You are a chemistry professor that helps identify the simplest, most representative molecule names from user queries. Always prefer well-known, simple examples that clearly demonstrate the concept.",
            )
            response = self.llm_service.generate(request)
            return response.content.strip()
        except Exception as e:
            self.logger.warning(f"Failed to interpret user query with LLM: {str(e)}")

            # Fallback 1: Try to use the input directly if it looks like a molecule name
            if len(user_input.split()) <= 3 and not any(
                char in user_input for char in "?!.,:;"
            ):
                self.logger.info(
                    f"Using user input directly as molecule name: '{user_input}'"
                )
                return user_input

            # Fallback 2: Use a default molecule if the input is complex
            self.logger.info(
                f"Using default molecule 'water' for complex query: '{user_input}'"
            )
            return "water"

    def fetch_sdf_for_cid(self, cid: int) -> Optional[str]:
        """
        Fetches the SDF string for a given PubChem CID using a direct REST call.
        Returns SDF as text, or None if something fails.
        """
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            return resp.text
        return None

    def _get_molecular_formula(self, compound_name: str) -> Optional[str]:
        """
        Use LLM to get the molecular formula for a compound name.

        Args:
            compound_name: Name of the chemical compound

        Returns:
            Molecular formula if LLM is confident, None otherwise
        """
        try:
            request = LLMRequest(
                system_prompt="You are a chemistry expert. Provide molecular formulas for chemical compounds.",
                user_prompt=self.FORMULA_EXTRACTION_PROMPT.format(query=compound_name),
                temperature=0.1,  # Low temperature for more deterministic output
            )

            response = self.llm_service.generate(request)
            formula = response.content.strip()

            if formula.lower() == "unknown":
                return None

            return formula

        except Exception as e:
            self.logger.error(f"Error getting molecular formula from LLM: {str(e)}")
            return None

    def _search_with_fallbacks(self, query: str) -> List[Dict[str, Any]]:
        """
        Search for a compound using multiple methods with fallbacks.

        Args:
            query: The compound name or identifier to search for

        Returns:
            List of compound data dictionaries
        """
        self.logger.info(f"Starting search with fallbacks for: {query}")

        # Step 1: Try direct search with normalized queries
        normalized_queries = self._normalize_query(query)
        for normalized_query in normalized_queries:
            try:
                results = self._search_pubchem_direct(normalized_query)
                if results:
                    self.logger.info(
                        f"Found results using direct search with query: {normalized_query}"
                    )
                    return results
            except Exception as e:
                self.logger.warning(
                    f"Direct search failed for {normalized_query}: {str(e)}"
                )

        # Step 2: Try REST search with normalized queries
        for normalized_query in normalized_queries:
            try:
                results = self._search_pubchem_rest(normalized_query)
                if results:
                    self.logger.info(
                        f"Found results using REST search with query: {normalized_query}"
                    )
                    return results
            except Exception as e:
                self.logger.warning(
                    f"REST search failed for {normalized_query}: {str(e)}"
                )

        # Step 3: Try getting formula from LLM and searching with that
        formula = self._get_molecular_formula(query)
        if formula:
            self.logger.info(f"Got formula from LLM: {formula}")
            try:
                # Try searching by formula using PubChemPy
                compounds = pcp.get_compounds(formula, "formula")
                if compounds:
                    self.logger.info(f"Found results using formula search: {formula}")
                    return compounds
            except Exception as e:
                self.logger.warning(f"Formula search failed for {formula}: {str(e)}")

        self.logger.warning(f"All search methods failed for query: {query}")
        return []

    def get_molecule_sdfs(self, user_input: str) -> List[Dict[str, Any]]:
        """
        Get SDF data for molecules based on user input.

        Args:
            user_input: User's query about a molecule

        Returns:
            List of dictionaries containing molecule data and SDF content
        """
        self.logger.info(f"Processing request for: {user_input}")

        # First interpret the user's query to get the molecule name
        molecule_name = self.interpret_user_query(user_input)
        if not molecule_name:
            self.logger.warning("Could not interpret user query into a molecule name")
            return []

        # Search for the molecule using our fallback pipeline
        compounds = self._search_with_fallbacks(molecule_name)
        if not compounds:
            self.logger.warning(f"No compounds found for {molecule_name}")
            return []

        # Process the compounds and get their SDF data
        results = []
        for compound in compounds:
            try:
                # Handle both Compound objects and dictionaries
                cid = compound.cid if hasattr(compound, "cid") else compound.get("cid")
                if not cid:
                    continue

                # First try to get 3D SDF
                sdf_data = None
                sdf_3d_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
                try:
                    response = requests.get(sdf_3d_url, timeout=30)
                    if response.status_code == 200:
                        sdf_data = response.text
                        self.logger.info(f"Successfully got 3D SDF for CID {cid}")
                    else:
                        self.logger.warning(
                            f"Failed to get 3D SDF for CID {cid} (status {response.status_code})"
                        )
                except Exception as e:
                    self.logger.warning(f"Error getting 3D SDF for CID {cid}: {str(e)}")

                # If 3D failed, try 2D
                if not sdf_data:
                    sdf_2d_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
                    try:
                        response = requests.get(sdf_2d_url, timeout=30)
                        if response.status_code == 200:
                            sdf_data = response.text
                            self.logger.info(f"Successfully got 2D SDF for CID {cid}")
                        else:
                            self.logger.warning(
                                f"Failed to get 2D SDF for CID {cid} (status {response.status_code})"
                            )
                    except Exception as e:
                        self.logger.warning(
                            f"Error getting 2D SDF for CID {cid}: {str(e)}"
                        )

                if sdf_data:
                    # Get compound attributes safely
                    name = (
                        compound.iupac_name
                        if hasattr(compound, "iupac_name")
                        else compound.get("iupac_name", str(cid))
                    )
                    formula = (
                        compound.molecular_formula
                        if hasattr(compound, "molecular_formula")
                        else compound.get("formula", "")
                    )

                    results.append(
                        {"name": name, "cid": cid, "formula": formula, "sdf": sdf_data}
                    )
                    self.logger.info(f"Successfully processed compound CID {cid}")
                else:
                    self.logger.error(f"Failed to get any SDF data for CID {cid}")

            except Exception as e:
                self.logger.error(f"Error processing compound: {str(e)}")
                continue

        if not results:
            self.logger.warning(
                f"No valid compounds with SDF data found for {molecule_name}"
            )

        return results

    def get_compound_details(self, cid: int) -> Optional[PubChemCompound]:
        """
        Get detailed information for a specific compound by CID.
        """
        try:
            compound = pcp.Compound.from_cid(cid)
            sdf_str = self.fetch_sdf_for_cid(cid)

            return PubChemCompound(
                name=compound.iupac_name,  # Use IUPAC name as the primary name
                cid=compound.cid,
                molecular_formula=compound.molecular_formula,
                molecular_weight=compound.molecular_weight,
                iupac_name=compound.iupac_name,
                sdf=sdf_str,
                canonical_smiles=compound.canonical_smiles,
                isomeric_smiles=compound.isomeric_smiles,
                elements=compound.elements,
                atoms=compound.atoms,
                bonds=compound.bonds,
                charge=compound.charge,
                synonyms=compound.synonyms,
            )
        except Exception as e:
            print(f"Error fetching compound details for CID {cid}: {str(e)}")
            return None

    def save_compound_details_to_json(
        self, cid: int, base_dir: str = "compound_data"
    ) -> str:
        """
        Save all available compound details to a JSON file in a subdirectory.

        Args:
            cid: PubChem Compound ID
            base_dir: Base directory for saving compound data (default: "compound_data")

        Returns:
            Path to the created JSON file

        Raises:
            ValueError: If compound details cannot be fetched
        """
        try:
            # Create the compound directory
            compound_dir = os.path.join(base_dir, f"cid_{cid}")
            os.makedirs(compound_dir, exist_ok=True)

            # Get compound details
            compound = pcp.Compound.from_cid(cid)
            sdf_str = self.fetch_sdf_for_cid(cid)

            # Get 3D SDF if available
            sdf_3d_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
            response = requests.get(sdf_3d_url, timeout=30)
            sdf_3d = response.text if response.status_code == 200 else None

            # Create a comprehensive data structure
            data = {
                "metadata": {
                    "cid": cid,
                    "timestamp": datetime.datetime.now().isoformat(),
                    "data_source": "PubChem API",
                },
                "basic_info": {
                    "name": compound.iupac_name,
                    "iupac_name": compound.iupac_name,
                    "molecular_formula": compound.molecular_formula,
                    "molecular_weight": compound.molecular_weight,
                    "exact_mass": getattr(compound, "exact_mass", None),
                    "monoisotopic_mass": getattr(compound, "monoisotopic_mass", None),
                    "charge": compound.charge,
                    "complexity": getattr(compound, "complexity", None),
                },
                "structure": {
                    "smiles": compound.isomeric_smiles,
                    "canonical_smiles": compound.canonical_smiles,
                    "sdf_2d": sdf_str,
                    "sdf_3d": sdf_3d,
                    "inchi": getattr(compound, "inchi", None),
                    "inchikey": getattr(compound, "inchikey", None),
                    "fingerprint": getattr(compound, "fingerprint", None),
                },
                "composition": {
                    "elements": compound.elements,
                    "atoms_count": len(compound.atoms) if compound.atoms else 0,
                    "bonds_count": len(compound.bonds) if compound.bonds else 0,
                    "rotatable_bond_count": getattr(
                        compound, "rotatable_bond_count", None
                    ),
                    "heavy_atom_count": getattr(compound, "heavy_atom_count", None),
                    "isotope_atom_count": getattr(compound, "isotope_atom_count", None),
                    "atom_stereo_count": getattr(compound, "atom_stereo_count", None),
                    "defined_atom_stereo_count": getattr(
                        compound, "defined_atom_stereo_count", None
                    ),
                    "undefined_atom_stereo_count": getattr(
                        compound, "undefined_atom_stereo_count", None
                    ),
                    "bond_stereo_count": getattr(compound, "bond_stereo_count", None),
                    "defined_bond_stereo_count": getattr(
                        compound, "defined_bond_stereo_count", None
                    ),
                    "undefined_bond_stereo_count": getattr(
                        compound, "undefined_bond_stereo_count", None
                    ),
                    "covalent_unit_count": getattr(
                        compound, "covalent_unit_count", None
                    ),
                },
                "properties": {
                    "h_bond_donor_count": getattr(compound, "h_bond_donor_count", None),
                    "h_bond_acceptor_count": getattr(
                        compound, "h_bond_acceptor_count", None
                    ),
                    "tpsa": getattr(compound, "tpsa", None),
                    "xlogp": getattr(compound, "xlogp", None),
                    "molecular_volume": getattr(compound, "molecular_volume", None),
                    "polarizability": getattr(compound, "polarizability", None),
                },
                "identifiers": {
                    "synonyms": compound.synonyms,
                    "mesh_entries": getattr(compound, "mesh_entries", None),
                    "record_type": getattr(compound, "record_type", None),
                },
            }

            # Add atom details if available
            if compound.atoms:
                data["structure"]["atoms"] = [
                    {
                        "number": atom.number if hasattr(atom, "number") else None,
                        "element": atom.element if hasattr(atom, "element") else None,
                        "x": getattr(atom, "x", None),
                        "y": getattr(atom, "y", None),
                        "z": getattr(atom, "z", None),
                        "charge": getattr(atom, "charge", None),
                    }
                    for atom in compound.atoms
                ]

            # Add bond details if available
            if compound.bonds:
                data["structure"]["bonds"] = [
                    {
                        "aid1": bond.aid1,
                        "aid2": bond.aid2,
                        "order": bond.order,
                        "style": getattr(bond, "style", None),
                    }
                    for bond in compound.bonds
                ]

            # Save to JSON file
            json_path = os.path.join(compound_dir, "compound_details.json")
            with open(json_path, "w") as f:
                json.dump(data, f, indent=2)

            self.logger.info(f"Saved compound details to {json_path}")
            return json_path

        except Exception as e:
            self.logger.error(f"Error saving compound details for CID {cid}: {str(e)}")
            raise

    def get_molecule_data(self, user_query: str) -> Dict[str, Any]:
        """
        Step A:
        Fetch raw molecule data (e.g., PDB block, name, elements, etc.) without generating HTML.

        Returns a dictionary with minimal data:
         {
           "pdb_data": str,
           "name": str,
           "cid": int,
           "formula": str,
           "atoms": [...],
           "bonds": [...],
           ...
         }
        """
        self.logger.info(f"[DEBUG] Fetching molecule data for: {user_query}")
        # 1) Interpret user query
        molecule_name = self.interpret_user_query(user_query)
        if not molecule_name:
            raise ValueError("Could not interpret user query into a molecule name.")

        # 2) Search for the molecule
        compounds = self._search_with_fallbacks(molecule_name)
        if not compounds:
            raise ValueError(f"No compounds found for {molecule_name}")

        # 3) Use the first compound
        compound = compounds[0]
        cid = compound.cid if hasattr(compound, "cid") else compound.get("cid")
        if not cid:
            raise ValueError("Compound has no valid CID.")

        # 4) Fetch 3D SDF or fallback to 2D
        sdf_data = None
        sdf_3d_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        try:
            response = requests.get(sdf_3d_url, timeout=30)
            if response.status_code == 200:
                sdf_data = response.text
            else:
                self.logger.warning(
                    f"Failed 3D SDF for CID {cid}, status {response.status_code}"
                )
        except Exception as e:
            self.logger.warning(f"Error getting 3D SDF for CID {cid}: {str(e)}")

        if not sdf_data:
            sdf_2d_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
            )
            try:
                response = requests.get(sdf_2d_url, timeout=30)
                if response.status_code == 200:
                    sdf_data = response.text
                else:
                    self.logger.warning(
                        f"Failed 2D SDF for CID {cid}, status {response.status_code}"
                    )
            except Exception as e:
                self.logger.warning(f"Error getting 2D SDF for CID {cid}: {str(e)}")

        if not sdf_data:
            raise ValueError(f"Unable to retrieve any SDF data for CID {cid}")

        # 5) Convert to PDB block (for easy visualization), gather info
        pdb_data = _sdf_to_pdb_block(sdf_data)
        name = (
            compound.iupac_name
            if hasattr(compound, "iupac_name")
            else compound.get("iupac_name", str(cid))
        )
        formula = (
            compound.molecular_formula
            if hasattr(compound, "molecular_formula")
            else compound.get("formula", "")
        )

        # 6) Return essential info as a dict
        return {
            "pdb_data": pdb_data,
            "name": name,
            "cid": cid,
            "formula": formula,
            "sdf": sdf_data,
            # Optionally add more details if desired (atoms, synonyms, etc.)
            # For advanced usage, see self.get_compound_details
        }

    def get_molecule_2d_info(self, user_query: str) -> Dict[str, Any]:
        """Fetch 2D structural information for a molecule.

        This is similar to :meth:`get_molecule_data` but instead of returning a
        PDB block intended for 3D visualization it provides atomic coordinates
        in two dimensions for drawing a 2D diagram.

        Args:
            user_query: User supplied text describing the molecule.

        Returns:
            Dictionary with keys ``atoms`` and ``bonds`` containing 2D
            coordinates and connectivity information as well as ``name``,
            ``cid`` and ``formula``.
        """

        self.logger.info(f"[DEBUG] Fetching 2D molecule info for: {user_query}")

        molecule_name = self.interpret_user_query(user_query)
        if not molecule_name:
            raise ValueError("Could not interpret user query into a molecule name.")

        compounds = self._search_with_fallbacks(molecule_name)
        if not compounds:
            raise ValueError(f"No compounds found for {molecule_name}")

        compound = compounds[0]
        cid = compound.cid if hasattr(compound, "cid") else compound.get("cid")
        if not cid:
            raise ValueError("Compound has no valid CID.")

        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        response = requests.get(sdf_url, timeout=30)
        if response.status_code != 200:
            raise ValueError(f"Failed to get 2D SDF for CID {cid}")

        sdf_data = response.text

        mol = Chem.MolFromMolBlock(sdf_data, sanitize=True, removeHs=False)
        if mol is None:
            raise ValueError("Unable to parse SDF data")

        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)

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

        name = (
            compound.iupac_name
            if hasattr(compound, "iupac_name")
            else compound.get("iupac_name", str(cid))
        )
        formula = (
            compound.molecular_formula
            if hasattr(compound, "molecular_formula")
            else compound.get("formula", "")
        )

        return {
            "atoms": atoms,
            "bonds": bonds,
            "name": name,
            "cid": cid,
            "formula": formula,
        }

    def get_molecules_2d_layout(
        self, queries: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Fetch 2D info for multiple molecules and attach layout boxes.

        Args:
            queries: List of dictionaries each containing ``query`` and ``box``
                keys describing the molecule to fetch and its desired placement
                rectangle.

        Returns:
            List of dictionaries with molecule data combined with the provided
            box information under the key ``box``.
        """

        layout = []
        for item in queries:
            q = item.get("query")
            box = item.get("box")
            molecule = self.get_molecule_2d_info(q)
            layout.append(
                {
                    **molecule,
                    "box": box,
                    "query": q,
                }
            )
        return layout

    def generate_visualization(self, molecule_data: Dict[str, Any]) -> str:
        """
        Step B:
        Given existing molecule data, run script agent & generate HTML visualization
        without refetching from PubChem or re-running earlier logic.

        Args:
            molecule_data: Dictionary containing at least 'pdb_data', 'name', 'cid', 'sdf', etc.

        Returns:
            str: The complete HTML for an interactive 3D visualization
        """
        try:
            pdb_data = molecule_data["pdb_data"]
            display_title = molecule_data.get("name", "Molecule")
            cid = molecule_data.get("cid")

            # Minimal creation of a more detailed structure for script generation:
            # If advanced data is needed, parse from 'sdf' or re-check the compound details
            molecule_dict = {
                "name": display_title,
                "cid": cid,
                "sdf": molecule_data.get("sdf"),
            }
            # If you want advanced atomic details, parse them from 'sdf' or previously stored data
            # For now we keep it simple

            # 1) Create script agent & generate script
            script_agent = ScriptAgent(llm_service=self.llm_service)
            script = script_agent.generate_script_from_molecule(
                display_title, f"User query for CID {cid}", molecule_dict
            )

            # 2) Validate and convert script if needed
            script = validate_and_convert_script(script)

            # 3) Produce HTML using the visualizer
            visualizer = MoleculeVisualizer()
            html = visualizer.generate_interactive_html(
                pdb_data=pdb_data, title=display_title, script_data=script
            )
            return html

        except KeyError as e:
            raise ValueError(f"Missing required field in molecule_data: {str(e)}")
        except Exception as e:
            self.logger.error(f"Error generating visualization: {str(e)}")
            raise ValueError(f"Could not generate HTML visualization: {str(e)}")

    def get_molecule_package(self, user_query: str) -> MoleculePackage:
        """
        Get a complete molecule visualization package from a user query.

        # Write a debug file to confirm this method is being called
        if DEBUG_PUBCHEM:
            debug_info = {
                'method': 'get_molecule_package',
                'user_query': user_query,
                'timestamp': datetime.datetime.now().isoformat()
            }
            write_debug_file('pubchem_package_start.json', json.dumps(debug_info, indent=2))

        self.logger.info(f"[DEBUG] Getting molecule package for query: {user_query}")

        try:
            # First, get the molecule SDFs
            search_result = self.get_molecule_sdfs(user_query)

            if not search_result:
                raise ValueError(f"No molecules found for query: {user_query}, search result: {search_result}")

            # Use the first compound
            compound = search_result[0]

            if compound['sdf'] is None:
                raise ValueError(f"No SDF data available for compound {compound['name']} (CID: {compound['cid']})")

            pdb_data = _sdf_to_pdb_block(compound['sdf'])

            # Get a display title (use name if available, otherwise ID)
            display_title = compound['name'] if compound['name'] else f"CID {compound['cid']}"
            self.logger.info(f"[DEBUG] Using molecule: {display_title}")

            if DEBUG_PUBCHEM:
                write_debug_file('pubchem_sdf.txt', compound['sdf'] or '')

            try:
                # Generate a display title for the molecule
                display_title = compound['name'] if compound['name'] else "Molecule"
                self.logger.info(f"Using display title: {display_title}")

                # Get compound details using PubChemPy for additional data
                compound_details = pcp.Compound.from_cid(compound['cid'])

                if compound_details.isomeric_smiles:
                    mol = Chem.MolFromSmiles(compound_details.isomeric_smiles)
                    if mol is None:
                        raise ValueError(f"Failed to parse SMILES for compound {compound_details.name} (CID: {compound_details.cid})")

                    smarts_pattern = Chem.MolToSmarts(mol)

                    # NEW: Identify functional groups using RDKit fragment functions

                molecule_data = {
                    'name': compound['name'],
                    'cid': compound['cid'],
                    'smiles': compound_details.isomeric_smiles if hasattr(compound_details, 'isomeric_smiles') else None,
                    'smarts_pattern': smarts_pattern if 'smarts_pattern' in locals() else None,
                    'iupac_name': compound_details.iupac_name if hasattr(compound_details, 'iupac_name') else None,
                    'molecular_formula': compound_details.molecular_formula if hasattr(compound_details, 'molecular_formula') else None,
                    'molecular_weight': compound_details.molecular_weight if hasattr(compound_details, 'molecular_weight') else None,
                    'elements': compound_details.elements if hasattr(compound_details, 'elements') else None,
                    'atoms': [
                        {
                            'number': atom.number if hasattr(atom, 'number') else None,
                            'element': atom.element if hasattr(atom, 'element') else None,
                            'x': getattr(atom, 'x', None),
                            'y': getattr(atom, 'y', None),
                            'z': getattr(atom, 'z', None),
                            'charge': getattr(atom, 'charge', None)
                        }
                        for atom in compound_details.atoms
                    ] if hasattr(compound_details, 'atoms') and compound_details.atoms else None,
                    'bonds': [
                        {
                            'aid1': bond.aid1 if hasattr(bond, 'aid1') else None,
                            'aid2': bond.aid2 if hasattr(bond, 'aid2') else None,
                            'order': bond.order if hasattr(bond, 'order') else None,
                            'style': getattr(bond, 'style', None)
                        }
                        for bond in compound_details.bonds
                    ] if hasattr(compound_details, 'bonds') and compound_details.bonds else None,
                    'charge': compound_details.charge if hasattr(compound_details, 'charge') else None,
                    'synonyms': compound_details.synonyms if hasattr(compound_details, 'synonyms') else None
                }

                self.logger.info(f"[DEBUG] Molecule data: {molecule_data}")


                # Create a script agent with the specified model if provided
                script_agent = ScriptAgent(
                    llm_service=self.llm_service
                )

                # Generate script using molecule data with elemental indices
                script = script_agent.generate_script_from_molecule(compound['name'], user_query, molecule_data)

                # Validate and convert script if needed
                script = validate_and_convert_script(script)

                # Create the visualization HTML
                visualizer = MoleculeVisualizer()
                html = visualizer.generate_interactive_html(
                    pdb_data=pdb_data,
                    title=display_title,
                    script_data=script
                )

                if DEBUG_PUBCHEM:
                    write_debug_file('pubchem_debug.json', json.dumps({
                        'query': user_query,
                        'pdb_data_length': len(pdb_data),
                        'html_length': len(html),
                        'sdf_length': len(compound['sdf']),
                        'timestamp': datetime.datetime.now().isoformat()
                    }))

                return MoleculePackage(
                    pdb_data=pdb_data,
                    html=html,
                    title=display_title
                )

            except Exception as e:
                self.logger.error(f"[ERROR] Failed to generate visualization: {str(e)}")
                self.logger.error(f"[ERROR] Error type: {type(e).__name__}")
                import traceback
                self.logger.error(f"[ERROR] Traceback:\n{traceback.format_exc()}")
                raise ValueError(f"Failed to generate visualization: {str(e)}")

        except Exception as e:
            self.logger.error(f"[ERROR] Error in get_molecule_package: {str(e)}")
            self.logger.error(f"[ERROR] Error type: {type(e).__name__}")
            import traceback
            self.logger.error(f"[ERROR] Traceback:\n{traceback.format_exc()}")
            raise

        Args:
            user_query: The user's query about a molecule

        Returns:
            MoleculePackage: A package containing PDB data, HTML, and title for visualization

        Raises:
            ValueError: If no molecules are found or other errors occur
        """
        # Write a debug file to confirm this method is being called
        if DEBUG_PUBCHEM:
            debug_info = {
                "method": "get_molecule_package",
                "user_query": user_query,
                "timestamp": datetime.datetime.now().isoformat(),
            }
            write_debug_file(
                "pubchem_package_start.json", json.dumps(debug_info, indent=2)
            )

        self.logger.info(f"[DEBUG] Getting molecule package for query: {user_query}")

        try:
            # First, get the molecule SDFs
            search_result = self.get_molecule_sdfs(user_query)

            if not search_result:
                raise ValueError(
                    f"No molecules found for query: {user_query}, search result: {search_result}"
                )

            # Use the first compound
            compound = search_result[0]

            if compound["sdf"] is None:
                raise ValueError(
                    f"No SDF data available for compound {compound['name']} (CID: {compound['cid']})"
                )

            pdb_data = _sdf_to_pdb_block(compound["sdf"])

            # Get a display title (use name if available, otherwise ID)
            display_title = (
                compound["name"] if compound["name"] else f"CID {compound['cid']}"
            )
            self.logger.info(f"[DEBUG] Using molecule: {display_title}")

            if DEBUG_PUBCHEM:
                write_debug_file("pubchem_sdf.txt", compound["sdf"] or "")

            try:
                # Generate a display title for the molecule
                display_title = compound["name"] if compound["name"] else "Molecule"
                self.logger.info(f"Using display title: {display_title}")

                # Get compound details using PubChemPy for additional data
                compound_details = pcp.Compound.from_cid(compound["cid"])

                if compound_details.isomeric_smiles:
                    mol = Chem.MolFromSmiles(compound_details.isomeric_smiles)
                    if mol is None:
                        raise ValueError(
                            f"Failed to parse SMILES for compound {compound_details.name} (CID: {compound_details.cid})"
                        )

                    smarts_pattern = Chem.MolToSmarts(mol)

                    # NEW: Identify functional groups using RDKit fragment functions

                molecule_data = {
                    "name": compound["name"],
                    "cid": compound["cid"],
                    "smiles": (
                        compound_details.isomeric_smiles
                        if hasattr(compound_details, "isomeric_smiles")
                        else None
                    ),
                    "smarts_pattern": (
                        smarts_pattern if "smarts_pattern" in locals() else None
                    ),
                    "iupac_name": (
                        compound_details.iupac_name
                        if hasattr(compound_details, "iupac_name")
                        else None
                    ),
                    "molecular_formula": (
                        compound_details.molecular_formula
                        if hasattr(compound_details, "molecular_formula")
                        else None
                    ),
                    "molecular_weight": (
                        compound_details.molecular_weight
                        if hasattr(compound_details, "molecular_weight")
                        else None
                    ),
                    "elements": (
                        compound_details.elements
                        if hasattr(compound_details, "elements")
                        else None
                    ),
                    "atoms": (
                        [
                            {
                                "number": (
                                    atom.number if hasattr(atom, "number") else None
                                ),
                                "element": (
                                    atom.element if hasattr(atom, "element") else None
                                ),
                                "x": getattr(atom, "x", None),
                                "y": getattr(atom, "y", None),
                                "z": getattr(atom, "z", None),
                                "charge": getattr(atom, "charge", None),
                            }
                            for atom in compound_details.atoms
                        ]
                        if hasattr(compound_details, "atoms") and compound_details.atoms
                        else None
                    ),
                    "bonds": (
                        [
                            {
                                "aid1": bond.aid1 if hasattr(bond, "aid1") else None,
                                "aid2": bond.aid2 if hasattr(bond, "aid2") else None,
                                "order": bond.order if hasattr(bond, "order") else None,
                                "style": getattr(bond, "style", None),
                            }
                            for bond in compound_details.bonds
                        ]
                        if hasattr(compound_details, "bonds") and compound_details.bonds
                        else None
                    ),
                    "charge": (
                        compound_details.charge
                        if hasattr(compound_details, "charge")
                        else None
                    ),
                    "synonyms": (
                        compound_details.synonyms
                        if hasattr(compound_details, "synonyms")
                        else None
                    ),
                }

                self.logger.info(f"[DEBUG] Molecule data: {molecule_data}")

                # Create a script agent with the specified model if provided
                script_agent = ScriptAgent(llm_service=self.llm_service)

                # Generate script using molecule data with elemental indices
                script = script_agent.generate_script_from_molecule(
                    compound["name"], user_query, molecule_data
                )

                # Validate and convert script if needed
                script = validate_and_convert_script(script)

                # Create the visualization HTML
                visualizer = MoleculeVisualizer()
                html = visualizer.generate_interactive_html(
                    pdb_data=pdb_data, title=display_title, script_data=script
                )

                if DEBUG_PUBCHEM:
                    write_debug_file(
                        "pubchem_debug.json",
                        json.dumps(
                            {
                                "query": user_query,
                                "pdb_data_length": len(pdb_data),
                                "html_length": len(html),
                                "sdf_length": len(compound["sdf"]),
                                "timestamp": datetime.datetime.now().isoformat(),
                            }
                        ),
                    )

                return MoleculePackage(
                    pdb_data=pdb_data, html=html, title=display_title
                )

            except Exception as e:
                self.logger.error(f"[ERROR] Failed to generate visualization: {str(e)}")
                self.logger.error(f"[ERROR] Error type: {type(e).__name__}")
                import traceback

                self.logger.error(f"[ERROR] Traceback:\n{traceback.format_exc()}")
                raise ValueError(f"Failed to generate visualization: {str(e)}")

        except Exception as e:
            self.logger.error(f"[ERROR] Error in get_molecule_package: {str(e)}")
            self.logger.error(f"[ERROR] Error type: {type(e).__name__}")
            import traceback

            self.logger.error(f"[ERROR] Traceback:\n{traceback.format_exc()}")
            raise

    def visualize_molecule(self, user_query: str) -> MoleculePackage:
        """
        Create a visualization package for a molecule based on user query.

        Args:
            user_query: User's query about a molecule

        Returns:
            MoleculePackage containing PDB data and visualization HTML
        """
        try:
            # Get molecule data
            search_result = self.get_molecule_sdfs(user_query)

            if not search_result:
                raise ValueError(
                    f"No molecules found for query: {user_query}, search result: {search_result}"
                )

            # Use the first compound
            compound = search_result[0]

            if compound["sdf"] is None:
                raise ValueError(
                    f"No SDF data available for compound {compound['name']} (CID: {compound['cid']})"
                )

            pdb_data = _sdf_to_pdb_block(compound["sdf"])

            # Get a display title (use name if available, otherwise ID)
            display_title = (
                compound["name"] if compound["name"] else f"CID {compound['cid']}"
            )
            self.logger.info(f"[DEBUG] Using molecule: {display_title}")

            if DEBUG_PUBCHEM:
                write_debug_file("pubchem_sdf.txt", compound["sdf"] or "")

            try:
                # Generate a display title for the molecule
                display_title = compound["name"] if compound["name"] else "Molecule"
                self.logger.info(f"Using display title: {display_title}")

                # Get compound details using PubChemPy for additional data
                compound_details = pcp.Compound.from_cid(compound["cid"])

                if compound_details.isomeric_smiles:
                    mol = Chem.MolFromSmiles(compound_details.isomeric_smiles)
                    if mol is None:
                        raise ValueError(
                            f"Failed to parse SMILES for compound {compound_details.name} (CID: {compound_details.cid})"
                        )

                    smarts_pattern = Chem.MolToSmarts(mol)

                    # NEW: Identify functional groups using RDKit fragment functions

                molecule_data = {
                    "name": compound["name"],
                    "cid": compound["cid"],
                    "smiles": (
                        compound_details.isomeric_smiles
                        if hasattr(compound_details, "isomeric_smiles")
                        else None
                    ),
                    "smarts_pattern": (
                        smarts_pattern if "smarts_pattern" in locals() else None
                    ),
                    "iupac_name": (
                        compound_details.iupac_name
                        if hasattr(compound_details, "iupac_name")
                        else None
                    ),
                    "molecular_formula": (
                        compound_details.molecular_formula
                        if hasattr(compound_details, "molecular_formula")
                        else None
                    ),
                    "molecular_weight": (
                        compound_details.molecular_weight
                        if hasattr(compound_details, "molecular_weight")
                        else None
                    ),
                    "elements": (
                        compound_details.elements
                        if hasattr(compound_details, "elements")
                        else None
                    ),
                    "atoms": (
                        [
                            {
                                "number": (
                                    atom.number if hasattr(atom, "number") else None
                                ),
                                "element": (
                                    atom.element if hasattr(atom, "element") else None
                                ),
                                "x": getattr(atom, "x", None),
                                "y": getattr(atom, "y", None),
                                "z": getattr(atom, "z", None),
                                "charge": getattr(atom, "charge", None),
                            }
                            for atom in compound_details.atoms
                        ]
                        if hasattr(compound_details, "atoms") and compound_details.atoms
                        else None
                    ),
                    "bonds": (
                        [
                            {
                                "aid1": bond.aid1 if hasattr(bond, "aid1") else None,
                                "aid2": bond.aid2 if hasattr(bond, "aid2") else None,
                                "order": bond.order if hasattr(bond, "order") else None,
                                "style": getattr(bond, "style", None),
                            }
                            for bond in compound_details.bonds
                        ]
                        if hasattr(compound_details, "bonds") and compound_details.bonds
                        else None
                    ),
                    "charge": (
                        compound_details.charge
                        if hasattr(compound_details, "charge")
                        else None
                    ),
                    "synonyms": (
                        compound_details.synonyms
                        if hasattr(compound_details, "synonyms")
                        else None
                    ),
                }

                self.logger.info(f"[DEBUG] Molecule data: {molecule_data}")

                # Create a script agent with the specified model if provided
                script_agent = ScriptAgent(llm_service=self.llm_service)

                # Generate script using molecule data with elemental indices
                script = script_agent.generate_script_from_molecule(
                    compound["name"], user_query, molecule_data
                )

                # Validate and convert script if needed
                script = validate_and_convert_script(script)

                # Create the visualization HTML
                visualizer = MoleculeVisualizer()
                html = visualizer.generate_interactive_html(
                    pdb_data=pdb_data, title=display_title, script_data=script
                )

                if DEBUG_PUBCHEM:
                    write_debug_file(
                        "pubchem_debug.json",
                        json.dumps(
                            {
                                "query": user_query,
                                "pdb_data_length": len(pdb_data),
                                "html_length": len(html),
                                "sdf_length": len(compound["sdf"]),
                                "timestamp": datetime.datetime.now().isoformat(),
                            }
                        ),
                    )

                return MoleculePackage(
                    pdb_data=pdb_data, html=html, title=display_title
                )

            except Exception as e:
                self.logger.error(
                    f"Error in script generation or visualization: {str(e)}"
                )
                raise

        except Exception as e:
            self.logger.error(f"Error in molecule visualization: {str(e)}")
            raise

    def render_layout_placeholder(self, layout: List[Dict[str, Any]]) -> None:
        """Placeholder for future 2-D rendering implementation."""
        pass
