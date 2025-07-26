"""Modern PubChem client using generate_structured for LLM orchestration."""

import logging
from typing import Any, NamedTuple, Optional

from pydantic import BaseModel

from api.llm.openai_provider import generate_structured

from ..services import PubChemSearchService
from .molecule_assembler import MoleculeAssembler

logger = logging.getLogger(__name__)


class MoleculeInterpretation(BaseModel):
    """Response model for molecule name interpretation."""

    molecule_name: str
    confidence: str  # "high", "medium", "low", "none"


class MolecularFormula(BaseModel):
    """Response model for molecular formula extraction."""

    formula: str
    confidence: str  # "high", "medium", "low", "none"


class MoleculePackage(NamedTuple):
    """Container for molecule visualization package data."""

    pdb_data: str
    html: str
    title: str


class PubChemClient:
    """Modern PubChem client with proper separation of concerns and LLM orchestration."""

    def __init__(
        self,
        search_service: Optional[PubChemSearchService] = None,
        molecule_assembler: Optional[MoleculeAssembler] = None,
    ):
        """Initialize the PubChem client with dependency injection.

        Args:
            search_service: PubChem search service instance
            molecule_assembler: Molecule assembler instance
        """
        self.search_service = search_service or PubChemSearchService()
        self.molecule_assembler = molecule_assembler or MoleculeAssembler(self.search_service)

    def interpret_user_query(self, user_input: str) -> str:
        """Use LLM to interpret user input into a molecule name or identifier.

        Args:
            user_input: The user's query about a molecule

        Returns:
            A molecule name or identifier, or the original input if interpretation fails
        """
        try:
            system_prompt = """You are a chemistry professor that helps identify the simplest, most representative molecule names from user queries. Always prefer well-known, simple examples that clearly demonstrate the concept."""

            user_prompt = f"""Given this user input is related to molecular structures or trying to learn something about the microscopic structures of molecules or some topic related, give us the molecule that would best help the user learn what they are trying to learn about.

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

Respond with the molecule name and your confidence level."""

            response = generate_structured(
                user_prompt=user_prompt,
                response_model=MoleculeInterpretation,
                system_prompt=system_prompt,
                model_name="gpt-4o-mini",
            )

            if response.confidence == "none":
                logger.warning(f"LLM could not interpret query: {user_input}")
                return self._fallback_interpretation(user_input)

            return response.molecule_name.strip()

        except Exception as e:
            logger.warning(f"Failed to interpret user query with LLM: {str(e)}")
            return self._fallback_interpretation(user_input)

    def _fallback_interpretation(self, user_input: str) -> str:
        """Fallback interpretation when LLM fails.

        Args:
            user_input: Original user input

        Returns:
            Best guess molecule name or default
        """
        # Fallback 1: Try to use the input directly if it looks like a molecule name
        if len(user_input.split()) <= 3 and not any(char in user_input for char in "?!.,:;"):
            logger.info(f"Using user input directly as molecule name: '{user_input}'")
            return user_input

        # Fallback 2: Use a default molecule if the input is complex
        logger.info(f"Using default molecule 'water' for complex query: '{user_input}'")
        return "water"

    def get_molecular_formula(self, compound_name: str) -> Optional[str]:
        """Use LLM to get the molecular formula for a compound name.

        Args:
            compound_name: Name of the chemical compound

        Returns:
            Molecular formula if LLM is confident, None otherwise
        """
        try:
            system_prompt = "You are a chemistry expert. Provide molecular formulas for chemical compounds."

            user_prompt = f"""Given a chemical compound name, provide its molecular formula.
Only respond with the formula and confidence level. If you're not certain, set confidence to 'none'.

Example inputs and outputs:
Input: "water" → H2O (high confidence)
Input: "glucose" → C6H12O6 (high confidence)
Input: "random chemical" → unknown (no confidence)

Input: "{compound_name}"
Output:"""

            response = generate_structured(
                user_prompt=user_prompt,
                response_model=MolecularFormula,
                system_prompt=system_prompt,
                model_name="gpt-4o-mini",
            )

            if response.confidence == "none" or response.formula.lower() == "unknown":
                return None

            return response.formula

        except Exception as e:
            logger.error(f"Error getting molecular formula from LLM: {str(e)}")
            return None

    def get_molecule_data(self, user_query: str) -> dict[str, Any]:
        """Fetch raw molecule data without generating HTML.

        Args:
            user_query: User's query about a molecule

        Returns:
            Dictionary with molecule data including PDB, name, CID, formula, etc.

        Raises:
            ValueError: If no molecule data can be retrieved
        """
        logger.debug(f"Fetching molecule data for: {user_query}")

        # Interpret user query
        molecule_name = self.interpret_user_query(user_query)
        if not molecule_name:
            raise ValueError("Could not interpret user query into a molecule name.")

        # Search for the molecule
        compounds = self.search_service.search_with_fallbacks(molecule_name)
        if not compounds:
            # Try searching by formula if direct search fails
            formula = self.get_molecular_formula(molecule_name)
            if formula:
                logger.info(f"Trying formula search with: {formula}")
                compounds = self.search_service.search_by_formula(formula)

            if not compounds:
                raise ValueError(f"No compounds found for {molecule_name}")

        # Use the first compound and assemble data
        cid = compounds[0].cid
        return self.molecule_assembler.assemble_molecule_data(cid)

    def get_molecule_2d_info(self, user_query: str) -> dict[str, Any]:
        """Fetch 2D structural information for a molecule.

        Args:
            user_query: User's query about a molecule

        Returns:
            Dictionary with 2D molecular structure data

        Raises:
            ValueError: If no 2D data can be retrieved
        """
        logger.debug(f"Fetching 2D molecule info for: {user_query}")

        # Interpret user query
        molecule_name = self.interpret_user_query(user_query)
        if not molecule_name:
            raise ValueError("Could not interpret user query into a molecule name.")

        # Search for the molecule
        compounds = self.search_service.search_with_fallbacks(molecule_name)
        if not compounds:
            raise ValueError(f"No compounds found for {molecule_name}")

        # Use the first compound and assemble 2D data
        cid = compounds[0].cid
        return self.molecule_assembler.assemble_2d_molecule_data(cid)

    def get_molecules_2d_layout(self, queries: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Fetch 2D info for multiple molecules and attach layout boxes.

        Args:
            queries: List of dictionaries with 'query' and 'box' keys

        Returns:
            List of dictionaries with molecule data and box information
        """
        return self.molecule_assembler.assemble_molecules_2d_layout(queries)

    def get_molecule_sdfs(self, user_input: str) -> list[dict[str, Any]]:
        """Get SDF data for molecules based on user input.

        Args:
            user_input: User's query about a molecule

        Returns:
            List of dictionaries containing molecule data and SDF content
        """
        logger.info(f"Processing SDF request for: {user_input}")

        # Interpret the user's query
        molecule_name = self.interpret_user_query(user_input)
        if not molecule_name:
            logger.warning("Could not interpret user query into a molecule name")
            return []

        # Search for the molecule
        compounds = self.search_service.search_with_fallbacks(molecule_name)
        if not compounds:
            logger.warning(f"No compounds found for {molecule_name}")
            return []

        # Process compounds and get their SDF data
        results = []
        for compound in compounds:
            try:
                cid = compound.cid

                # Try to get 3D SDF first, fallback to 2D
                sdf_data = self.search_service.fetch_sdf(cid, "3d")
                if not sdf_data:
                    sdf_data = self.search_service.fetch_sdf(cid, "2d")

                if sdf_data:
                    results.append(
                        {
                            "name": getattr(compound, "iupac_name", str(cid)),
                            "cid": cid,
                            "formula": getattr(compound, "molecular_formula", ""),
                            "sdf": sdf_data,
                        }
                    )
                    logger.info(f"Successfully processed compound CID {cid}")
                else:
                    logger.error(f"Failed to get any SDF data for CID {cid}")

            except Exception as e:
                logger.error(f"Error processing compound: {str(e)}")
                continue

        return results
