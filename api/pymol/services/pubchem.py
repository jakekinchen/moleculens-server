"""PubChem search and data retrieval service."""

import logging
import re
import urllib.parse
from typing import Optional

import pubchempy as pcp
import requests

logger = logging.getLogger(__name__)


class PubChemSearchService:
    """Service for searching and retrieving data from PubChem."""

    def __init__(self, timeout: int = 30):
        """Initialize the PubChem search service.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout

    def normalize_query(self, query: str) -> list[str]:
        """Normalize a chemical query string to improve search results.

        Args:
            query: Raw query string

        Returns:
            List of normalized query strings to try
        """
        query = query.strip()
        variations = [query]

        # Handle special cases for chemical notation
        if "[" in query and "]" in query:
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
                        "4,7,13,16,21,24-Hexaoxa-1,10-diazabicyclo[8.8.8]hexacosane",
                    ]
                )

        # Handle complex names with multiple parts
        if " and " in query.lower() or " with " in query.lower():
            parts = re.split(r"\s+(?:and|with)\s+", query, flags=re.IGNORECASE)
            if parts:
                variations.append(parts[0].strip())

        # Handle specific chemical classes
        if "phosphazene" in query.lower():
            variations.extend(
                [
                    "Hexachlorocyclotriphosphazene",
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

    def search_by_name_rest(self, query: str, max_results: int = 5) -> list[pcp.Compound]:
        """Search PubChem using the REST API directly.

        Args:
            query: Search query
            max_results: Maximum number of results to return

        Returns:
            List of PubChem Compound objects
        """
        logger.debug(f"Searching PubChem REST API for: {query}")

        encoded_query = urllib.parse.quote(query)
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_query}/cids/JSON"

        try:
            response = requests.get(search_url, timeout=self.timeout)
            if response.status_code != 200:
                logger.warning(f"PubChem REST search failed with status {response.status_code}")
                return []

            data = response.json()
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
                    compound = pcp.Compound.from_cid(cid)
                    compounds.append(compound)
                except Exception as e:
                    logger.warning(f"Failed to get details for CID {cid}: {str(e)}")
                    continue

            return compounds

        except requests.exceptions.RequestException as e:
            logger.error(f"Network error in PubChem REST search: {str(e)}")
            return []

    def search_by_name_direct(self, query: str) -> list[pcp.Compound]:
        """Search PubChem using direct compound search.

        Args:
            query: Search query

        Returns:
            List of PubChem Compound objects
        """
        logger.debug(f"Direct PubChem search for: {query}")

        try:
            # First try exact name match
            search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(query)}/JSON"
            response = requests.get(search_url, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                if "PC_Compounds" in data:
                    cid = data["PC_Compounds"][0]["id"]["id"]["cid"]
                    logger.debug(f"Found exact match CID: {cid}")
                    compound = pcp.Compound.from_cid(int(cid))
                    return [compound]

            # If exact match fails, try the autocomplete API
            autocomplete_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{urllib.parse.quote(query)}/json"
            )
            response = requests.get(autocomplete_url, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                if "dictionary_terms" in data and "compound" in data["dictionary_terms"]:
                    for suggestion in data["dictionary_terms"]["compound"]:
                        if "cid" in suggestion:
                            cid = suggestion["cid"]
                            logger.debug(f"Found suggested CID: {cid}")
                            compound = pcp.Compound.from_cid(int(cid))
                            return [compound]

            return []

        except requests.exceptions.RequestException as e:
            logger.error(f"Network error in direct PubChem search: {str(e)}")
            return []

    def search_by_formula(self, formula: str) -> list[pcp.Compound]:
        """Search PubChem by molecular formula.

        Args:
            formula: Molecular formula

        Returns:
            List of PubChem Compound objects
        """
        try:
            compounds = pcp.get_compounds(formula, "formula")
            return compounds if compounds else []
        except Exception as e:
            logger.error(f"Error in formula search for {formula}: {str(e)}")
            return []

    def search_with_fallbacks(self, query: str) -> list[pcp.Compound]:
        """Search for a compound using multiple methods with fallbacks.

        Args:
            query: The compound name or identifier to search for

        Returns:
            List of PubChem Compound objects
        """
        logger.info(f"Starting search with fallbacks for: {query}")

        # Step 1: Try direct search with normalized queries
        normalized_queries = self.normalize_query(query)
        for normalized_query in normalized_queries:
            try:
                results = self.search_by_name_direct(normalized_query)
                if results:
                    logger.info(f"Found results using direct search with query: {normalized_query}")
                    return results
            except Exception as e:
                logger.warning(f"Direct search failed for {normalized_query}: {str(e)}")

        # Step 2: Try REST search with normalized queries
        for normalized_query in normalized_queries:
            try:
                results = self.search_by_name_rest(normalized_query)
                if results:
                    logger.info(f"Found results using REST search with query: {normalized_query}")
                    return results
            except Exception as e:
                logger.warning(f"REST search failed for {normalized_query}: {str(e)}")

        logger.warning(f"All search methods failed for query: {query}")
        return []

    def fetch_sdf(self, cid: int, record_type: str = "3d") -> Optional[str]:
        """Fetch SDF data for a given PubChem CID.

        Args:
            cid: PubChem Compound ID
            record_type: "3d" or "2d" for the type of structure

        Returns:
            SDF data as string, or None if fetch fails
        """
        if record_type == "3d":
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        else:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"

        try:
            response = requests.get(url, timeout=self.timeout)
            if response.status_code == 200:
                return response.text
            else:
                logger.warning(f"Failed to fetch {record_type} SDF for CID {cid} (status {response.status_code})")
                return None
        except requests.exceptions.RequestException as e:
            logger.error(f"Network error fetching {record_type} SDF for CID {cid}: {str(e)}")
            return None

    def get_compound_details(self, cid: int) -> Optional[pcp.Compound]:
        """Get detailed information for a specific compound by CID.

        Args:
            cid: PubChem Compound ID

        Returns:
            PubChem Compound object or None if fetch fails
        """
        try:
            return pcp.Compound.from_cid(cid)
        except Exception as e:
            logger.error(f"Error fetching compound details for CID {cid}: {str(e)}")
            return None
