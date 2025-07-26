"""RCSB PDB and AlphaFold structure retrieval service."""

import logging
from typing import Optional

import requests

logger = logging.getLogger(__name__)


class RCSBService:
    """Service for retrieving structures from RCSB PDB and AlphaFold."""

    BASE_URL = "https://files.rcsb.org/download/"
    AF_BASE_URL = "https://alphafold.ebi.ac.uk/files/"
    DATA_API_URL = "https://data.rcsb.org/rest/v1/core/"

    def __init__(self, timeout: int = 30):
        """Initialize the RCSB service.

        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout

    def _check_format(self, fmt: str) -> str:
        """Validate and normalize file format.

        Args:
            fmt: File format string

        Returns:
            Normalized format string

        Raises:
            ValueError: If format is not supported
        """
        fmt = fmt.lower()
        if fmt not in {"pdb", "cif"}:
            raise ValueError("format must be 'pdb' or 'cif'")
        return fmt

    def fetch_structure(self, identifier: str, file_format: str = "pdb") -> Optional[str]:
        """Download a structure file from RCSB.

        Args:
            identifier: PDB identifier
            file_format: "pdb" for PDB format or "cif" for mmCIF

        Returns:
            File contents as string, or None if fetch fails
        """
        try:
            fmt = self._check_format(file_format)
            extension = "cif" if fmt == "cif" else "pdb"
            url = f"{self.BASE_URL}{identifier.upper()}.{extension}"

            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()
            return response.text

        except requests.exceptions.RequestException as e:
            logger.error(f"Network error fetching structure {identifier}: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Error fetching structure {identifier}: {str(e)}")
            return None

    def fetch_alphafold_model(self, uniprot_id: str, file_format: str = "pdb") -> Optional[str]:
        """Retrieve a predicted structure from AlphaFold DB.

        Args:
            uniprot_id: UniProt accession (e.g., P68871)
            file_format: "pdb" or "cif"

        Returns:
            Model file contents as string, or None if fetch fails
        """
        try:
            fmt = self._check_format(file_format)
            extension = "cif" if fmt == "cif" else "pdb"
            url = f"{self.AF_BASE_URL}AF-{uniprot_id}-F1-model_v4.{extension}"

            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()
            return response.text

        except requests.exceptions.RequestException as e:
            logger.error(f"Network error fetching AlphaFold model {uniprot_id}: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Error fetching AlphaFold model {uniprot_id}: {str(e)}")
            return None

    def get_structure_metadata(self, identifier: str) -> Optional[dict]:
        """Get metadata for a structure from RCSB API.

        Args:
            identifier: PDB identifier

        Returns:
            Structure metadata as dictionary, or None if fetch fails
        """
        try:
            url = f"{self.DATA_API_URL}entry/{identifier.upper()}"
            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()
            return response.json()

        except requests.exceptions.RequestException as e:
            logger.error(f"Network error fetching metadata for {identifier}: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Error fetching metadata for {identifier}: {str(e)}")
            return None
