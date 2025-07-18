import requests

class RCSBAgent:
    """Retrieve experimental and computed structures plus metadata from RCSB and AlphaFold."""

    BASE_URL = "https://files.rcsb.org/download/"
    AF_BASE_URL = "https://alphafold.ebi.ac.uk/files/"
    DATA_API_URL = "https://data.rcsb.org/rest/v1/core/"

    def _check_format(self, fmt: str) -> str:
        fmt = fmt.lower()
        if fmt not in {"pdb", "cif"}:
            raise ValueError("format must be 'pdb' or 'cif'")
        return fmt

    def fetch_structure(self, identifier: str, file_format: str = "pdb") -> str:
        """Download a structure file from RCSB.

        Parameters
        ----------
        identifier: PDB or UniProt identifier.
        file_format: "pdb" for PDB format or "cif" for mmCIF.

        Returns
        -------
        str
            File contents as a string.
        """
        fmt = self._check_format(file_format)
        extension = "cif" if fmt == "cif" else "pdb"
        url = f"{self.BASE_URL}{identifier.upper()}.{extension}"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.text

    def fetch_alphafold_model(self, uniprot_id: str, file_format: str = "pdb") -> str:
        """
        Retrieve a predicted structure from AlphaFold DB.
        uniprot_id: UniProt accession, e.g. P68871.
        file_format: 'pdb' or 'cif'.
        Returns the contents of the model file.
        """
        fmt = self._check_format(file_format)
        suffix = "model_v4.cif" if fmt == "cif" else "model_v4.pdb"
        filename = f"AF-{uniprot_id}-F1-{suffix}"
        resp = requests.get(f"{self.AF_BASE_URL}{filename}", timeout=30)
        resp.raise_for_status()
        return resp.text

    def fetch_entry_metadata(self, identifier: str) -> dict:
        """
        Retrieve JSON metadata for an entry via the RCSB Data API.
        identifier: PDB ID.
        Returns a dict with metadata.
        """
        resp = requests.get(f"{self.DATA_API_URL}entry/{identifier}", timeout=30)
        resp.raise_for_status()
        return resp.json()
