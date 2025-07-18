import requests

class RCSBAgent:
    """Simple agent to retrieve structure files from the RCSB PDB."""

    BASE_URL = "https://files.rcsb.org/download/"

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
        fmt = file_format.lower()
        if fmt not in {"pdb", "cif"}:
            raise ValueError("format must be 'pdb' or 'cif'")
        url = f"{self.BASE_URL}{identifier.upper()}.{ 'cif' if fmt == 'cif' else 'pdb' }"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.text
