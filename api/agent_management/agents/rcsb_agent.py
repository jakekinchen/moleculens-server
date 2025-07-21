import requests


class RCSBAgent:
    """Retrieve experimental and computed structures plus metadata from RCSB
    and AlphaFold."""

    BASE_URL = "https://files.rcsb.org/download/"
    AF_BASE_URL = "https://alphafold.ebi.ac.uk/files/"
    DATA_API_URL = "https://data.rcsb.org/rest/v1/core/"
    SEQ_COORD_URL = "https://seq.rcsb.org/rest/v1/entry/{identifier}"
    GRAPHQL_URL = "https://data.rcsb.org/graphql"
    ESM_BASE_URL = "https://esmatlas.com/api/v1/prediction/"
    UPLOAD_URL = "https://www.rcsb.org/api/v1/molstar/upload"
    ALIGN_URL = "https://www.rcsb.org/api/v1/align"

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
        """Retrieve a predicted structure from AlphaFold DB.

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
        """Retrieve JSON metadata for an entry via the RCSB Data API.

        identifier: PDB ID.
        Returns a dict with metadata.
        """
        resp = requests.get(f"{self.DATA_API_URL}entry/{identifier}", timeout=30)
        resp.raise_for_status()
        return resp.json()

    def fetch_sequence_annotations(self, identifier: str) -> dict:
        """Retrieve residue-level annotations from the Sequence Coordinates
        Service."""
        url = self.SEQ_COORD_URL.format(identifier=identifier)
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.json()

    def fetch_graphql_model(self, identifier: str, model_id: str) -> dict:
        """Call the GraphQL API to obtain a computed structure model."""
        query = """
        query($id: String!, $model: String!) {
            entry(entry_id: $id) {
                rcsb_id
                struct {
                    title
                }
            }
        }
        """
        payload = {"query": query, "variables": {"id": identifier, "model": model_id}}
        resp = requests.post(self.GRAPHQL_URL, json=payload, timeout=30)
        resp.raise_for_status()
        return resp.json()

    def fetch_esmf_model(self, uniprot_id: str, file_format: str = "pdb") -> str:
        """Retrieve a predicted structure from the ESMFold API."""
        fmt = self._check_format(file_format)
        ext = "cif" if fmt == "cif" else "pdb"
        url = f"{self.ESM_BASE_URL}{uniprot_id}.{ext}"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.text

    def upload_structure(self, file_bytes: bytes, filename: str = "upload.pdb") -> str:
        """Upload a user-provided structure and return the shareable ID."""
        files = {"file": (filename, file_bytes)}
        resp = requests.post(self.UPLOAD_URL, files=files, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        return data.get("id") or data.get("url")

    def fetch_pairwise_alignment(self, id_one: str, id_two: str) -> dict:
        """Align two structures or accessions using the RCSB alignment API."""
        payload = {"id1": id_one, "id2": id_two}
        resp = requests.post(self.ALIGN_URL, json=payload, timeout=30)
        resp.raise_for_status()
        return resp.json()

    def fetch_group_entries(self, group_id: str) -> dict:
        """Retrieve entries that belong to a specific RCSB group."""
        url = f"{self.DATA_API_URL}group/{group_id}"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.json()

    def fetch_feature_annotations(self, identifier: str) -> dict:
        """Get domain and motif annotations for a structure or accession."""
        url = f"{self.DATA_API_URL}feature/{identifier}"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.json()
