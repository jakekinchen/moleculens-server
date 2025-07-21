"""Monkey‑patch supplying ``RCSBAgent.fetch_sequence_coordinates`` without
altering the original *rcsb_agent.py* file.

The RCSB Sequence Coordinates Service returns JSON of the form:

{
  "coordset": [
    {
      "seq_id": 1,
      "cartn_x": 12.345,
      "cartn_y": 34.567,
      "cartn_z": 78.901
    },
    ...
  ]
}

We collapse this into ``{str(seq_id): [x, y, z]}``.
"""

from typing import Dict, List

import requests

from .rcsb_agent import RCSBAgent


def _fetch_sequence_coordinates(
    self: "RCSBAgent", identifier: str
) -> Dict[str, List[float]]:
    """Retrieve residue‑level Cartesian coordinates.

    Parameters
    ----------
    identifier : str
        PDB ID, polymer‑entity ID or UniProt accession.

    Returns
    -------
    dict
        Mapping ``{residue_number: [x, y, z]}``.

    Raises
    ------
    ValueError
        If the service replies with an error or the payload is empty.
    """
    url = f"https://data.rcsb.org/rest/v1/sequence/coordinates/{identifier}"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code != 200:
            raise ValueError(
                f"Sequence Coordinates Service responded {resp.status_code} for '{identifier}'"
            )
        payload = resp.json()
    except Exception as exc:
        raise ValueError(
            f"Failed to contact Sequence Coordinates Service: {exc}"
        ) from exc

    coords: Dict[str, List[float]] = {}
    for entry in payload.get("coordset", []):
        seq_id = entry.get("seq_id")
        x, y, z = entry.get("cartn_x"), entry.get("cartn_y"), entry.get("cartn_z")
        if seq_id is None or None in (x, y, z):
            continue
        coords[str(seq_id)] = [float(x), float(y), float(z)]

    if not coords:
        raise ValueError(f"No coordinates returned for '{identifier}'.")
    return coords


# Attach method to class at import‑time
setattr(RCSBAgent, "fetch_sequence_coordinates", _fetch_sequence_coordinates)
