from fastapi import HTTPException

from api.pymol.services.rcsb import RCSBService

from .routes import router  # existing APIRouter instance


@router.get(
    "/sequence-coordinates/{identifier}",
    response_model=dict[str, list[float]],
    summary="Residue‑level Cartesian coordinates for a structure",
)
def get_sequence_coordinates(identifier: str):
    """Fetch residue‑to‑coordinate mapping from the RCSB
    Sequence Coordinates Service."""
    service = RCSBService()
    try:
        return service.fetch_sequence_coordinates(identifier)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
