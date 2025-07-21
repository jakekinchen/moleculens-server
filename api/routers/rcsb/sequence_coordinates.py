from typing import Dict, List

from fastapi import HTTPException

from api.agent_management.agent_factory import AgentFactory

# Ensure the monkey‑patch is loaded
from api.agent_management.agents import rcsb_sequence_patch  # noqa: F401

from .routes import router  # existing APIRouter instance


@router.get(
    "/sequence-coordinates/{identifier}",
    response_model=Dict[str, List[float]],
    summary="Residue‑level Cartesian coordinates for a structure",
)
def get_sequence_coordinates(identifier: str):
    """
    Fetch residue‑to‑coordinate mapping from the RCSB Sequence Coordinates Service.
    """
    agent = AgentFactory.create_rcsb_agent()
    try:
        return agent.fetch_sequence_coordinates(identifier)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc))
