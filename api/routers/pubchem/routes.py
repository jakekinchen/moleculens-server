"""
PubChem API routes for molecular structure retrieval and search.
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from typing import Optional

from agent_management.llm_service import LLMService
from agent_management.agent_factory import AgentFactory
from agent_management.models import PubChemSearchResult, PubChemCompound
from dependencies.use_llm import use_llm

router = APIRouter(
    prefix="/pubchem",
    tags=["pubchem"],
    responses={404: {"description": "Not found"}},
)

@router.get("/search/", response_model=PubChemSearchResult)
async def search_molecules(
    query: str = Query(..., description="Query string to search for molecules"),
    llm_service: LLMService = Depends(use_llm)
):
    """
    Search for molecules in PubChem based on a text query.
    The query will be interpreted by an LLM to extract molecule names.
    """
    try:
        agent = AgentFactory.create_pubchem_agent()
        return agent.get_molecule_sdfs(query)
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error searching PubChem: {str(e)}"
        )

@router.get("/compound/{cid}", response_model=Optional[PubChemCompound])
async def get_compound(
    cid: int = Query(..., description="PubChem Compound ID (CID)"),
    llm_service: LLMService = Depends(use_llm)
):
    """
    Get detailed information about a specific compound by its PubChem CID.
    """
    try:
        agent = AgentFactory.create_pubchem_agent()
        result = agent.get_compound_details(cid)
        if result is None:
            raise HTTPException(
                status_code=404,
                detail=f"Compound with CID {cid} not found"
            )
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error fetching compound details: {str(e)}"
        ) 