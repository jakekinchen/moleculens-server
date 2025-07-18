from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Literal
from agent_management.agents.rcsb_agent import RCSBAgent

router = APIRouter(
    prefix="/rcsb",
    tags=["RCSB"],
    responses={404: {"description": "Not found"}},
)

class StructureRequest(BaseModel):
    identifier: str
    format: Literal["pdb", "cif"] = "pdb"

class StructureResponse(BaseModel):
    data: str

@router.post("/fetch-structure/", response_model=StructureResponse)
def fetch_structure(request: StructureRequest) -> StructureResponse:
    try:
        agent = RCSBAgent()
        content = agent.fetch_structure(request.identifier, request.format)
        return StructureResponse(data=content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
