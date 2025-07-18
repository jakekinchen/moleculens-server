from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Literal, Dict, Any
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

class ModelRequest(BaseModel):
    uniprot_id: str
    format: Literal["pdb", "cif"] = "pdb"

class MetadataResponse(BaseModel):
    metadata: Dict[str, Any]

@router.post("/fetch-structure/", response_model=StructureResponse)
def fetch_structure(request: StructureRequest) -> StructureResponse:
    try:
        agent = RCSBAgent()
        content = agent.fetch_structure(request.identifier, request.format)
        return StructureResponse(data=content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/fetch-model/", response_model=StructureResponse)
def fetch_model(request: ModelRequest) -> StructureResponse:
    try:
        agent = RCSBAgent()
        content = agent.fetch_alphafold_model(request.uniprot_id, request.format)
        return StructureResponse(data=content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/entry/{identifier}", response_model=MetadataResponse)
def entry_metadata(identifier: str) -> MetadataResponse:
    try:
        agent = RCSBAgent()
        meta = agent.fetch_entry_metadata(identifier)
        return MetadataResponse(metadata=meta)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
