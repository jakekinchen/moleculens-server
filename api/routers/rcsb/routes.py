from typing import Any, Dict, Literal

from agent_management.agents.rcsb_agent import RCSBAgent
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

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


class AnnotationResponse(BaseModel):
    annotations: Dict[str, Any]


class GraphQLRequest(BaseModel):
    identifier: str
    model_id: str


class UploadRequest(BaseModel):
    data: str
    filename: str = "upload.pdb"


class UploadResponse(BaseModel):
    upload_id: str


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


@router.get("/annotations/{identifier}", response_model=AnnotationResponse)
def sequence_annotations(identifier: str) -> AnnotationResponse:
    try:
        agent = RCSBAgent()
        data = agent.fetch_sequence_annotations(identifier)
        return AnnotationResponse(annotations=data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/computed-model/", response_model=MetadataResponse)
def computed_model(request: GraphQLRequest) -> MetadataResponse:
    try:
        agent = RCSBAgent()
        data = agent.fetch_graphql_model(request.identifier, request.model_id)
        return MetadataResponse(metadata=data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/fetch-esm-model/", response_model=StructureResponse)
def fetch_esm_model(request: ModelRequest) -> StructureResponse:
    try:
        agent = RCSBAgent()
        content = agent.fetch_esmf_model(request.uniprot_id, request.format)
        return StructureResponse(data=content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/upload-structure/", response_model=UploadResponse)
def upload_structure(request: UploadRequest) -> UploadResponse:
    try:
        agent = RCSBAgent()
        upload_id = agent.upload_structure(request.data.encode(), request.filename)
        return UploadResponse(upload_id=upload_id)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
