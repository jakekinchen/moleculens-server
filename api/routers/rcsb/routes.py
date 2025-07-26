from typing import Any, Literal

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from api.pymol.services.rcsb import RCSBService

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
    metadata: dict[str, Any]


class AnnotationResponse(BaseModel):
    annotations: dict[str, Any]


class GraphQLRequest(BaseModel):
    identifier: str
    model_id: str


class UploadRequest(BaseModel):
    data: str
    filename: str = "upload.pdb"


class UploadResponse(BaseModel):
    upload_id: str


class AlignmentRequest(BaseModel):
    identifier1: str
    identifier2: str


@router.post("/fetch-structure/", response_model=StructureResponse)
def fetch_structure(request: StructureRequest) -> StructureResponse:
    try:
        service = RCSBService()
        content = service.fetch_structure(request.identifier, request.format)
        return StructureResponse(data=content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.post("/fetch-model/", response_model=StructureResponse)
def fetch_model(request: ModelRequest) -> StructureResponse:
    try:
        service = RCSBService()
        content = service.fetch_alphafold_model(request.uniprot_id, request.format)
        return StructureResponse(data=content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.get("/entry/{identifier}", response_model=MetadataResponse)
def entry_metadata(identifier: str) -> MetadataResponse:
    try:
        service = RCSBService()
        meta = service.fetch_entry_metadata(identifier)
        return MetadataResponse(metadata=meta)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.get("/annotations/{identifier}", response_model=AnnotationResponse)
def sequence_annotations(identifier: str) -> AnnotationResponse:
    try:
        service = RCSBService()
        data = service.fetch_sequence_annotations(identifier)
        return AnnotationResponse(annotations=data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.post("/computed-model/", response_model=MetadataResponse)
def computed_model(request: GraphQLRequest) -> MetadataResponse:
    try:
        service = RCSBService()
        data = service.fetch_graphql_model(request.identifier, request.model_id)
        return MetadataResponse(metadata=data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.post("/fetch-esm-model/", response_model=StructureResponse)
def fetch_esm_model(request: ModelRequest) -> StructureResponse:
    try:
        service = RCSBService()
        content = service.fetch_esmf_model(request.uniprot_id, request.format)
        return StructureResponse(data=content)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.post("/upload-structure/", response_model=UploadResponse)
def upload_structure(request: UploadRequest) -> UploadResponse:
    try:
        service = RCSBService()
        upload_id = service.upload_structure(request.data.encode(), request.filename)
        return UploadResponse(upload_id=upload_id)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.post("/align/", response_model=MetadataResponse)
def pairwise_alignment(request: AlignmentRequest) -> MetadataResponse:
    try:
        service = RCSBService()
        data = service.fetch_pairwise_alignment(request.identifier1, request.identifier2)
        return MetadataResponse(metadata=data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.get("/group/{group_id}", response_model=MetadataResponse)
def group_entries(group_id: str) -> MetadataResponse:
    try:
        service = RCSBService()
        data = service.fetch_group_entries(group_id)
        return MetadataResponse(metadata=data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@router.get("/feature-annotations/{identifier}", response_model=AnnotationResponse)
def feature_annotations(identifier: str) -> AnnotationResponse:
    try:
        service = RCSBService()
        data = service.fetch_feature_annotations(identifier)
        return AnnotationResponse(annotations=data)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e
