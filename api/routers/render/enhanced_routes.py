"""Enhanced render routes for client-friendly 2D/3D data."""

import base64
from typing import Literal, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from .routes import RenderRequest
from .routes import render as base_render

router = APIRouter(prefix="/render", tags=["Enhanced Render"])


class MoleculeRenderRequest(BaseModel):
    """Simplified request for molecule rendering."""

    molecule_name: str
    render_type: Literal["2d_transparent", "3d_pdb", "both"] = "both"
    size: Literal["small", "medium", "large"] = "medium"
    quality: Literal["fast", "high", "publication"] = "high"


class MoleculeRenderResponse(BaseModel):
    """Response with both 2D and 3D data."""

    molecule_name: str
    png_base64: Optional[str] = None
    png_url: Optional[str] = None
    pdb_data: Optional[str] = None
    pdb_url: Optional[str] = None
    metadata: dict = {}


@router.post("/molecule", response_model=MoleculeRenderResponse)
async def render_molecule(request: MoleculeRenderRequest) -> MoleculeRenderResponse:
    """Render molecule in 2D transparent PNG and/or 3D PDB format."""

    # Size configurations
    size_configs = {"small": (256, 256), "medium": (512, 512), "large": (1024, 1024)}

    # Quality configurations
    quality_configs = {
        "fast": {"ray_trace": False, "dpi": 72, "antialias": False},
        "high": {"ray_trace": True, "dpi": 150, "antialias": True},
        "publication": {"ray_trace": True, "dpi": 300, "antialias": True, "ray_trace_mode": "poster"},
    }

    resolution = size_configs[request.size]
    quality_opts = quality_configs[request.quality]

    response_data = MoleculeRenderResponse(molecule_name=request.molecule_name)

    try:
        # Get 2D transparent PNG
        if request.render_type in ["2d_transparent", "both"]:
            png_request = RenderRequest(
                description=f"Show {request.molecule_name} molecule with transparent background",
                format="image",
                transparent_background=True,
                resolution=resolution,
                **quality_opts,
            )

            png_response = await base_render(png_request)

            # Handle response format
            if hasattr(png_response, "body"):
                # Direct file response
                response_data.png_base64 = base64.b64encode(png_response.body).decode()
            elif isinstance(png_response, dict) and "url" in png_response:
                # Large file URL
                response_data.png_url = png_response["url"]
                response_data.metadata.update(png_response.get("metadata", {}))

        # Get 3D PDB data
        if request.render_type in ["3d_pdb", "both"]:
            pdb_request = RenderRequest(description=f"Load {request.molecule_name} molecule structure", format="model")

            pdb_response = await base_render(pdb_request)

            # Handle response format
            if hasattr(pdb_response, "body"):
                # Direct file response
                response_data.pdb_data = pdb_response.body.decode()
            elif isinstance(pdb_response, dict) and "url" in pdb_response:
                # Large file URL
                response_data.pdb_url = pdb_response["url"]
                response_data.metadata.update(pdb_response.get("metadata", {}))

        return response_data

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to render molecule: {str(e)}") from e


class BatchMoleculeRequest(BaseModel):
    """Request for batch molecule rendering."""

    molecules: list[str]
    render_type: Literal["2d_transparent", "3d_pdb", "both"] = "2d_transparent"
    size: Literal["small", "medium", "large"] = "small"
    quality: Literal["fast", "high"] = "fast"


@router.post("/batch", response_model=list[MoleculeRenderResponse])
async def render_molecules_batch(request: BatchMoleculeRequest) -> list[MoleculeRenderResponse]:
    """Render multiple molecules efficiently."""

    results = []
    for molecule in request.molecules:
        try:
            molecule_request = MoleculeRenderRequest(
                molecule_name=molecule, render_type=request.render_type, size=request.size, quality=request.quality
            )
            result = await render_molecule(molecule_request)
            results.append(result)
        except Exception as e:
            # Add failed result
            results.append(MoleculeRenderResponse(molecule_name=molecule, metadata={"error": str(e)}))

    return results


class ProteinRenderRequest(BaseModel):
    """Request for protein structure rendering."""

    pdb_id: str
    render_type: Literal["2d_transparent", "3d_pdb", "both"] = "both"
    representation: Literal["cartoon", "surface", "stick", "ribbon"] = "cartoon"
    size: Literal["small", "medium", "large"] = "medium"
    quality: Literal["fast", "high", "publication"] = "high"


@router.post("/protein", response_model=MoleculeRenderResponse)
async def render_protein(request: ProteinRenderRequest) -> MoleculeRenderResponse:
    """Render protein structure with specific representation."""

    size_configs = {"small": (256, 256), "medium": (512, 512), "large": (1024, 1024)}

    quality_configs = {
        "fast": {"ray_trace": False, "dpi": 72},
        "high": {"ray_trace": True, "dpi": 150},
        "publication": {"ray_trace": True, "dpi": 300, "ray_trace_mode": "poster"},
    }

    resolution = size_configs[request.size]
    quality_opts = quality_configs[request.quality]

    response_data = MoleculeRenderResponse(molecule_name=request.pdb_id)

    try:
        # Get 2D representation
        if request.render_type in ["2d_transparent", "both"]:
            description = f"Load PDB {request.pdb_id} and show with {request.representation} representation, transparent background"

            png_request = RenderRequest(
                description=description,
                format="image",
                transparent_background=True,
                resolution=resolution,
                **quality_opts,
            )

            png_response = await base_render(png_request)

            if hasattr(png_response, "body"):
                response_data.png_base64 = base64.b64encode(png_response.body).decode()
            elif isinstance(png_response, dict) and "url" in png_response:
                response_data.png_url = png_response["url"]

        # Get 3D PDB data
        if request.render_type in ["3d_pdb", "both"]:
            # Use RCSB endpoint for known structures
            from ..rcsb.routes import StructureRequest, fetch_structure

            structure_request = StructureRequest(identifier=request.pdb_id, format="pdb")

            structure_response = fetch_structure(structure_request)
            response_data.pdb_data = structure_response.data

        return response_data

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to render protein: {str(e)}") from e
