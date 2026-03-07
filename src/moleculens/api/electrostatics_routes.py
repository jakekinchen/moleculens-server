"""FastAPI route handlers for electrostatics endpoints."""

import hashlib
import json
from typing import Any

from fastapi import APIRouter, HTTPException

from moleculens.api.schemas import (
    DipoleResponse,
    ElectrostaticsJobResponse,
    ElectrostaticsMeta,
    ElectrostaticsRequest,
    ElectrostaticsResult,
    ElectrostaticsTimings,
    ESPMeshResponse,
    PotentialRangeResponse,
    SurfaceResponse,
    UnitsResponse,
)
from moleculens.core import get_logger, settings
from moleculens.db import Job, JobQueue, JobStatus

logger = get_logger(__name__)

router = APIRouter(prefix="/v1/electrostatics", tags=["electrostatics"])
job_queue = JobQueue()


def _missing_electrostatics_result_error(cache_key: str) -> str:
    return f"Cached result missing for cache key {cache_key[:12]}"


def _load_or_invalidate_electrostatics_result(job: Job) -> ElectrostaticsResult | None:
    """Load a completed electrostatics result or invalidate the stale cache entry."""
    result = _load_electrostatics_result(job.cache_key, job.result_meta)
    if result is not None:
        return result

    error_message = _missing_electrostatics_result_error(job.cache_key)
    logger.warning(
        "Completed electrostatics job missing result payload; invalidating cache",
        job_id=job.id[:12],
        cache_key=job.cache_key[:12],
    )
    job_queue.invalidate_cache_key(job.cache_key, error_message)
    return None


def _compute_electrostatics_cache_key(
    sdf_content: str,
    charge: int,
    multiplicity: int,
    method: str,
    surface_params: dict[str, Any],
    potential_params: dict[str, Any],
    client_cache_key: str | None = None,
    inchi_key: str | None = None,
) -> str:
    """Compute stable cache key for electrostatics computation."""
    if client_cache_key:
        mol_id = client_cache_key
    elif inchi_key:
        mol_id = inchi_key
    else:
        mol_id = hashlib.sha256(sdf_content.strip().encode()).hexdigest()[:16]

    params = {
        "mol_id": mol_id,
        "type": "electrostatics",
        "charge": charge,
        "multiplicity": multiplicity,
        "method": method.lower(),
        "surface": surface_params,
        "potential": potential_params,
    }

    content = json.dumps(params, sort_keys=True)
    return hashlib.sha256(content.encode()).hexdigest()


@router.post("/compute", response_model=ElectrostaticsJobResponse)
async def compute_electrostatics(request: ElectrostaticsRequest) -> ElectrostaticsJobResponse:
    """Submit an electrostatics computation job.

    Computes:
    - Molecular surface mesh colored by electrostatic potential (ESP/MEP)
    - Dipole moment vector and magnitude

    If a cached result exists, returns immediately with cached=true.
    Otherwise, queues the job and returns job_id for polling.
    """
    # Build surface and potential params dicts
    surface_params = {
        "grid_spacing": request.surface.grid_spacing,
        "probe_radius": request.surface.probe_radius,
        "padding": request.surface.padding,
    }
    potential_params = {
        "units": request.potential.units,
        "softening_epsilon": request.potential.softening_epsilon,
        "clamp_percentiles": request.potential.clamp_percentiles,
    }

    # Compute cache key
    cache_key = _compute_electrostatics_cache_key(
        sdf_content=request.sdf_content,
        charge=request.charge,
        multiplicity=request.multiplicity,
        method=request.method,
        surface_params=surface_params,
        potential_params=potential_params,
        client_cache_key=request.client_cache_key,
        inchi_key=request.inchi_key,
    )

    try:
        # Submit job (reuses orbital job queue infrastructure)
        job_id, is_cached = job_queue.submit_job(
            sdf_content=request.sdf_content,
            method=request.method,
            basis="electrostatics",  # Use as job type marker
            grid_spacing=request.surface.grid_spacing,
            isovalue=0.0,  # Not used for electrostatics
            orbitals=["electrostatics"],  # Marker for job type
            inchi_key=request.inchi_key,
            geometry_hash=cache_key,  # Use our computed cache key
        )
    except Exception as e:
        logger.error("Failed to submit electrostatics job", error=str(e))
        raise HTTPException(status_code=400, detail=str(e)) from e

    # Get job to build response
    job = job_queue.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=500, detail="Job not found after creation")

    # Store additional params in request_params for worker
    # This is a workaround - ideally we'd have a separate table for electrostatics
    job_queue._update_job_params(
        job_id,
        {
            "job_type": "electrostatics",
            "sdf_content": request.sdf_content,
            "charge": request.charge,
            "multiplicity": request.multiplicity,
            "method": request.method,
            "surface": surface_params,
            "potential": potential_params,
        },
    )

    response = ElectrostaticsJobResponse(
        cached=is_cached,
        jobId=job.id,
        status=job.status.value,  # type: ignore[arg-type]
        cacheKey=cache_key,
        createdAt=job.created_at.isoformat() if job.created_at else None,
    )

    # If cached or done, load result
    if job.status == JobStatus.DONE:
        result = _load_or_invalidate_electrostatics_result(job)
        if result is None:
            retry_job_id, retry_is_cached = job_queue.submit_job(
                sdf_content=request.sdf_content,
                method=request.method,
                basis="electrostatics",
                grid_spacing=request.surface.grid_spacing,
                isovalue=0.0,
                orbitals=["electrostatics"],
                inchi_key=request.inchi_key,
                geometry_hash=cache_key,
            )
            retry_job = job_queue.get_job(retry_job_id)
            if retry_job is None:
                raise HTTPException(status_code=500, detail="Job not found after cache recovery")

            job_queue._update_job_params(
                retry_job_id,
                {
                    "job_type": "electrostatics",
                    "sdf_content": request.sdf_content,
                    "charge": request.charge,
                    "multiplicity": request.multiplicity,
                    "method": request.method,
                    "surface": surface_params,
                    "potential": potential_params,
                },
            )

            response = ElectrostaticsJobResponse(
                cached=retry_is_cached,
                jobId=retry_job.id,
                status=retry_job.status.value,  # type: ignore[arg-type]
                cacheKey=retry_job.cache_key,
                createdAt=retry_job.created_at.isoformat() if retry_job.created_at else None,
            )

            if retry_job.status == JobStatus.DONE:
                retry_result = _load_or_invalidate_electrostatics_result(retry_job)
                if retry_result is None:
                    return ElectrostaticsJobResponse(
                        cached=False,
                        jobId=retry_job.id,
                        status="error",
                        cacheKey=retry_job.cache_key,
                        createdAt=retry_job.created_at.isoformat()
                        if retry_job.created_at
                        else None,
                        computeTimeMs=retry_job.compute_time_ms,
                        errorMessage=_missing_electrostatics_result_error(retry_job.cache_key),
                    )

                response.result = retry_result
                response.compute_time_ms = retry_job.compute_time_ms

            return response

        response.result = result
        response.compute_time_ms = job.compute_time_ms

    return response


@router.get("/jobs/{job_id}", response_model=ElectrostaticsJobResponse)
async def get_electrostatics_job(job_id: str) -> ElectrostaticsJobResponse:
    """Get status and result of an electrostatics computation job."""
    job = job_queue.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    response = ElectrostaticsJobResponse(
        cached=False,
        jobId=job.id,
        status=job.status.value,  # type: ignore[arg-type]
        cacheKey=job.cache_key,
        createdAt=job.created_at.isoformat() if job.created_at else None,
        computeTimeMs=job.compute_time_ms,
    )

    if job.status == JobStatus.DONE:
        result = _load_or_invalidate_electrostatics_result(job)
        if result is None:
            return ElectrostaticsJobResponse(
                cached=False,
                jobId=job.id,
                status="error",
                cacheKey=job.cache_key,
                createdAt=job.created_at.isoformat() if job.created_at else None,
                computeTimeMs=job.compute_time_ms,
                errorMessage=_missing_electrostatics_result_error(job.cache_key),
            )

        response.result = result
    elif job.status == JobStatus.ERROR:
        response.error_message = job.error_message

    return response


def _load_electrostatics_result(
    cache_key: str, result_meta: dict | None
) -> ElectrostaticsResult | None:
    """Load electrostatics result from cache."""
    if result_meta is None:
        return None

    cache_dir = settings.molecule_cache_dir / cache_key
    result_file = cache_dir / "electrostatics_result.json"

    if not result_file.exists():
        logger.warning("Electrostatics result file not found", cache_key=cache_key[:12])
        return None

    try:
        with result_file.open() as f:
            data = json.load(f)
    except Exception as e:
        logger.error("Failed to load electrostatics result", cache_key=cache_key[:12], error=str(e))
        return None

    # Build response from cached data
    surface_data = data.get("surface", {})
    mesh_data = surface_data.get("mesh", {})
    range_data = surface_data.get("range", {})
    units_data = surface_data.get("units", {})

    surface = SurfaceResponse(
        mesh=ESPMeshResponse(
            vertices=mesh_data.get("vertices", ""),
            normals=mesh_data.get("normals", ""),
            indices=mesh_data.get("indices", ""),
            potential=mesh_data.get("potential", ""),
            vertexCount=mesh_data.get("vertexCount", 0),
            triangleCount=mesh_data.get("triangleCount", 0),
        ),
        range=PotentialRangeResponse(
            min=range_data.get("min", 0),
            max=range_data.get("max", 0),
            p05=range_data.get("p05"),
            p95=range_data.get("p95"),
        ),
        units=UnitsResponse(
            coords=units_data.get("coords", "angstrom"),
            potential=units_data.get("potential", "kcal/mol/e"),
        ),
    )

    dipole_data = data.get("dipole", {})
    dipole = DipoleResponse(
        originAngstrom=dipole_data.get("origin_angstrom", [0, 0, 0]),
        vectorDebye=dipole_data.get("vector_debye", [0, 0, 0]),
        magnitudeDebye=dipole_data.get("magnitude_debye", 0),
        convention=dipole_data.get("convention", "physics"),
    )

    meta_data = data.get("meta", {})
    timings = meta_data.get("timings_ms", {})
    meta = ElectrostaticsMeta(
        chargeModel=meta_data.get("charge_model", "unknown"),
        method=meta_data.get("method", "unknown"),
        timingsMs=ElectrostaticsTimings(
            total=timings.get("total", 0),
            surface=timings.get("surface"),
            charges=timings.get("charges"),
            potential=timings.get("potential"),
        ),
        sdfHash=meta_data.get("sdf_hash", ""),
        molecularCharge=meta_data.get("molecular_charge", 0),
        multiplicity=meta_data.get("multiplicity", 1),
    )

    return ElectrostaticsResult(
        surface=surface,
        dipole=dipole,
        meta=meta,
    )
