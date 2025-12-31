"""FastAPI route handlers."""

import json

from fastapi import APIRouter, HTTPException, Response
from fastapi.responses import FileResponse

from moleculens import __version__
from moleculens.api.schemas import (
    ComputeRequest,
    ComputeResult,
    ComputeResultMeta,
    HealthResponse,
    JobResponse,
    MeshDataResponse,
    OrbitalResponse,
)
from moleculens.core import get_logger, settings
from moleculens.db import JobQueue, JobStatus, get_db_engine

logger = get_logger(__name__)

router = APIRouter()
job_queue = JobQueue()


@router.get("/health", response_model=HealthResponse)
async def health_check() -> HealthResponse:
    """Check service health."""
    # Check database connectivity
    db_ok = False
    try:
        engine = get_db_engine()
        with engine.connect() as conn:
            conn.execute(__import__("sqlalchemy").text("SELECT 1"))
        db_ok = True
    except Exception as e:
        logger.warning("Database health check failed", error=str(e))

    # Check cache directory
    cache_ok = settings.molecule_cache_dir.exists()

    status = "healthy"
    if not db_ok:
        status = "unhealthy"
    elif not cache_ok:
        status = "degraded"

    return HealthResponse(
        status=status,  # type: ignore[arg-type]
        version=__version__,
        database=db_ok,
        cacheDir=cache_ok,
    )


@router.post("/v1/orbitals/compute", response_model=JobResponse)
async def compute_orbitals(request: ComputeRequest) -> JobResponse:
    """Submit an orbital computation job.

    If a cached result exists, returns immediately with cached=true.
    Otherwise, queues the job and returns job_id for polling.
    """
    try:
        job_id, is_cached = job_queue.submit_job(
            sdf_content=request.sdf_content,
            method=request.method,
            basis=request.basis,
            grid_spacing=request.grid_spacing,
            isovalue=request.isovalue,
            orbitals=request.orbitals,
            inchi_key=request.inchi_key,
        )
    except Exception as e:
        logger.error("Failed to submit job", error=str(e))
        raise HTTPException(status_code=400, detail=str(e)) from e

    # Get job to build response
    job = job_queue.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=500, detail="Job not found after creation")

    response = JobResponse(
        cached=is_cached,
        jobId=job.id,
        status=job.status.value,  # type: ignore[arg-type]
        cacheKey=job.cache_key,
        createdAt=job.created_at.isoformat() if job.created_at else None,
    )

    # If cached or done, load result
    if job.status == JobStatus.DONE:
        response.result = _load_result(job.cache_key, job.result_meta)
        response.compute_time_ms = job.compute_time_ms

    return response


@router.get("/v1/orbitals/jobs/{job_id}", response_model=JobResponse)
async def get_job_status(job_id: str) -> JobResponse:
    """Get status and result of a computation job."""
    job = job_queue.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    response = JobResponse(
        cached=False,
        jobId=job.id,
        status=job.status.value,  # type: ignore[arg-type]
        cacheKey=job.cache_key,
        createdAt=job.created_at.isoformat() if job.created_at else None,
        computeTimeMs=job.compute_time_ms,
    )

    if job.status == JobStatus.DONE:
        response.result = _load_result(job.cache_key, job.result_meta)
    elif job.status == JobStatus.ERROR:
        response.error_message = job.error_message

    return response


@router.get("/v1/orbitals/artifacts/{cache_key}/{filename}")
async def get_artifact(cache_key: str, filename: str) -> Response:
    """Download a cached artifact file.

    Useful for downloading raw cube files or other artifacts.
    """
    # Validate filename to prevent path traversal
    if ".." in filename or "/" in filename or "\\" in filename:
        raise HTTPException(status_code=400, detail="Invalid filename")

    artifact_path = job_queue.get_artifact_path(cache_key, filename)
    if artifact_path is None:
        raise HTTPException(status_code=404, detail="Artifact not found")

    # Determine content type
    content_type = "application/octet-stream"
    if filename.endswith(".json"):
        content_type = "application/json"
    elif filename.endswith(".cube"):
        content_type = "chemical/x-gaussian-cube"

    return FileResponse(
        path=artifact_path,
        media_type=content_type,
        filename=filename,
    )


def _load_result(cache_key: str, result_meta: dict | None) -> ComputeResult | None:
    """Load computation result from cache."""
    if result_meta is None:
        return None

    cache_dir = settings.molecule_cache_dir / cache_key
    result_file = cache_dir / "result.json"

    if not result_file.exists():
        logger.warning("Result file not found", cache_key=cache_key[:12])
        return None

    try:
        with result_file.open() as f:
            data = json.load(f)
    except Exception as e:
        logger.error("Failed to load result", cache_key=cache_key[:12], error=str(e))
        return None

    # Build response from cached data
    orbitals = {}
    for name, orb_data in data.get("orbitals", {}).items():
        positive = None
        negative = None

        if orb_data.get("positive"):
            positive = MeshDataResponse(
                vertices=orb_data["positive"]["vertices"],
                normals=orb_data["positive"]["normals"],
                indices=orb_data["positive"]["indices"],
                vertexCount=orb_data["positive"]["vertexCount"],
                triangleCount=orb_data["positive"]["triangleCount"],
            )
        if orb_data.get("negative"):
            negative = MeshDataResponse(
                vertices=orb_data["negative"]["vertices"],
                normals=orb_data["negative"]["normals"],
                indices=orb_data["negative"]["indices"],
                vertexCount=orb_data["negative"]["vertexCount"],
                triangleCount=orb_data["negative"]["triangleCount"],
            )

        orbitals[name] = OrbitalResponse(
            positive=positive,
            negative=negative,
            energyEv=orb_data.get("energyEv"),
            isovalue=orb_data.get("isovalue", 0.05),
        )

    density = None
    if data.get("density"):
        density = MeshDataResponse(
            vertices=data["density"]["vertices"],
            normals=data["density"]["normals"],
            indices=data["density"]["indices"],
            vertexCount=data["density"]["vertexCount"],
            triangleCount=data["density"]["triangleCount"],
        )

    meta = data.get("meta", {})

    return ComputeResult(
        orbitals=orbitals,
        density=density,
        meta=ComputeResultMeta(
            method=meta.get("method", "unknown"),
            basis=meta.get("basis", "unknown"),
            gridSpacingAngstrom=meta.get("gridSpacingAngstrom", 0.25),
            psi4NumThreads=meta.get("psi4NumThreads", 1),
            computeTimeMs=meta.get("computeTimeMs"),
        ),
    )
