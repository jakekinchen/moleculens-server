"""FastAPI route handlers for conformer endpoints."""

import asyncio
import hashlib
import json
import time
from typing import Any

from fastapi import APIRouter, HTTPException, Response
from fastapi.responses import FileResponse

from moleculens.api.schemas import (
    ConformerJobResponse,
    ConformerMeta,
    ConformerRequest,
    ConformerTimingMs,
)
from moleculens.core import get_logger, settings
from moleculens.db import JobQueue, JobStatus

logger = get_logger(__name__)

router = APIRouter(prefix="/v1/conformers", tags=["conformers"])
job_queue = JobQueue()


def _compute_conformer_cache_key(
    molblock2d: str,
    geom_version: str,
    params: dict[str, Any],
) -> str:
    payload = (
        molblock2d.encode()
        + b"\n"
        + geom_version.encode()
        + b"\n"
        + json.dumps(params, sort_keys=True, separators=(",", ":")).encode()
    )
    return hashlib.sha256(payload).hexdigest()


def _job_status_to_response(job_status: JobStatus) -> str:
    if job_status == JobStatus.DONE:
        return "done"
    if job_status == JobStatus.ERROR:
        return "failed"
    return "pending"


def _load_conformer_result(
    cache_key: str,
    result_meta: dict[str, Any] | None,
) -> tuple[str, ConformerMeta | None, ConformerTimingMs | None] | None:
    if result_meta is None:
        return None

    cache_dir = settings.molecule_cache_dir / cache_key
    result_file = cache_dir / "structure.sdf"

    if not result_file.exists():
        logger.warning("Conformer result file not found", cache_key=cache_key[:12])
        return None

    sdf3d = result_file.read_text()

    timings_data = result_meta.get("timings_ms", {})
    timing = None
    if timings_data:
        timing = ConformerTimingMs(
            embed=timings_data.get("embed"),
            opt=timings_data.get("opt"),
            total=timings_data.get("total"),
        )

    meta = ConformerMeta(
        method=result_meta.get("method", "unknown"),
        opt=result_meta.get("opt", "none"),
        max_attempts=result_meta.get("max_attempts", 0),
        max_opt_iters=result_meta.get("max_opt_iters", 0),
        add_hs=result_meta.get("add_hs", False),
        geom_version=result_meta.get("geom_version", "unknown"),
        quality=result_meta.get("quality", "unknown"),
        has_metal=result_meta.get("has_metal", False),
        embed_status=result_meta.get("embed_status", "unknown"),
        opt_status=result_meta.get("opt_status", "unknown"),
    )

    return sdf3d, meta, timing


def _build_conformer_response(job: Any, cache_key: str) -> ConformerJobResponse:
    response = ConformerJobResponse(
        status=_job_status_to_response(job.status),
        cache_key=cache_key,
        job_id=job.id,
    )

    if job.status == JobStatus.DONE:
        loaded = _load_conformer_result(cache_key, job.result_meta)
        if loaded:
            sdf3d, meta, timing = loaded
            response.sdf3d = sdf3d
            response.meta = meta
            response.timing_ms = timing
        else:
            response.status = "failed"
            response.error_message = "Conformer result missing"
    elif job.status == JobStatus.ERROR:
        response.error_message = job.error_message

    return response


def _build_conformer_params(request: ConformerRequest) -> dict[str, Any]:
    return {
        "method": request.params.method,
        "opt": request.params.opt,
        "max_attempts": request.params.max_attempts,
        "max_opt_iters": request.params.max_opt_iters,
        "add_hs": request.params.add_hs,
    }


def _job_is_terminal(job: Any) -> bool:
    return job.status in (JobStatus.DONE, JobStatus.ERROR)


async def _wait_for_job(job_id: str, wait_ms: int) -> None:
    if wait_ms <= 0:
        return

    deadline = time.monotonic() + (wait_ms / 1000.0)
    while time.monotonic() < deadline:
        job = job_queue.get_job(job_id)
        if job is None or _job_is_terminal(job):
            return

        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return

        await asyncio.sleep(min(0.1, remaining))


@router.post("/compute", response_model=ConformerJobResponse)
async def compute_conformer(request: ConformerRequest) -> ConformerJobResponse:
    """Submit a conformer computation job."""
    params = _build_conformer_params(request)
    cache_key = _compute_conformer_cache_key(
        molblock2d=request.molblock2d,
        geom_version=request.geom_version,
        params=params,
    )

    try:
        job_id, _ = job_queue.submit_job_with_cache_key(
            cache_key=cache_key,
            request_params={
                "job_type": "conformer",
                "molblock2d": request.molblock2d,
                "params": params,
                "geom_version": request.geom_version,
            },
        )
    except Exception as e:
        logger.error("Failed to submit conformer job", error=str(e))
        raise HTTPException(status_code=400, detail=str(e)) from e

    job = job_queue.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=500, detail="Job not found after creation")

    if not _job_is_terminal(job):
        await _wait_for_job(job_id, request.wait_ms)
        job = job_queue.get_job(job_id)
        if job is None:
            raise HTTPException(status_code=500, detail="Job not found after wait")

    return _build_conformer_response(job, cache_key)


@router.get("/jobs/{job_id}", response_model=ConformerJobResponse)
async def get_conformer_job(job_id: str) -> ConformerJobResponse:
    """Get status and result of a conformer computation job."""
    job = job_queue.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")

    return _build_conformer_response(job, job.cache_key)


@router.get("/artifacts/{cache_key}/{filename}")
async def get_conformer_artifact(cache_key: str, filename: str) -> Response:
    """Download a cached conformer artifact."""
    if ".." in filename or "/" in filename or "\\" in filename:
        raise HTTPException(status_code=400, detail="Invalid filename")

    artifact_path = job_queue.get_artifact_path(cache_key, filename)
    if artifact_path is None:
        raise HTTPException(status_code=404, detail="Artifact not found")

    content_type = "application/octet-stream"
    if filename.endswith(".sdf"):
        content_type = "chemical/x-mdl-sdfile"

    return FileResponse(
        path=artifact_path,
        media_type=content_type,
        filename=filename,
    )
