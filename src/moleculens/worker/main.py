"""Worker process main loop.

Claims jobs from Postgres queue and runs Psi4 computations.
"""

import signal
import time
import uuid

from moleculens.compute.psi4_runner import Psi4ComputationError, run_psi4_computation
from moleculens.core import get_logger, settings, setup_logging
from moleculens.db import JobQueue, init_db

logger = get_logger(__name__)

# Worker ID for this process
WORKER_ID = f"worker-{uuid.uuid4().hex[:8]}"

# Graceful shutdown flag
_shutdown_requested = False


def _signal_handler(signum: int, frame: object) -> None:
    """Handle shutdown signals."""
    global _shutdown_requested
    logger.info("Shutdown signal received", signal=signum)
    _shutdown_requested = True


def worker_loop() -> None:
    """Main worker loop - claim and process jobs."""
    global _shutdown_requested

    # Set up signal handlers
    signal.signal(signal.SIGTERM, _signal_handler)
    signal.signal(signal.SIGINT, _signal_handler)

    job_queue = JobQueue()
    poll_interval = settings.worker_poll_seconds

    logger.info(
        "Worker started",
        worker_id=WORKER_ID,
        poll_interval=poll_interval,
        psi4_threads=settings.psi4_num_threads,
        psi4_memory_mb=settings.psi4_memory_mb,
    )

    while not _shutdown_requested:
        try:
            # Try to claim a job
            job = job_queue.claim_next_job(WORKER_ID)

            if job is None:
                # No jobs available, wait and retry
                time.sleep(poll_interval)
                continue

            logger.info(
                "Processing job",
                job_id=job.id[:12],
                cache_key=job.cache_key[:12],
            )

            # Extract parameters
            params = job.request_params
            sdf_content = params["sdf_content"]
            method = params["method"]
            basis = params["basis"]
            grid_spacing = params["grid_spacing"]
            isovalue = params["isovalue"]
            orbitals = params["orbitals"]

            # Set up output directory
            output_dir = settings.molecule_cache_dir / job.cache_key
            output_dir.mkdir(parents=True, exist_ok=True)

            # Set up scratch directory for this job
            scratch_dir = settings.psi_scratch / job.id
            scratch_dir.mkdir(parents=True, exist_ok=True)

            try:
                # Run computation
                result = run_psi4_computation(
                    sdf_content=sdf_content,
                    method=method,
                    basis=basis,
                    grid_spacing=grid_spacing,
                    isovalue=isovalue,
                    orbitals=orbitals,
                    output_dir=output_dir,
                    scratch_dir=scratch_dir,
                )

                # Mark job complete
                job_queue.complete_job(
                    job_id=job.id,
                    result_meta=result.meta,
                    artifact_paths=result.artifact_files,
                    compute_time_ms=result.compute_time_ms,
                )

                logger.info(
                    "Job completed successfully",
                    job_id=job.id[:12],
                    compute_time_ms=result.compute_time_ms,
                )

            except Psi4ComputationError as e:
                logger.error("Computation failed", job_id=job.id[:12], error=str(e))
                job_queue.fail_job(job.id, str(e))

            except Exception as e:
                logger.exception("Unexpected error processing job", job_id=job.id[:12])
                job_queue.fail_job(job.id, f"Unexpected error: {e}")

            finally:
                # Clean up job-specific scratch directory
                try:
                    import shutil

                    if scratch_dir.exists():
                        shutil.rmtree(scratch_dir)
                except Exception as e:
                    logger.warning("Failed to clean scratch dir", error=str(e))

        except Exception:
            logger.exception("Error in worker loop")
            time.sleep(poll_interval * 5)  # Back off on errors

    logger.info("Worker shutting down", worker_id=WORKER_ID)


def run() -> None:
    """Entry point for moleculens-worker command."""
    setup_logging(settings.log_level)

    logger.info("Initializing worker", worker_id=WORKER_ID)

    # Ensure directories exist
    settings.ensure_directories()

    # Initialize database
    init_db()

    # Run worker loop
    worker_loop()


if __name__ == "__main__":
    run()
