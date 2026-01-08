"""Worker process main loop.

Claims jobs from Postgres queue and runs Psi4 computations.
"""

import signal
import time
import uuid

from moleculens.compute.conformers import (
    ConformerComputationError,
    run_conformer_computation,
)
from moleculens.compute.electrostatics import (
    ElectrostaticsComputationError,
    run_electrostatics_computation,
)
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

            # Detect job type - check explicit type first, then fallback markers
            job_type = params.get("job_type")
            if job_type is None:
                # Check for electrostatics marker set by electrostatics_routes
                if params.get("basis") == "electrostatics" or params.get("orbitals") == [
                    "electrostatics"
                ]:
                    job_type = "electrostatics"
                else:
                    job_type = "orbital"

            # Set up output directory
            output_dir = settings.molecule_cache_dir / job.cache_key
            output_dir.mkdir(parents=True, exist_ok=True)

            # Set up scratch directory for this job
            scratch_dir = settings.psi_scratch / job.id
            scratch_dir.mkdir(parents=True, exist_ok=True)

            try:
                if job_type == "electrostatics":
                    # Run electrostatics computation
                    surface_params = params.get("surface", {})
                    potential_params = params.get("potential", {})

                    result = run_electrostatics_computation(
                        sdf_content=params["sdf_content"],
                        charge=params.get("charge", 0),
                        multiplicity=params.get("multiplicity", 1),
                        method=params.get("method", "xtb-gfn2"),
                        grid_spacing=surface_params.get("grid_spacing", 0.25),
                        probe_radius=surface_params.get("probe_radius", 1.4),
                        padding=surface_params.get("padding", 3.0),
                        softening_epsilon=potential_params.get("softening_epsilon", 0.3),
                        clamp_percentiles=tuple(
                            potential_params.get("clamp_percentiles", [5.0, 95.0])
                        ),
                        output_dir=output_dir,
                    )
                elif job_type == "conformer":
                    conformer_params = params.get("params", {})
                    molblock2d = params.get("molblock2d")
                    if not molblock2d:
                        raise ConformerComputationError("Missing molblock2d for conformer job")

                    result = run_conformer_computation(
                        molblock2d=molblock2d,
                        method=conformer_params.get("method", "etkdg_v3"),
                        opt=conformer_params.get("opt", "none"),
                        max_attempts=conformer_params.get("max_attempts", 10),
                        max_opt_iters=conformer_params.get("max_opt_iters", 200),
                        add_hs=conformer_params.get("add_hs", False),
                        geom_version=params.get("geom_version", "unknown"),
                        output_dir=output_dir,
                    )
                else:
                    # Run Psi4 orbital computation
                    result = run_psi4_computation(
                        sdf_content=params["sdf_content"],
                        method=params["method"],
                        basis=params["basis"],
                        grid_spacing=params["grid_spacing"],
                        isovalue=params["isovalue"],
                        orbitals=params["orbitals"],
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
                    job_type=job_type,
                    compute_time_ms=result.compute_time_ms,
                )

            except (
                Psi4ComputationError,
                ElectrostaticsComputationError,
                ConformerComputationError,
            ) as e:
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
