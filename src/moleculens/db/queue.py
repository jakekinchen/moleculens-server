"""Postgres-based job queue implementation.

Uses SELECT ... FOR UPDATE SKIP LOCKED for efficient job claiming.
"""

import hashlib
import json
import uuid
from collections.abc import Generator
from contextlib import contextmanager
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

from sqlalchemy import create_engine, select, text
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session, sessionmaker

from moleculens.core import get_logger, settings
from moleculens.db.models import Base, CacheEntry, Job, JobStatus

logger = get_logger(__name__)

# Module-level engine and session factory
_engine: Engine | None = None
_session_factory: sessionmaker[Session] | None = None


def get_db_engine() -> Engine:
    """Get or create the database engine."""
    global _engine
    if _engine is None:
        _engine = create_engine(
            settings.database_url,
            pool_size=5,
            max_overflow=10,
            pool_pre_ping=True,
        )
    return _engine


def get_db_session() -> Session:
    """Get a new database session."""
    global _session_factory
    if _session_factory is None:
        _session_factory = sessionmaker(bind=get_db_engine())
    return _session_factory()


@contextmanager
def db_session() -> Generator[Session, None, None]:
    """Context manager for database sessions."""
    session = get_db_session()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def init_db() -> None:
    """Initialize database schema."""
    engine = get_db_engine()
    Base.metadata.create_all(engine)
    logger.info("Database initialized")


def compute_cache_key(
    geometry_hash: str,
    method: str,
    basis: str,
    grid_spacing: float,
    isovalue: float,
    orbitals: list[str],
    inchi_key: str | None = None,
) -> str:
    """Compute stable cache key from computation parameters.

    Args:
        geometry_hash: Hash of molecular geometry
        method: Computational method (e.g., 'scf', 'b3lyp')
        basis: Basis set (e.g., 'sto-3g', '6-31g*')
        grid_spacing: Grid spacing in Angstrom
        isovalue: Isosurface value
        orbitals: List of orbital types to compute
        inchi_key: Optional InChIKey for deduplication

    Returns:
        64-character hex cache key
    """
    # Use InChIKey if provided, otherwise geometry hash
    mol_id = inchi_key if inchi_key else geometry_hash

    # Build deterministic string
    params = {
        "mol_id": mol_id,
        "method": method.lower(),
        "basis": basis.lower().replace("*", "star"),
        "grid_spacing": f"{grid_spacing:.4f}",
        "isovalue": f"{isovalue:.4f}",
        "orbitals": sorted(orbitals),
    }

    content = json.dumps(params, sort_keys=True)
    return hashlib.sha256(content.encode()).hexdigest()


class JobQueue:
    """Postgres-based job queue with caching support."""

    def __init__(self, cache_dir: Path | None = None):
        """Initialize job queue.

        Args:
            cache_dir: Directory for artifact caching
        """
        self.cache_dir = cache_dir or settings.molecule_cache_dir

    def _active_job_is_stale(self, job: Job) -> bool:
        """Return whether an active job has exceeded the stale threshold."""
        if job.status not in {JobStatus.PENDING, JobStatus.QUEUED, JobStatus.RUNNING}:
            return False

        reference_time = job.started_at or job.created_at
        if reference_time is None:
            return False
        if reference_time.tzinfo is None:
            reference_time = reference_time.replace(tzinfo=UTC)

        return reference_time <= datetime.now(UTC) - timedelta(seconds=settings.job_stale_seconds)

    def _mark_active_job_stale(self, job: Job) -> None:
        """Mark an active job as stale so a fresh replacement can be created."""
        original_status = job.status.value
        job.status = JobStatus.ERROR
        job.completed_at = datetime.now(UTC)
        job.error_message = (
            f"Stale {original_status} job auto-cleared after "
            f"{settings.job_stale_seconds}s without completion"
        )

        logger.warning(
            "Auto-cleared stale active job",
            job_id=job.id[:12],
            cache_key=job.cache_key[:12],
            previous_status=original_status,
        )

    def _get_reusable_active_job_id(self, session: Session, cache_key: str) -> str | None:
        """Return a fresh active job for a cache key, clearing stale ones first."""
        active_jobs = (
            session.execute(
                select(Job)
                .where(Job.cache_key == cache_key)
                .where(Job.status.in_([JobStatus.PENDING, JobStatus.QUEUED, JobStatus.RUNNING]))
                .order_by(Job.created_at.desc())
            )
            .scalars()
            .all()
        )

        for active_job in active_jobs:
            if self._active_job_is_stale(active_job):
                self._mark_active_job_stale(active_job)
                continue

            logger.info(
                "Found active job for cache key",
                cache_key=cache_key[:12],
                job_id=active_job.id[:12],
                status=active_job.status.value,
            )
            return active_job.id

        return None

    def _required_artifact_paths(self, artifact_paths: list[str]) -> list[str]:
        """Return the minimum artifact set required for a cache entry to be reusable."""
        result_manifests = [path for path in artifact_paths if path.endswith("result.json")]
        return result_manifests or artifact_paths

    def _cache_artifacts_exist(self, cache_key: str, artifact_paths: list[str]) -> bool:
        """Check whether the required cache artifacts still exist on disk."""
        required_paths = self._required_artifact_paths(artifact_paths)
        if not required_paths:
            return False

        cache_root = self.cache_dir / cache_key
        return all((cache_root / path).exists() for path in required_paths)

    def _invalidate_cache_key(self, session: Session, cache_key: str, error_message: str) -> None:
        """Invalidate a cache entry and mark stale completed jobs as failed."""
        cache_entry = session.get(CacheEntry, cache_key)
        if cache_entry is not None:
            session.delete(cache_entry)

        stale_jobs = (
            session.execute(
                select(Job).where(Job.cache_key == cache_key).where(Job.status == JobStatus.DONE)
            )
            .scalars()
            .all()
        )

        invalidated_jobs = 0
        for stale_job in stale_jobs:
            stale_job.status = JobStatus.ERROR
            stale_job.error_message = error_message
            stale_job.completed_at = datetime.now(UTC)
            invalidated_jobs += 1

        logger.warning(
            "Invalidated stale cache entry",
            cache_key=cache_key[:12],
            invalidated_jobs=invalidated_jobs,
            reason=error_message[:200],
        )

    def invalidate_cache_key(self, cache_key: str, error_message: str) -> None:
        """Public wrapper for invalidating a corrupt cache entry."""
        with db_session() as session:
            self._invalidate_cache_key(session, cache_key, error_message)

    def _reuse_cached_job_if_valid(
        self, session: Session, cache_key: str
    ) -> tuple[str, bool] | None:
        """Return a reusable cached job if its artifacts still exist."""
        cache_entry = session.get(CacheEntry, cache_key)
        if cache_entry is None:
            return None

        cache_error = f"Cached result missing for cache key {cache_key[:12]}"
        if not self._cache_artifacts_exist(cache_key, cache_entry.artifact_paths):
            self._invalidate_cache_key(session, cache_key, cache_error)
            return None

        cache_entry.hit_count += 1
        cache_entry.last_accessed_at = datetime.now(UTC)

        existing_job = session.execute(
            select(Job)
            .where(Job.cache_key == cache_key)
            .where(Job.status == JobStatus.DONE)
            .limit(1)
        ).scalar_one_or_none()

        if existing_job is None:
            return None

        if existing_job.result_meta is None:
            existing_job.result_meta = cache_entry.result_meta

        logger.info(
            "Returning cached result",
            cache_key=cache_key[:12],
            job_id=existing_job.id[:12],
        )
        return existing_job.id, True

    def submit_job(
        self,
        sdf_content: str,
        method: str,
        basis: str,
        grid_spacing: float,
        isovalue: float,
        orbitals: list[str],
        inchi_key: str | None = None,
        geometry_hash: str | None = None,
    ) -> tuple[str, bool]:
        """Submit a new computation job or return cached result.

        Args:
            sdf_content: SDF file content
            method: Computational method
            basis: Basis set
            grid_spacing: Grid spacing in Angstrom
            isovalue: Isosurface value
            orbitals: List of orbital types
            inchi_key: Optional InChIKey
            geometry_hash: Pre-computed geometry hash

        Returns:
            Tuple of (job_id, is_cached)
        """
        # Compute cache key
        if geometry_hash is None:
            # Parse SDF to get geometry hash
            from moleculens.compute.sdf_parser import parse_sdf

            mol = parse_sdf(sdf_content)
            geometry_hash = mol.geometry_hash()

        cache_key = compute_cache_key(
            geometry_hash=geometry_hash,
            method=method,
            basis=basis,
            grid_spacing=grid_spacing,
            isovalue=isovalue,
            orbitals=orbitals,
            inchi_key=inchi_key,
        )

        with db_session() as session:
            cached_job = self._reuse_cached_job_if_valid(session, cache_key)
            if cached_job is not None:
                return cached_job

            active_job_id = self._get_reusable_active_job_id(session, cache_key)
            if active_job_id:
                return active_job_id, False

            # Create new job
            job_id = str(uuid.uuid4())
            job = Job(
                id=job_id,
                cache_key=cache_key,
                status=JobStatus.PENDING,
                request_params={
                    "sdf_content": sdf_content,
                    "method": method,
                    "basis": basis,
                    "grid_spacing": grid_spacing,
                    "isovalue": isovalue,
                    "orbitals": orbitals,
                    "inchi_key": inchi_key,
                },
            )
            session.add(job)

            logger.info("Created new job", job_id=job_id[:12], cache_key=cache_key[:12])
            return job_id, False

    def submit_job_with_cache_key(
        self,
        cache_key: str,
        request_params: dict[str, Any],
    ) -> tuple[str, bool]:
        """Submit a new computation job with an explicit cache key.

        Args:
            cache_key: Pre-computed cache key
            request_params: Job request parameters to store

        Returns:
            Tuple of (job_id, is_cached)
        """
        with db_session() as session:
            cached_job = self._reuse_cached_job_if_valid(session, cache_key)
            if cached_job is not None:
                return cached_job

            active_job_id = self._get_reusable_active_job_id(session, cache_key)
            if active_job_id:
                return active_job_id, False

            # Create new job
            job_id = str(uuid.uuid4())
            job = Job(
                id=job_id,
                cache_key=cache_key,
                status=JobStatus.PENDING,
                request_params=request_params,
            )
            session.add(job)

            logger.info("Created new job", job_id=job_id[:12], cache_key=cache_key[:12])
            return job_id, False

    def get_job(self, job_id: str) -> Job | None:
        """Get job by ID."""
        with db_session() as session:
            job = session.get(Job, job_id)
            if job:
                # Detach from session for returning
                session.expunge(job)
            return job

    def get_job_by_cache_key(self, cache_key: str) -> Job | None:
        """Get most recent job for a cache key."""
        with db_session() as session:
            job = session.execute(
                select(Job)
                .where(Job.cache_key == cache_key)
                .order_by(Job.created_at.desc())
                .limit(1)
            ).scalar_one_or_none()
            if job:
                session.expunge(job)
            return job

    def claim_next_job(self, worker_id: str) -> Job | None:
        """Claim the next pending job for processing.

        Uses SELECT ... FOR UPDATE SKIP LOCKED for safe concurrent access.

        Args:
            worker_id: Unique identifier for the claiming worker

        Returns:
            Claimed job or None if no jobs available
        """
        with db_session() as session:
            while True:
                # Use raw SQL for FOR UPDATE SKIP LOCKED
                result = session.execute(
                    text("""
                        SELECT id FROM jobs
                        WHERE status = 'pending'
                        ORDER BY created_at ASC
                        LIMIT 1
                        FOR UPDATE SKIP LOCKED
                    """)
                ).fetchone()

                if result is None:
                    return None

                job_id = result[0]
                job = session.get(Job, job_id)
                if job is None:
                    continue

                if self._active_job_is_stale(job):
                    self._mark_active_job_stale(job)
                    session.flush()
                    continue

                job.status = JobStatus.RUNNING
                job.started_at = datetime.now(UTC)
                job.worker_id = worker_id
                job.attempt_count += 1
                session.flush()

                logger.info(
                    "Claimed job",
                    job_id=job_id[:12],
                    worker_id=worker_id[:12],
                    attempt=job.attempt_count,
                )

                session.expunge(job)
                return job

    def complete_job(
        self,
        job_id: str,
        result_meta: dict[str, Any],
        artifact_paths: list[str],
        compute_time_ms: float,
    ) -> None:
        """Mark job as completed and create cache entry.

        Args:
            job_id: Job ID
            result_meta: Result metadata
            artifact_paths: List of artifact file paths
            compute_time_ms: Computation time in milliseconds
        """
        with db_session() as session:
            job = session.get(Job, job_id)
            if job is None:
                logger.error("Job not found for completion", job_id=job_id[:12])
                return

            job.status = JobStatus.DONE
            job.completed_at = datetime.now(UTC)
            job.result_meta = result_meta
            job.compute_time_ms = compute_time_ms

            # Create or update cache entry
            cache_entry = session.get(CacheEntry, job.cache_key)
            if cache_entry is None:
                cache_entry = CacheEntry(
                    cache_key=job.cache_key,
                    params_hash=job.cache_key,  # Already a hash
                    result_meta=result_meta,
                    artifact_paths=artifact_paths,
                )
                session.add(cache_entry)
            else:
                cache_entry.result_meta = result_meta
                cache_entry.artifact_paths = artifact_paths

            logger.info(
                "Job completed",
                job_id=job_id[:12],
                cache_key=job.cache_key[:12],
                compute_time_ms=compute_time_ms,
            )

    def fail_job(self, job_id: str, error_message: str) -> None:
        """Mark job as failed.

        Args:
            job_id: Job ID
            error_message: Error description
        """
        with db_session() as session:
            job = session.get(Job, job_id)
            if job is None:
                logger.error("Job not found for failure", job_id=job_id[:12])
                return

            job.status = JobStatus.ERROR
            job.completed_at = datetime.now(UTC)
            job.error_message = error_message

            logger.error(
                "Job failed",
                job_id=job_id[:12],
                error=error_message[:200],
            )

    def _update_job_params(self, job_id: str, params: dict[str, Any]) -> None:
        """Update job request parameters.

        Used to add additional parameters after job creation.

        Args:
            job_id: Job ID
            params: Parameters to merge into request_params
        """
        with db_session() as session:
            job = session.get(Job, job_id)
            if job is None:
                return

            # Merge new params into existing
            current_params = dict(job.request_params or {})
            current_params.update(params)
            job.request_params = current_params

    def get_cache_entry(self, cache_key: str) -> CacheEntry | None:
        """Get cache entry by key."""
        with db_session() as session:
            entry = session.get(CacheEntry, cache_key)
            if entry:
                session.expunge(entry)
            return entry

    def get_artifact_path(self, cache_key: str, filename: str) -> Path | None:
        """Get full path to a cached artifact.

        Args:
            cache_key: Cache key
            filename: Artifact filename

        Returns:
            Full path to artifact or None if not found
        """
        cache_entry = self.get_cache_entry(cache_key)
        if cache_entry is None:
            return None

        # Check if filename is in the artifact list
        if filename not in cache_entry.artifact_paths:
            return None

        full_path = self.cache_dir / cache_key / filename
        if full_path.exists():
            return full_path

        return None
