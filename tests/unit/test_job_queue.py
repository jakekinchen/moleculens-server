"""Unit tests for stale cache invalidation and job recovery."""

from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import SimpleNamespace

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from moleculens.api import electrostatics_routes, routes
from moleculens.db import queue as queue_module
from moleculens.db.models import Base, CacheEntry, Job, JobStatus
from moleculens.db.queue import JobQueue


@pytest.fixture
def isolated_job_queue(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> JobQueue:
    """Provide an isolated SQLite-backed job queue and cache directory."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    engine = create_engine(f"sqlite:///{tmp_path / 'queue.db'}")
    monkeypatch.setattr(queue_module, "_engine", engine)
    monkeypatch.setattr(queue_module, "_session_factory", sessionmaker(bind=engine))
    monkeypatch.setattr(queue_module.settings, "molecule_cache_dir", cache_dir)

    Base.metadata.create_all(engine)

    routes.job_queue = JobQueue(cache_dir=cache_dir)
    return JobQueue(cache_dir=cache_dir)


def _seed_done_job(
    job_queue: JobQueue,
    *,
    geometry_hash: str,
    result_meta: dict | None,
    artifact_paths: list[str],
) -> tuple[str, str]:
    """Create a completed job + cache entry for the given cache key."""
    job_id, is_cached = job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        geometry_hash=geometry_hash,
    )
    assert is_cached is False

    job = job_queue.get_job(job_id)
    assert job is not None
    cache_key = job.cache_key

    with queue_module.db_session() as session:
        seeded_job = session.get(Job, job_id)
        assert seeded_job is not None
        seeded_job.status = JobStatus.DONE
        seeded_job.result_meta = result_meta
        seeded_job.compute_time_ms = 123.0

        session.add(
            CacheEntry(
                cache_key=cache_key,
                params_hash=cache_key,
                result_meta={"method": "scf", "basis": "sto-3g"},
                artifact_paths=artifact_paths,
            )
        )

    return job_id, cache_key


def _age_active_job(
    job_id: str,
    *,
    status: JobStatus,
    seconds_ago: int,
) -> None:
    """Age an active job so stale-job recovery paths can be exercised."""
    with queue_module.db_session() as session:
        job = session.get(Job, job_id)
        assert job is not None
        timestamp = datetime.now(UTC) - timedelta(seconds=seconds_ago)
        job.status = status
        job.created_at = timestamp
        job.started_at = timestamp if status == JobStatus.RUNNING else None
        job.completed_at = None
        job.error_message = None


def test_submit_job_reuses_valid_cache_and_restores_missing_result_meta(
    isolated_job_queue: JobQueue,
) -> None:
    job_id, cache_key = _seed_done_job(
        isolated_job_queue,
        geometry_hash="geom-valid",
        result_meta=None,
        artifact_paths=["result.json"],
    )

    cache_root = isolated_job_queue.cache_dir / cache_key
    cache_root.mkdir()
    (cache_root / "result.json").write_text('{"orbitals": {}, "meta": {"method": "scf"}}')

    reused_job_id, is_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        geometry_hash="geom-valid",
    )

    assert is_cached is True
    assert reused_job_id == job_id

    reused_job = isolated_job_queue.get_job(job_id)
    assert reused_job is not None
    assert reused_job.result_meta == {"method": "scf", "basis": "sto-3g"}


def test_submit_job_invalidates_stale_cache_and_requeues(
    isolated_job_queue: JobQueue,
) -> None:
    stale_job_id, cache_key = _seed_done_job(
        isolated_job_queue,
        geometry_hash="geom-stale",
        result_meta={"method": "scf", "basis": "sto-3g"},
        artifact_paths=["result.json"],
    )

    new_job_id, is_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        geometry_hash="geom-stale",
    )

    assert is_cached is False
    assert new_job_id != stale_job_id

    stale_job = isolated_job_queue.get_job(stale_job_id)
    assert stale_job is not None
    assert stale_job.status == JobStatus.ERROR
    assert stale_job.error_message == f"Cached result missing for cache key {cache_key[:12]}"

    new_job = isolated_job_queue.get_job(new_job_id)
    assert new_job is not None
    assert new_job.status == JobStatus.PENDING

    assert isolated_job_queue.get_cache_entry(cache_key) is None


@pytest.mark.asyncio
async def test_get_job_status_reports_error_for_done_job_with_missing_result(
    isolated_job_queue: JobQueue,
) -> None:
    stale_job_id, cache_key = _seed_done_job(
        isolated_job_queue,
        geometry_hash="geom-poll-stale",
        result_meta={"method": "scf", "basis": "sto-3g"},
        artifact_paths=["result.json"],
    )

    response = await routes.get_job_status(stale_job_id)

    assert response.status == "error"
    assert response.error_message == f"Cached result missing for cache key {cache_key[:12]}"
    assert response.result is None

    stale_job = isolated_job_queue.get_job(stale_job_id)
    assert stale_job is not None
    assert stale_job.status == JobStatus.ERROR
    assert isolated_job_queue.get_cache_entry(cache_key) is None


def test_submit_job_does_not_reuse_by_inchi_key_alone(isolated_job_queue: JobQueue) -> None:
    first_job_id, first_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        inchi_key="SAME-INCHI",
        geometry_hash="geom-a",
    )
    second_job_id, second_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        inchi_key="SAME-INCHI",
        geometry_hash="geom-b",
    )

    assert first_cached is False
    assert second_cached is False
    assert first_job_id != second_job_id


def test_submit_job_reuses_explicit_cache_identity(isolated_job_queue: JobQueue) -> None:
    first_job_id, first_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        cache_identity="client-cache-key",
        geometry_hash="geom-a",
    )
    second_job_id, second_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        cache_identity="client-cache-key",
        geometry_hash="geom-b",
    )

    assert first_cached is False
    assert second_cached is False
    assert first_job_id == second_job_id


@pytest.mark.asyncio
async def test_compute_orbitals_forwards_client_cache_key(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    class StubQueue:
        def submit_job(self, **kwargs: object) -> tuple[str, bool]:
            captured.update(kwargs)
            return "job-1", False

        def get_job(self, job_id: str) -> SimpleNamespace:
            return SimpleNamespace(
                id=job_id,
                status=JobStatus.PENDING,
                cache_key="server-cache-key",
                created_at=datetime.now(UTC),
                compute_time_ms=None,
            )

    monkeypatch.setattr(routes, "job_queue", StubQueue())

    request = routes.ComputeRequest(
        sdf_content="water\n  test\n\n  1  0  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$",
        clientCacheKey="client-cache-key",
    )

    response = await routes.compute_orbitals(request)

    assert response.cache_key == "server-cache-key"
    assert captured["cache_identity"] == "client-cache-key"


@pytest.mark.asyncio
async def test_compute_electrostatics_returns_persisted_cache_key(
    monkeypatch: pytest.MonkeyPatch,
    water_sdf: str,
) -> None:
    captured: dict[str, object] = {}

    class StubQueue:
        def submit_job_with_cache_key(
            self, *, cache_key: str, request_params: dict[str, object]
        ) -> tuple[str, bool]:
            captured["cache_key"] = cache_key
            captured["request_params"] = request_params
            return "job-1", False

        def get_job(self, job_id: str) -> SimpleNamespace:
            return SimpleNamespace(
                id=job_id,
                status=JobStatus.PENDING,
                cache_key="persisted-cache-key",
                created_at=datetime.now(UTC),
                compute_time_ms=None,
            )

    monkeypatch.setattr(electrostatics_routes, "job_queue", StubQueue())

    request = electrostatics_routes.ElectrostaticsRequest(sdfContent=water_sdf)
    response = await electrostatics_routes.compute_electrostatics(request)

    assert response.cache_key == "persisted-cache-key"
    assert captured["cache_key"] == captured["request_params"]["cache_identity"]


def test_submit_job_replaces_stale_pending_job(
    isolated_job_queue: JobQueue,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(queue_module.settings, "job_stale_seconds", 60)

    stale_job_id, is_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        geometry_hash="geom-pending-stale",
    )
    assert is_cached is False

    _age_active_job(stale_job_id, status=JobStatus.PENDING, seconds_ago=120)

    new_job_id, is_cached = isolated_job_queue.submit_job(
        sdf_content="not-used",
        method="scf",
        basis="sto-3g",
        grid_spacing=0.25,
        isovalue=0.05,
        orbitals=["homo"],
        geometry_hash="geom-pending-stale",
    )

    assert is_cached is False
    assert new_job_id != stale_job_id

    stale_job = isolated_job_queue.get_job(stale_job_id)
    assert stale_job is not None
    assert stale_job.status == JobStatus.ERROR
    assert stale_job.error_message == "Stale pending job auto-cleared after 60s without completion"


def test_submit_job_with_cache_key_replaces_stale_running_job(
    isolated_job_queue: JobQueue,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(queue_module.settings, "job_stale_seconds", 60)

    stale_job_id, is_cached = isolated_job_queue.submit_job_with_cache_key(
        cache_key="cache-key-running-stale",
        request_params={"kind": "electrostatics"},
    )
    assert is_cached is False

    _age_active_job(stale_job_id, status=JobStatus.RUNNING, seconds_ago=120)

    new_job_id, is_cached = isolated_job_queue.submit_job_with_cache_key(
        cache_key="cache-key-running-stale",
        request_params={"kind": "electrostatics"},
    )

    assert is_cached is False
    assert new_job_id != stale_job_id

    stale_job = isolated_job_queue.get_job(stale_job_id)
    assert stale_job is not None
    assert stale_job.status == JobStatus.ERROR
    assert stale_job.error_message == "Stale running job auto-cleared after 60s without completion"
