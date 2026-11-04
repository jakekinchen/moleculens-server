"""Actual PostgreSQL qualification for the bounded legacy drain/rollback worker."""

import os
import subprocess
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit

import psycopg2
import pytest
from psycopg2 import sql
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from moleculens.compute.conformers import run_conformer_computation
from moleculens.db import queue as module
from moleculens.db.models import JobStatus


def test_legacy_admission_and_attempt_fence_on_postgres(tmp_path, monkeypatch):
    dsn = os.environ.get("TEST_DATABASE_URL")
    binary = Path(os.environ.get("GO_BINARY", "build/moleculens-linux"))
    if not dsn or not binary.is_file():
        pytest.skip(
            "TEST_DATABASE_URL and GO_BINARY are required for real PostgreSQL qualification"
        )
    admin = psycopg2.connect(dsn)
    admin.autocommit = True
    name = "mllegacy_" + uuid.uuid4().hex
    with admin.cursor() as cursor:
        cursor.execute(
            sql.SQL("CREATE DATABASE {} TEMPLATE template0 ENCODING 'UTF8'").format(
                sql.Identifier(name)
            )
        )
    uri = urlsplit(dsn)
    isolated = urlunsplit(uri._replace(path="/" + name))
    engine = None
    try:
        subprocess.run(
            [str(binary.resolve()), "migrate"],
            env={**os.environ, "DATABASE_URL": isolated},
            check=True,
        )
        engine = create_engine(isolated, pool_size=8, max_overflow=24)
        monkeypatch.setattr(module, "_engine", engine)
        monkeypatch.setattr(module, "_session_factory", sessionmaker(bind=engine))
        queue = module.JobQueue(tmp_path)
        sdf = (Path(__file__).parents[1] / "contracts/water-v2000.mol").read_text()
        params = {"job_type": "conformer", "molblock2d": sdf, "params": {"seed": 42}}
        with ThreadPoolExecutor(max_workers=16) as threads:
            results = list(
                threads.map(
                    lambda _: queue.submit_job_with_cache_key("same-key", params), range(32)
                )
            )
        assert len({r[0] for r in results}) == 1
        with ThreadPoolExecutor(max_workers=8) as threads:
            claims = list(threads.map(lambda i: queue.claim_next_job(f"worker-{i}"), range(8)))
        claimed = [job for job in claims if job is not None]
        assert len(claimed) == 1
        job = claimed[0]
        assert queue.complete_job(job.id, {}, [], 1, attempt_token="stale") is False
        assert queue.fail_job(job.id, "late", attempt_token="stale") is False
        output = tmp_path / job.cache_key / job.attempt_token
        result = run_conformer_computation(sdf, seed=42, output_dir=output)
        assert (
            queue.complete_job(
                job.id,
                result.meta,
                result.artifact_files,
                result.compute_time_ms,
                attempt_token=job.attempt_token,
            )
            is True
        )
        assert queue.fail_job(job.id, "late", attempt_token=job.attempt_token) is False
        assert queue.get_job(job.id).status == JobStatus.DONE
        assert queue.get_artifact_path(job.cache_key, "structure.sdf") == output / "structure.sdf"
    finally:
        if engine is not None:
            engine.dispose()
        with admin.cursor() as cursor:
            cursor.execute(sql.SQL("DROP DATABASE {} WITH(FORCE)").format(sql.Identifier(name)))
        admin.close()
