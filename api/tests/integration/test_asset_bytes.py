import os
import subprocess
import time
from pathlib import Path

import pytest
import requests


@pytest.mark.integration
@pytest.mark.slow
def test_render_model_bytes():
    """Build the Docker image, run the server, render a model and ensure
    the returned payload is non-empty (>1 kB).
    This is a smoke-test proving the full containerized stack works.
    """

    image = "moleculens-test:latest"

    # Build image (use cache when available)
    subprocess.run(
        ["docker", "build", "-t", image, Path(__file__).resolve().parents[3]],
        check=True,
    )

    # Bind container port 8000 to a random free host port ("0").
    container_id = (
        subprocess.check_output(["docker", "run", "-d", "-p", "0:8000", image])
        .decode()
        .strip()
    )

    # Discover which host port Docker chose.
    host_port = (
        subprocess.check_output(["docker", "port", container_id, "8000/tcp"])
        .decode()
        .split(":")[-1]
        .strip()
    )

    def _cleanup() -> None:
        subprocess.run(
            ["docker", "rm", "-f", container_id],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    try:
        # Wait until /health reports healthy (max 120 s)
        deadline = time.time() + 120
        while time.time() < deadline:
            try:
                res = requests.get(f"http://localhost:{host_port}/health", timeout=5)
                if res.status_code == 200:
                    break
            except requests.RequestException:
                pass
            time.sleep(2)
        else:
            pytest.fail("API inside container did not become healthy within 120 s")

        payload = {
            "description": "fetch 1ubq; hide everything; show cartoon; orient",
            "format": "model",
        }
        response = requests.post(
            f"http://localhost:{host_port}/render", json=payload, timeout=180
        )
        response.raise_for_status()

        size = len(response.content)
        assert size > 1024, f"Expected model >1 kB, got {size} bytes"

    finally:
        _cleanup()
