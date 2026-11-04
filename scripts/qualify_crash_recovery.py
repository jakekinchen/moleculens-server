"""Kill and recover an explicitly selected disposable qualification worker."""

import argparse
import json
import os
import subprocess
import time
import urllib.request
import uuid
from pathlib import Path

import psycopg2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--docker-context", required=True)
    parser.add_argument("--worker", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = args.url
    docker = ["docker", "--context", args.docker_context]

    def call(path, data=None):
        req = urllib.request.Request(
            base + path,
            data=json.dumps(data).encode() if data is not None else None,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=10) as response:
            return json.load(response)

    sdf = (
        (Path(__file__).parents[1] / "tests/contracts/water-v2000.mol")
        .read_text()
        .replace("water", "worker crash qualification " + uuid.uuid4().hex, 1)
    )
    r = call("/v1/orbitals/compute", {"sdfContent": sdf, "orbitals": ["homo", "lumo", "density"]})
    job = r["jobId"]
    start = time.monotonic()
    while r["status"] != "running":
        assert r["status"] != "error" and time.monotonic() - start < 20, r
        time.sleep(0.02)
        r = call("/v1/orbitals/jobs/" + job)
    with psycopg2.connect(os.environ["TEST_DATABASE_URL"]) as conn, conn.cursor() as cur:
        cur.execute("SELECT attempt,token FROM moleculens_v2.jobs WHERE id=%s", (job,))
        first = cur.fetchone()
    subprocess.run(docker + ["kill", args.worker], check=True, capture_output=True)
    state = json.loads(
        subprocess.check_output(
            docker + ["inspect", args.worker, "--format", "{{json .State}}"], text=True
        )
    )
    assert state["Running"] is False and state["Pid"] == 0, state
    subprocess.run(docker + ["start", args.worker], check=True, capture_output=True)
    while True:
        r = call("/v1/orbitals/jobs/" + job)
        if r["status"] in ("done", "error"):
            break
        assert time.monotonic() - start < 90, r
        time.sleep(0.2)
    assert r["status"] == "done", r
    with psycopg2.connect(os.environ["TEST_DATABASE_URL"]) as conn, conn.cursor() as cur:
        cur.execute("SELECT attempt,token,manifest FROM moleculens_v2.jobs WHERE id=%s", (job,))
        second = cur.fetchone()
    assert first[0] == 1 and second[0] == 2 and first[1] != second[1]
    assert second[2]["token"] == second[1]
    receipt = {
        "jobId": job,
        "originalAttempt": first[0],
        "finalAttempt": second[0],
        "originalToken": first[1],
        "publishedToken": second[1],
        "containerAfterSIGKILL": {
            "running": state["Running"],
            "pid": state["Pid"],
            "exitCode": state["ExitCode"],
        },
        "finalStatus": r["status"],
        "samePublicJobId": r["jobId"] == job,
        "wallSeconds": time.monotonic() - start,
        "engineBuild": second[2]["provenance"]["build"],
        "passed": True,
    }
    args.output.write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
