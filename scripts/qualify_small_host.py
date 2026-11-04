"""Exercise the disposable complete stack; requires a VM resource sampler alongside it."""

import argparse
import concurrent.futures
import json
import re
import subprocess
import threading
import time
import urllib.error
import urllib.request
import uuid
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--docker-context", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    docker = ["docker", "--context", args.docker_context]
    receipts = {
        "cases": [],
        "scope": "local 2-CPU/4-GiB VM; Linux amd64 images under Rosetta; not an OVH CPU benchmark",
    }
    busy_reads = []

    def save():
        args.output.write_text(json.dumps(receipts, indent=2) + "\n")

    def call(path, body=None, retry_busy=False):
        req = urllib.request.Request(
            args.url + path,
            data=None if body is None else json.dumps(body).encode(),
            headers={"Content-Type": "application/json"},
        )
        deadline = time.monotonic() + 2
        while True:
            try:
                with urllib.request.urlopen(req, timeout=25) as response:
                    return json.load(response)
            except urllib.error.HTTPError as error:
                detail = json.load(error)
                if (
                    retry_busy
                    and error.code == 503
                    and detail.get("detail", {}).get("code") == "busy"
                    and time.monotonic() < deadline
                ):
                    busy_reads.append(time.monotonic())
                    time.sleep(0.1)
                    continue
                raise RuntimeError({"status": error.code, "detail": detail}) from error

    def fixture(name):
        text = (
            Path(__file__).parents[1] / "tests/contracts/small-host" / (name + ".mol")
        ).read_text()
        return text.replace(name, name + " " + uuid.uuid4().hex, 1)

    def submit(path, body):
        response = call(path + "/compute", body)
        return response["jobId"]

    def poll(path, job, seconds=330):
        start = time.monotonic()
        while time.monotonic() - start < seconds:
            response = call(path + "/jobs/" + job, retry_busy=True)
            if response["status"] in ("done", "error", "failed"):
                return response
            time.sleep(0.25)
        raise RuntimeError("job deadline " + job)

    def row(job):
        assert re.fullmatch(r"g_[0-9a-f]+", job)
        sql = (
            "SELECT json_build_object('attempt',attempt,'token',token,'owner',owner,'state',state)::text FROM moleculens_v2.jobs WHERE id='"
            + job
            + "'"
        )
        return json.loads(
            subprocess.check_output(
                docker
                + [
                    "exec",
                    "moleculens-vps1-postgres-1",
                    "psql",
                    "-U",
                    "moleculens",
                    "-d",
                    "moleculens",
                    "-Atc",
                    sql,
                ],
                text=True,
            )
        )

    def record(name, response, start):
        entry = {
            "case": name,
            "jobId": response.get("jobId"),
            "status": response["status"],
            "wallSeconds": time.monotonic() - start,
            "computeTimeMs": response.get("computeTimeMs"),
        }
        if response["status"] == "done":
            if "result" in response:
                result = response["result"]
                entry["meta"] = result.get("meta")
                if result.get("density"):
                    entry["densityVertices"] = result["density"]["vertexCount"]
                if result.get("surface"):
                    entry["surfaceVertices"] = result["surface"]["mesh"]["vertexCount"]
            else:
                entry["meta"] = response.get("meta")
                assert "V2000" in response["sdf3d"]
        else:
            entry["error"] = response
        receipts["cases"].append(entry)
        save()
        assert response["status"] == "done", entry
        return entry

    print(
        "Starting caffeine B3LYP/6-31G(d) with concurrent naming, conformer and two pollers",
        flush=True,
    )
    start = time.monotonic()
    job = submit(
        "/v1/orbitals",
        {
            "sdfContent": fixture("caffeine"),
            "method": "b3lyp",
            "basis": "6-31g(d)",
            "orbitals": ["homo", "lumo", "density"],
            "gridSpacing": 0.25,
        },
    )
    while call("/v1/orbitals/jobs/" + job)["status"] == "pending":
        assert time.monotonic() - start < 15
        time.sleep(0.05)
    done = threading.Event()
    latencies = []

    def observer():
        while not done.is_set():
            now = time.monotonic()
            call("/v1/orbitals/jobs/" + job, retry_busy=True)
            latencies.append(time.monotonic() - now)
            time.sleep(0.2)

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        observers = [executor.submit(observer) for _ in range(2)]
        try:
            name_start = time.monotonic()
            naming = call(
                "/v1/nomenclature/explain",
                {"smiles": "O=C(O)c1ccccc1", "molblock": fixture("benzoic"), "verifyOpsin": True},
            )
            naming_seconds = time.monotonic() - name_start
            assert naming["verified"] and naming["opsinStatus"] == "matched", naming
            conformer_start = time.monotonic()
            short_job = submit(
                "/v1/conformers",
                {
                    "molblock2d": fixture("caffeine"),
                    "geomVersion": "etkdgv3",
                    "seed": 42,
                    "forceField": "mmff",
                    "waitMs": 0,
                },
            )
            conformer = poll("/v1/conformers", short_job, 45)
            record("caffeine-conformer-during-psi4", conformer, conformer_start)
            long_state = call("/v1/orbitals/jobs/" + job)["status"]
            receipts["responsiveness"] = {
                "namingSeconds": naming_seconds,
                "naming": naming["displayName"],
                "opsinStatus": naming["opsinStatus"],
                "longJobStateAfterShortWork": long_state,
                "distinctWorkerOwners": row(job)["owner"] != row(short_job)["owner"],
            }
            save()
            assert long_state == "running", "long calculation did not overlap the short operations"
            assert naming_seconds < 9
            response = poll("/v1/orbitals", job)
            record("caffeine-24-atoms-b3lyp-6-31g(d)-orbitals-density", response, start)
        finally:
            done.set()
            for observer_result in observers:
                observer_result.result()
    receipts["responsiveness"].update(
        {
            "pollCount": len(latencies),
            "maximumPollSeconds": max(latencies),
            "meanPollSeconds": sum(latencies) / len(latencies),
            "explicitBusyReadResponsesRetried": len(busy_reads),
            "busyRetryBudgetSeconds": 2,
        }
    )
    assert max(latencies) < 5
    save()

    for name, method, spacing in [
        ("dense-caffeine-gasteiger", "gasteiger", 0.1),
        ("caffeine-xtb", "xtb-gfn2", 0.25),
    ]:
        print("Starting " + name, flush=True)
        now = time.monotonic()
        esp_job = submit(
            "/v1/electrostatics",
            {
                "sdfContent": fixture("caffeine"),
                "method": method,
                "surface": {"gridSpacing": spacing, "probeRadius": 0.9, "padding": 1.2},
            },
        )
        record(name, poll("/v1/electrostatics", esp_job), now)

    print("Starting interrupted aspirin calculation", flush=True)
    now = time.monotonic()
    interrupted = submit(
        "/v1/orbitals",
        {
            "sdfContent": fixture("aspirin"),
            "method": "b3lyp",
            "basis": "6-31g(d)",
            "orbitals": ["homo", "lumo", "density"],
        },
    )
    while call("/v1/orbitals/jobs/" + interrupted)["status"] != "running":
        assert time.monotonic() - now < 15
        time.sleep(0.05)
    time.sleep(1)
    first = row(interrupted)
    assert first["state"] == "running"
    worker = "moleculens-vps1-chemistry-1"
    subprocess.run(docker + ["update", "--restart=no", worker], check=True, capture_output=True)
    subprocess.run(docker + ["kill", worker], check=True, capture_output=True)
    state = json.loads(
        subprocess.check_output(
            docker + ["inspect", worker, "--format", "{{json .State}}"], text=True
        )
    )
    assert not state["Running"] and state["Pid"] == 0
    subprocess.run(docker + ["start", worker], check=True, capture_output=True)
    subprocess.run(
        docker + ["update", "--restart=unless-stopped", worker], check=True, capture_output=True
    )
    response = poll("/v1/orbitals", interrupted)
    record("aspirin-21-atoms-after-SIGKILL", response, now)
    second = row(interrupted)
    assert first["attempt"] == 1 and second["attempt"] == 2 and first["token"] != second["token"]
    receipts["interruption"] = {
        "jobId": interrupted,
        "before": first,
        "after": second,
        "reapedPID": state["Pid"],
        "exitCode": state["ExitCode"],
        "recovery": "manual container restart, then automatic fenced lease recovery",
    }
    receipts["passed"] = True
    save()
    print(json.dumps(receipts, indent=2))


if __name__ == "__main__":
    main()
