"""Run as root inside the disposable VM; measure actual cgroups and VM headroom."""

import argparse
import json
import os
import subprocess
import time
from contextlib import suppress
from pathlib import Path


def command(*args):
    return subprocess.check_output(args, text=True).strip()


def tree_bytes(path):
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            with suppress(FileNotFoundError):
                total += (Path(root) / name).stat().st_size
    return total


def measure(output, stop, project):
    names = command(
        "docker", "ps", "-q", "--filter", f"label=com.docker.compose.project={project}"
    ).splitlines()
    containers = json.loads(command("docker", "inspect", *names))
    groups = {}
    for item in containers:
        pid = item["State"]["Pid"]
        suffix = Path(f"/proc/{pid}/cgroup").read_text().strip().split(":", 2)[2]
        groups[item["Name"].lstrip("/")] = Path("/sys/fs/cgroup") / suffix.lstrip("/")
    volumes = {
        name: command(
            "docker", "volume", "inspect", f"{project}_{name}", "--format", "{{.Mountpoint}}"
        )
        for name in ["artifacts", "postgres"]
    }
    info = {
        line.split(":")[0]: int(line.split()[1]) * 1024
        for line in Path("/proc/meminfo").read_text().splitlines()
    }
    receipt = {
        "os": Path("/etc/os-release").read_text(),
        "cpus": os.cpu_count(),
        "memoryTotalBytes": info["MemTotal"],
        "swapTotalBytes": info["SwapTotal"],
        "sampleIntervalSeconds": 0.5,
        "minimumMemAvailableBytes": info["MemAvailable"],
        "maximumSimultaneousCgroupCurrentBytes": 0,
        "containerPeakBytes": {},
        "containerOOMEvents": {},
        "volumePeakLogicalBytes": {},
        "maximumSingleAttemptLogicalBytes": 0,
        "samples": 0,
    }
    start = time.monotonic()
    while not stop.exists() and time.monotonic() - start < 1200:
        info = {
            line.split(":")[0]: int(line.split()[1]) * 1024
            for line in Path("/proc/meminfo").read_text().splitlines()
        }
        receipt["minimumMemAvailableBytes"] = min(
            receipt["minimumMemAvailableBytes"], info["MemAvailable"]
        )
        current = 0
        for name, group in groups.items():
            try:
                current += int((group / "memory.current").read_text())
                peak = int((group / "memory.peak").read_text())
                receipt["containerPeakBytes"][name] = max(
                    receipt["containerPeakBytes"].get(name, 0), peak
                )
                events = dict(
                    line.split() for line in (group / "memory.events").read_text().splitlines()
                )
                receipt["containerOOMEvents"][name] = max(
                    receipt["containerOOMEvents"].get(name, 0),
                    int(events["oom"]),
                    int(events["oom_kill"]),
                )
            except FileNotFoundError:
                pass
        receipt["maximumSimultaneousCgroupCurrentBytes"] = max(
            receipt["maximumSimultaneousCgroupCurrentBytes"], current
        )
        for name, path in volumes.items():
            size = tree_bytes(path)
            receipt["volumePeakLogicalBytes"][name] = max(
                receipt["volumePeakLogicalBytes"].get(name, 0), size
            )
        for attempt in Path(volumes["artifacts"]).glob("*/*"):
            if attempt.is_dir():
                receipt["maximumSingleAttemptLogicalBytes"] = max(
                    receipt["maximumSingleAttemptLogicalBytes"], tree_bytes(attempt)
                )
        receipt["samples"] += 1
        receipt["elapsedSeconds"] = time.monotonic() - start
        output.write_text(json.dumps(receipt, indent=2) + "\n")
        time.sleep(0.5)
    receipt["diskUsageBytes"] = command("df", "-B1", "/", "/var/lib/docker")
    receipt["finalSwapUsedBytes"] = info["SwapTotal"] - info["SwapFree"]
    output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stop", type=Path, required=True)
    parser.add_argument("--project", default="moleculens-vps1")
    args = parser.parse_args()
    measure(args.output, args.stop, args.project)
