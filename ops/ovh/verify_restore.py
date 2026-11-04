"""Restore a fresh candidate dump into a disposable DB and verify retained artifacts."""

import argparse
import hashlib
import json
import os
import subprocess
import uuid
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--dump", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--legacy-manifest", type=Path, help="Optional verified historical file/hash inventory"
)
args = parser.parse_args()
assert os.geteuid() == 0, "Run as root on the new OVH host."
os.chdir("/opt/moleculens")
dc = ["docker", "compose", "--env-file", ".env", "-f", "compose.vps1.yml", "exec", "-T", "postgres"]


def pg(database, sql):
    return subprocess.check_output(
        dc + ["psql", "-U", "moleculens", "-d", database, "-Atc", sql], stdin=subprocess.DEVNULL
    )


database = "restorecheck_" + uuid.uuid4().hex[:16]
created = False
try:
    subprocess.run(
        dc + ["createdb", "-U", "moleculens", "-T", "template0", "-E", "UTF8", database],
        stdin=subprocess.DEVNULL,
        check=True,
    )
    created = True
    with args.dump.open("rb") as source:
        subprocess.run(
            dc
            + ["pg_restore", "-U", "moleculens", "-d", database, "--exit-on-error", "--no-owner"],
            stdin=source,
            check=True,
        )
    sql = "COPY (SELECT row_to_json(t)::text FROM moleculens_v2.jobs t ORDER BY id) TO STDOUT"
    original = pg("moleculens", sql)
    restored = pg(database, sql)
    assert original == restored, "Restored job records differ."
    legacy_tables = {}
    for table, order in [("jobs", "id"), ("cache_entries", "cache_key")]:
        query = (
            f"COPY (SELECT row_to_json(t)::text FROM public.{table} t ORDER BY {order}) TO STDOUT"
        )
        live, copy = pg("moleculens", query), pg(database, query)
        assert live == copy, f"Restored historical {table} records differ."
        legacy_tables[table] = {
            "rows": len(copy.splitlines()),
            "sha256": hashlib.sha256(copy).hexdigest(),
        }
    volume = Path(
        subprocess.check_output(
            [
                "docker",
                "volume",
                "inspect",
                "moleculens-vps1_artifacts",
                "--format",
                "{{.Mountpoint}}",
            ],
            text=True,
        ).strip()
    ).resolve()
    manifests = pg(
        database,
        "SELECT manifest::text FROM moleculens_v2.jobs WHERE state='done' AND available ORDER BY id",
    )
    files_checked = 0
    jobs_checked = 0
    for line in manifests.splitlines():
        manifest = json.loads(line)
        attempt = (volume / manifest["key"] / manifest["token"]).resolve()
        assert attempt.is_relative_to(volume), "Manifest escapes artifact root."
        for entry in manifest["files"]:
            assert Path(entry["name"]).name == entry["name"]
            path = attempt / entry["name"]
            assert not path.is_symlink() and path.is_file()
            assert path.stat().st_size == entry["size"]
            with path.open("rb") as data:
                assert hashlib.file_digest(data, "sha256").hexdigest() == entry["sha256"]
            files_checked += 1
        jobs_checked += 1
    assert jobs_checked > 0 and files_checked > 0
    historical_files = 0
    if args.legacy_manifest:
        inventory = json.loads(args.legacy_manifest.read_text())
        legacy = Path(
            subprocess.check_output(
                [
                    "docker",
                    "volume",
                    "inspect",
                    "moleculens-vps1_legacy",
                    "--format",
                    "{{.Mountpoint}}",
                ],
                text=True,
            ).strip()
        ).resolve()
        expected = {entry["name"]: entry for entry in inventory["files"]}
        actual = set()
        for path in legacy.rglob("*"):
            assert not path.is_symlink(), "Historical artifact symlink prohibited."
            if path.is_dir():
                continue
            assert path.is_file() and path.resolve().is_relative_to(legacy)
            relative = str(path.relative_to(legacy))
            entry = expected[relative]
            assert path.stat().st_size == entry["size"]
            with path.open("rb") as source:
                assert hashlib.file_digest(source, "sha256").hexdigest() == entry["sha256"]
            actual.add(relative)
        assert actual == set(expected), "Historical file inventory differs."
        historical_files = len(actual)
    with args.dump.open("rb") as data:
        dump_hash = hashlib.file_digest(data, "sha256").hexdigest()
    receipt = {
        "dump": str(args.dump),
        "dumpSHA256": dump_hash,
        "jobRowsByteIdentical": True,
        "jobRowsSHA256": hashlib.sha256(original).hexdigest(),
        "completedJobsVerified": jobs_checked,
        "artifactFilesSHA256Verified": files_checked,
        "scope": "Logical dump restored into a separate disposable database; retained artifact bytes match restored manifests. Full OVH snapshot restore is not tested.",
    }
    receipt["historicalTables"] = legacy_tables
    receipt["historicalArtifactFilesSHA256Verified"] = historical_files
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
finally:
    if created:
        subprocess.run(
            dc + ["dropdb", "-U", "moleculens", database], stdin=subprocess.DEVNULL, check=True
        )
