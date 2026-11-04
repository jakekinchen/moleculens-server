"""Fetch the versioned native xTB release and verify its published SHA-256."""

import hashlib
import tarfile
import urllib.request
from pathlib import Path

URL = "https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz"
SHA256 = "62a8d18778286e815292ee53d76ce447daf460a4dea3782c0f25cbac7019b5df"
archive = Path("/tmp/xtb-release.tar.xz")
with urllib.request.urlopen(URL, timeout=60) as response:
    archive.write_bytes(response.read(100 * 1024 * 1024))
if hashlib.sha256(archive.read_bytes()).hexdigest() != SHA256:
    raise RuntimeError("Native xTB archive checksum mismatch")
with tarfile.open(archive) as package:
    package.extractall("/opt", filter="data")
Path("/opt/xtb-dist").rename("/opt/xtb")
archive.unlink()
