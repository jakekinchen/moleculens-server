"""Fail closed on missing/corrupt outputs; availability is not execution state."""

import base64
import gzip
import io
import json
from pathlib import Path

import numpy as np

from moleculens.compute.sdf_parser import parse_sdf

MAX_RESULT = 32 * 1024 * 1024


def validate_mesh(mesh: dict, potential: bool = False) -> None:
    n, t = mesh["vertexCount"], mesh["triangleCount"]
    if (
        not isinstance(n, int)
        or not isinstance(t, int)
        or not 0 < n <= 2_000_000
        or not 0 < t <= 4_000_000
    ):
        raise ValueError("Invalid mesh counts")
    counts = {"vertices": n * 3, "normals": n * 3, "indices": t * 3}
    if potential:
        counts["potential"] = n
    for name, count in counts.items():
        encoded = mesh[name]
        if not isinstance(encoded, str) or len(encoded) > MAX_RESULT:
            raise ValueError("Invalid mesh buffer size")
        compressed = base64.b64decode(encoded, validate=True)
        with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as stream:
            raw = stream.read(count * 4 + 1)
        if len(raw) != count * 4:
            raise ValueError("Mesh buffer length disagrees with counts")
        values = np.frombuffer(raw, dtype="<u4" if name == "indices" else "<f4")
        if name == "indices":
            if np.any(values >= n):
                raise ValueError("Mesh index out of range")
        elif not np.isfinite(values).all():
            raise ValueError("Non-finite mesh value")


def load_validated(kind: str, root: Path, meta: dict, params: dict | None = None) -> dict:
    params = params or {}
    if kind == "conformer":
        sdf = (root / "structure.sdf").read_text()
        parse_sdf(sdf)
        required = {"quality", "method", "opt", "geom_version", "embed_status", "opt_status"}
        if not required <= meta.keys():
            raise ValueError("Conformer metadata incomplete")
        return {"sdf3d": sdf}
    path = root / ("result.json" if kind == "orbital" else "electrostatics_result.json")
    if path.is_symlink() or path.stat().st_size > MAX_RESULT:
        raise ValueError("Invalid result file")
    data = json.loads(
        path.read_text(), parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value))
    )
    if not isinstance(data, dict) or not isinstance(data.get("meta"), dict) or not data["meta"]:
        raise ValueError("Result metadata absent")
    if kind == "orbital":
        if not {"method", "basis", "gridSpacingAngstrom", "psi4NumThreads"} <= data["meta"].keys():
            raise ValueError("Orbital metadata incomplete")
        orbs = data["orbitals"]
        if not isinstance(orbs, dict) or (not orbs and data.get("density") is None):
            raise ValueError("Orbital result empty")
        for name in params.get("orbitals", []):
            if name != "density" and name not in orbs:
                raise ValueError("Requested orbital absent")
        for orb in orbs.values():
            if not np.isfinite(orb["energyEv"]) or not np.isfinite(orb["isovalue"]):
                raise ValueError("Invalid orbital scalar")
            for sign in ("positive", "negative"):
                if orb.get(sign) is not None:
                    validate_mesh(orb[sign])
        if data.get("density") is not None:
            validate_mesh(data["density"])
    else:
        validate_mesh(data["surface"]["mesh"], potential=True)
        dipole = data["dipole"]
        for name in ("origin_angstrom", "vector_debye"):
            v = np.asarray(dipole[name])
            if v.shape != (3,) or not np.isfinite(v).all():
                raise ValueError("Invalid dipole")
        if not np.isclose(np.linalg.norm(dipole["vector_debye"]), dipole["magnitude_debye"]):
            raise ValueError("Dipole magnitude mismatch")
        if not {"method", "charge_model", "timings_ms", "sdf_hash"} <= data["meta"].keys():
            raise ValueError("ESP metadata incomplete")
    return data
