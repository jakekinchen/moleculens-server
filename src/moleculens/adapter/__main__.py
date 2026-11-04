"""Read a frozen request; write validated-result inputs and an integrity manifest.

The Go supervisor supplies a private output directory and no database credentials.
Each invocation exits after one attempt, including nomenclature. It owns no queue.
"""

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path


def provenance():
    packages = {}
    for name in ("numpy", "rdkit", "scikit-image", "openclatura", "py2opsin"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "unavailable"
    try:
        import psi4

        packages["psi4"] = psi4.__version__
    except ImportError:
        packages["psi4"] = ""
    source = hashlib.sha256()
    for path in sorted(Path(__file__).parents[1].rglob("*.py")):
        source.update(str(path.relative_to(Path(__file__).parents[1])).encode())
        source.update(path.read_bytes())
    packages["sourceSha256"] = source.hexdigest()
    for label, command in (("java", ["java", "-version"]), ("xtb", ["xtb", "--version"])):
        try:
            proc = subprocess.run(command, capture_output=True, text=True, timeout=5)
            packages[label] = (
                (proc.stdout + proc.stderr).strip()[:2000]
                if proc.returncode == 0
                else "unavailable"
            )
        except (OSError, subprocess.SubprocessError):
            packages[label] = "unavailable"
    packages.update(
        python=platform.python_version(),
        platform=platform.platform(),
        architecture=platform.machine(),
        build=os.environ.get("ENGINE_BUILD", "local-unqualified"),
    )
    return packages


def _json(path, obj):
    with path.open("x", encoding="utf-8") as f:
        json.dump(obj, f, allow_nan=False, sort_keys=True, separators=(",", ":"))
        f.flush()
        os.fsync(f.fileno())


def compute(request, output):
    p = request["params"]
    kind = request["kind"]
    if kind == "orbital":
        from moleculens.compute.psi4_runner import run_psi4_computation

        result = run_psi4_computation(
            request["input"],
            p["method"],
            p["basis"],
            p["gridSpacing"],
            p["isovalue"],
            p["orbitals"],
            output,
            output / "scratch",
            int(p["charge"]),
            int(p["multiplicity"]),
            p["densityIsovalue"],
        )
        from moleculens.api.schemas import ComputeResult

        raw = json.loads((output / "result.json").read_text())
        payload = ComputeResult.model_validate(raw).model_dump(by_alias=True)
    elif kind == "electrostatics":
        from moleculens.compute.electrostatics import run_electrostatics_computation

        result = run_electrostatics_computation(
            request["input"],
            int(p["charge"]),
            int(p["multiplicity"]),
            p["method"],
            p["gridSpacing"],
            p["probeRadius"],
            p["padding"],
            p["softeningEpsilon"],
            tuple(p["clampPercentiles"]),
            output,
        )
        from moleculens.api.schemas import ElectrostaticsMeta, ElectrostaticsResult

        meta = result.meta
        payload = ElectrostaticsResult(
            surface=result.surface,
            dipole=result.dipole,
            meta=ElectrostaticsMeta(
                chargeModel=meta["charge_model"],
                method=meta["method"],
                timingsMs=meta["timings_ms"],
                sdfHash=meta["sdf_hash"],
                molecularCharge=meta["molecular_charge"],
                multiplicity=meta["multiplicity"],
            ),
        ).model_dump(by_alias=True)
    elif kind == "conformer":
        from moleculens.compute.conformers import run_conformer_computation

        result = run_conformer_computation(
            request["input"],
            p["method"],
            p["opt"],
            int(p["maxAttempts"]),
            int(p["maxOptIters"]),
            p["addHs"],
            p["geomVersion"],
            output,
            int(p["seed"]),
        )
        from moleculens.api.schemas import ConformerMeta

        payload = {
            "sdf3d": (output / "structure.sdf").read_text(),
            "meta": ConformerMeta.model_validate(result.meta).model_dump(by_alias=True),
            "timingMs": result.meta["timings_ms"],
        }
    elif kind == "nomenclature":
        # This pure route has no database dependency; its response preserves
        # the existing explanation/mapping wire contract and typed errors.
        from moleculens.api.nomenclature_routes import explain
        from moleculens.api.schemas import NomenclatureRequest

        payload = explain(
            NomenclatureRequest(
                smiles=request["input"],
                molblock=p["molblock"] or None,
                verifyOpsin=p["verifyOpsin"],
            )
        ).model_dump(by_alias=True)
        result = None
    else:
        raise ValueError("Unsupported operation")
    _json(output / "response.json", payload)
    return ["response.json"] + (result.artifact_files if result else [])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe", type=Path)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--token")
    args = parser.parse_args()
    if args.probe:
        _json(args.probe, provenance())
        return
    start = time.monotonic()
    request = json.loads(args.request.read_text())
    if request["protocol"] != 1:
        raise ValueError("Unsupported adapter protocol")
    if hashlib.sha256(request["input"].encode()).hexdigest() != request["inputSha256"]:
        raise ValueError("Input digest mismatch")
    try:
        names = compute(request, args.output)
        files = []
        for name in sorted(set(names)):
            if Path(name).name != name:
                raise ValueError("Artifact name must be a basename")
            data = (args.output / name).read_bytes()
            files.append(
                {"name": name, "size": len(data), "sha256": hashlib.sha256(data).hexdigest()}
            )
        _json(
            args.output / "manifest.json",
            {
                "protocol": 1,
                "kind": request["kind"],
                "key": request["key"],
                "token": args.token,
                "compatibility": request["compatibility"],
                "inputSha256": request["inputSha256"],
                "provenance": provenance(),
                "files": files,
                "payload": "response.json",
                "elapsedMs": (time.monotonic() - start) * 1000,
            },
        )
    except Exception as exc:
        detail = getattr(exc, "detail", None)
        _json(
            args.output / "error.json",
            {
                "code": detail.get("code", "chemistry_failed")
                if isinstance(detail, dict)
                else "chemistry_failed",
                "httpStatus": getattr(exc, "status_code", 422),
                "message": (
                    str(detail.get("message", exc)) if isinstance(detail, dict) else str(exc)
                )[:1800],
            },
        )
        raise
    finally:
        if "psi4" in sys.modules:
            sys.modules["psi4"].core.clean()
            sys.modules["psi4"].core.clean_options()
        shutil.rmtree(args.output / "scratch", ignore_errors=True)


if __name__ == "__main__":
    main()
