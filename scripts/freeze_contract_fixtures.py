"""Independent descriptor vectors and consumer requests; no service imports."""

import hashlib
import json
import struct
from pathlib import Path

ROOT = Path(__file__).parents[1]
sdf = (ROOT / "tests/contracts/water-v2000.mol").read_text()


def canonical(v):
    if isinstance(v, str):
        encoded = json.dumps(v, ensure_ascii=False, separators=(",", ":"))
        for char in "<>&\u2028\u2029":
            encoded = encoded.replace(char, f"\\u{ord(char):04x}")
        return "s" + encoded
    if isinstance(v, bool):
        return "b1" if v else "b0"
    if isinstance(v, (int, float)):
        return "f" + struct.pack(">d", float(v)).hex()
    if isinstance(v, list):
        return "[" + "".join(canonical(x) + "," for x in v) + "]"
    if isinstance(v, dict):
        return "{" + "".join(canonical(k) + ":" + canonical(v[k]) + "," for k in sorted(v)) + "}"
    if v is None:
        return "n"
    raise TypeError(v)


cases = [
    (
        "orbital",
        "psi4-2",
        {
            "sdf_content": sdf,
            "method": "scf",
            "basis": "sto-3g",
            "grid_spacing": 0.25,
            "isovalue": 0.02,
            "density_isovalue": 0.01,
            "orbitals": ["lumo", "homo", "density"],
            "charge": 0,
            "multiplicity": 1,
            "client_cache_key": "consumer-hint",
        },
        {
            "method": "scf",
            "basis": "sto-3g",
            "gridSpacing": 0.25,
            "isovalue": 0.02,
            "densityIsovalue": 0.01,
            "orbitals": ["density", "homo", "lumo"],
            "charge": 0,
            "multiplicity": 1,
        },
        "error",
        sdf,
    ),
    (
        "electrostatics",
        "esp-2",
        {
            "sdfContent": sdf,
            "method": "gasteiger",
            "no_com": True,
            "no_reorient": True,
            "surface": {"gridSpacing": 0.25, "probeRadius": 0.9, "padding": 1.2},
        },
        {
            "method": "gasteiger",
            "charge": 0,
            "multiplicity": 1,
            "gridSpacing": 0.25,
            "probeRadius": 0.9,
            "padding": 1.2,
            "units": "kcal/mol/e",
            "softeningEpsilon": 0.3,
            "clampPercentiles": [5, 95],
        },
        "error",
        sdf,
    ),
    (
        "conformer",
        "rdkit-conformer-2",
        {
            "molblock2d": sdf,
            "geomVersion": "etkdgv3",
            "waitMs": 2000,
            "seed": 42,
            "optimize": True,
            "forceField": "mmff",
            "numConformers": 1,
        },
        {
            "method": "etkdg_v3",
            "opt": "mmff",
            "geomVersion": "etkdgv3",
            "seed": 42,
            "maxAttempts": 10,
            "maxOptIters": 200,
            "addHs": False,
        },
        "preserve-original-planar",
        sdf,
    ),
    (
        "nomenclature",
        "openclatura-1",
        {"smiles": "CCO", "verifyOpsin": True},
        {"molblock": "", "verifyOpsin": True},
        "error",
        "CCO",
    ),
]
fixtures = []
for kind, compatibility, request, params, fallback, text in cases:
    descriptor = {
        "protocol": 1,
        "kind": kind,
        "compatibility": compatibility,
        "inputSha256": hashlib.sha256(text.encode()).hexdigest(),
        "params": params,
        "fallback": fallback,
    }
    fixtures.append(
        {
            "kind": kind,
            "request": request,
            "params": params,
            "inputSha256": descriptor["inputSha256"],
            "key": hashlib.sha256(canonical(descriptor).encode()).hexdigest(),
        }
    )
(ROOT / "tests/contracts/requests.json").write_text(json.dumps(fixtures, indent=2) + "\n")
