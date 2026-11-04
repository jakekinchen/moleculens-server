"""End-to-end native qualification against an explicitly supplied isolated API."""

import argparse
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


def call(base, path, payload=None):
    req = urllib.request.Request(
        base + path,
        data=json.dumps(payload).encode() if payload is not None else None,
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=15) as response:
        return json.load(response)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    water = (Path(__file__).parents[1] / "tests/contracts/water-v2000.mol").read_text()
    cases = [
        (
            "orbital",
            "/v1/orbitals",
            {
                "sdf_content": water,
                "method": "scf",
                "basis": "sto-3g",
                "grid_spacing": 0.25,
                "isovalue": 0.05,
                "density_isovalue": 0.01,
                "orbitals": ["homo", "lumo", "density"],
                "charge": 0,
                "multiplicity": 1,
            },
        ),
        (
            "electrostatics",
            "/v1/electrostatics",
            {
                "sdfContent": water,
                "method": "gasteiger",
                "surface": {"gridSpacing": 0.25, "probeRadius": 0.9, "padding": 1.2},
            },
        ),
        ("electrostatics-xtb", "/v1/electrostatics", {"sdfContent": water, "method": "xtb-gfn2"}),
        (
            "conformer",
            "/v1/conformers",
            {
                "molblock2d": water,
                "geomVersion": "etkdgv3",
                "seed": 42,
                "numConformers": 1,
                "optimize": True,
                "forceField": "mmff",
                "waitMs": 2000,
            },
        ),
    ]
    receipts = {"api": args.url, "cases": []}
    for name, path, payload in cases:
        start = time.monotonic()
        response = call(args.url, path + "/compute", payload)
        job = response["jobId"]
        states = [response["status"]]
        while response["status"] not in ("done", "error", "failed"):
            if time.monotonic() - start > 150:
                raise RuntimeError(f"{name} deadline")
            time.sleep(0.2)
            response = call(args.url, path + "/jobs/" + job)
            if states[-1] != response["status"]:
                states.append(response["status"])
        record = {
            "kind": name,
            "jobId": job,
            "cacheKey": response["cacheKey"],
            "states": states,
            "wallSeconds": time.monotonic() - start,
            "computeTimeMs": response.get("computeTimeMs"),
        }
        receipts["cases"].append(record)
        args.output.write_text(json.dumps(receipts, indent=2))
        assert response["status"] == "done", response
        if name == "orbital":
            result = response["result"]
            assert set(result["orbitals"]) == {"homo", "lumo"}
            assert result["density"]["vertexCount"] > 0
            assert result["orbitals"]["homo"]["energyEv"] < result["orbitals"]["lumo"]["energyEv"]
            record["energiesEv"] = {k: v["energyEv"] for k, v in result["orbitals"].items()}
            record["densityVertices"] = result["density"]["vertexCount"]
        elif name.startswith("electrostatics"):
            result = response["result"]
            assert result["surface"]["mesh"]["vertexCount"] > 0
            assert result["dipole"]["magnitudeDebye"] > 0
            record["chargeModel"] = result["meta"]["chargeModel"]
            record["dipoleDebye"] = result["dipole"]["magnitudeDebye"]
        else:
            assert "V2000" in response["sdf3d"]
            assert response["meta"]["opt"] == "mmff"
            record["quality"] = response["meta"]["quality"]
        cached = call(args.url, path + "/compute", payload)
        assert cached["jobId"] == job and cached["cached"] is True
        record["cacheReuseVerified"] = True
        args.output.write_text(json.dumps(receipts, indent=2))
    naming = call(args.url, "/v1/nomenclature/explain", {"smiles": "CCO", "verifyOpsin": True})
    assert naming["displayName"] and naming["verified"], naming
    receipts["nomenclature"] = {
        k: naming[k]
        for k in ("displayName", "verified", "opsinStatus", "engineVersion", "rdkitVersion")
    }
    args.output.write_text(json.dumps(receipts, indent=2))
    print(json.dumps(receipts, indent=2))


if __name__ == "__main__":
    main()
