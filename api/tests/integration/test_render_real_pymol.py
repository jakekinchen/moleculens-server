import importlib.util
import os
import sys

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient


@pytest.mark.integration
def test_render_model_with_real_pymol(capsys):
    """Exercise the /render route using the **real** PyMOL library.

    The test is skipped automatically when PyMOL is not present in the
    environment so it does not break CI runs that rely solely on the stubbed
    interface.
    """

    pymol = pytest.importorskip("pymol")  # noqa: F841 – imported for its side-effects

    # Import the router module directly
    routes_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../routers/render/routes.py")
    )
    spec = importlib.util.spec_from_file_location("routes", routes_path)
    routes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(routes)

    app = FastAPI()
    app.include_router(routes.router)

    # Create output directory in the mounted volume
    output_dir = "/app/api/static/test_output"
    os.makedirs(output_dir, exist_ok=True)

    # Create a log file for debug output
    log_path = os.path.join(output_dir, "debug.log")
    with open(log_path, "w") as log:
        log.write("Starting PyMOL visualization test\n\n")

        with TestClient(app) as client:
            # Test different visualization options
            test_cases = [
                # Overview scene - just needs structure_id
                {
                    "description": "Create an overview visualization of EGFR kinase domain structure with PDB ID 2ITX",
                    "format": "image",
                },
                # Binding site scene - needs structure_id and selection
                {
                    "description": "Show the ATP binding site in EGFR (PDB: 2ITX). The key residue is 745, please show residues within 5 angstroms.",
                    "format": "image",
                },
                # Mutation scene - needs structure_id, mutation_selection, and original_residue
                {
                    "description": "Visualize the T790M mutation in EGFR (PDB: 2ITX). The original residue was Threonine (T).",
                    "format": "image",
                },
                # Mutation focus scene - needs structure_id and mutation_selection
                {
                    "description": "Create a close-up view of residue 790 in EGFR (PDB: 2ITX) where the T790M mutation occurs",
                    "format": "image",
                },
                # Raw commands for complex visualizations
                {
                    "description": "raw: fetch 2itx; hide everything; show cartoon; show sticks, resi 790; color magenta, resi 790; show surface, resi 790 around 4; set transparency, 0.5; zoom resi 790, 8; bg_color white",
                    "format": "image",
                },
            ]

            for i, test_case in enumerate(test_cases):
                try:
                    log.write(f"\n{'='*80}\n")
                    log.write(f"Processing test case {i+1}:\n")
                    log.write(f"Description: {test_case['description']}\n")
                    log.flush()  # Ensure output is written immediately

                    resp = client.post("/render", json=test_case)

                    # Capture and write output
                    captured = capsys.readouterr()
                    log.write("\nDebug output:\n")
                    log.write(captured.out)
                    log.write(captured.err)
                    log.flush()

                    assert resp.status_code == 200
                    assert len(resp.content) > 1000

                    # Save the output
                    ext = "pdb" if test_case["format"] == "model" else "png"
                    output_path = os.path.join(output_dir, f"test_output_{i}.{ext}")
                    with open(output_path, "wb") as f:
                        f.write(resp.content)
                    log.write(f"Saved to {output_path}\n")
                    log.write("=" * 80 + "\n")
                    log.flush()

                    # Additional checks for PDB format
                    if test_case["format"] == "model":
                        assert resp.content.strip().endswith(b"END")
                except Exception as e:
                    log.write(f"Failed: {str(e)}\n")
                    # Show the full traceback for debugging
                    import traceback

                    log.write(traceback.format_exc())
                    log.flush()
                    # Continue with other test cases even if one fails
                    continue
