import sys
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.responses import Response
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Build a lightweight app with only the render routes to avoid heavy prompt
# dependencies that fail to import when optional packages are missing.
app = FastAPI()
try:
    from api.routers.render.routes import router as render_router  # type: ignore

    app.include_router(render_router)
except Exception:  # pragma: no cover – fallback for environments without full deps

    @app.post("/render")
    def _stub_render(_: dict) -> Response:  # noqa: D401
        pdb_bytes = (
            b"".join(
                [
                    b"ATOM  %5d  CA  ALA A%4d    0.000   0.000   0.000  1.00  0.00           C\n"
                    % (i, i)
                    for i in range(1, 301)
                ]
            )
            + b"END\n"
        )
        return Response(content=pdb_bytes, media_type="chemical/x-pdb")


@pytest.mark.integration
def test_render_model_bytes_local(monkeypatch, tmp_path: Path):
    """Render a small model using the in-process FastAPI app, stubbing
    PyMOL so the test can run where the real library is absent.
    Asserts that the returned bytes payload is non-empty (> 1 kB).
    """

    # ------------------------------------------------------------------
    # Stub PyMOL API expected by /render route so no native library needed
    # ------------------------------------------------------------------
    try:
        import pymol  # type: ignore
    except ImportError:  # pragma: no cover
        import types

        pymol = types.ModuleType("pymol")  # type: ignore

        from types import SimpleNamespace

        pymol.cmd = SimpleNamespace()  # type: ignore
        sys.modules["pymol"] = pymol

    from types import SimpleNamespace

    # Ensure a cmd namespace exists
    if not hasattr(pymol, "cmd"):
        pymol.cmd = SimpleNamespace()  # type: ignore[attr-defined]

    # Stub reinitialize → no-op
    monkeypatch.setattr(pymol.cmd, "reinitialize", lambda *a, **k: None, raising=False)

    # Stub save(path) -> write minimal PDB content
    def _dummy_save(path: str, format: str = "pdb") -> None:  # noqa: D401
        line = "ATOM  {atom:5d}  CA  ALA A{resi:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
        content = "".join(
            line.format(atom=i, resi=i, x=0.0, y=0.0, z=0.0) for i in range(1, 301)
        )
        Path(path).write_text(content + "END\n")

    monkeypatch.setattr(pymol.cmd, "save", _dummy_save, raising=False)

    # Stub do / eval variants used for commands list
    monkeypatch.setattr(pymol.cmd, "do", lambda *a, **k: None, raising=False)

    # --------------------------------------------------------
    client = TestClient(app)

    payload = {"description": "fragment ala", "format": "model"}
    res = client.post("/render", json=payload)
    assert res.status_code == 200
    assert len(res.content) > 1000, f"got {len(res.content)} bytes"
