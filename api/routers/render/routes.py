import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import threading
import types
from hashlib import sha256
from pathlib import Path
from typing import Any, Dict, Literal

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel

try:
    from api.agent_management import pymol_translator
except Exception:  # pragma: no cover - fallback for test loaders
    translator_path = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), "../../agent_management/pymol_translator.py"
        )
    )
    package_root = os.path.dirname(translator_path)

    api_root = os.path.dirname(os.path.dirname(translator_path))
    api_pkg = sys.modules.setdefault("api", types.ModuleType("api"))
    api_pkg.__path__ = [api_root]
    agent_pkg = sys.modules.setdefault(
        "api.agent_management", types.ModuleType("api.agent_management")
    )
    agent_pkg.__path__ = [package_root]
    api_pkg.agent_management = agent_pkg  # type: ignore[attr-defined]

    spec = importlib.util.spec_from_file_location(
        "api.agent_management.pymol_translator", translator_path
    )
    assert spec is not None, "Could not create module spec for pymol_translator"
    pymol_translator = importlib.util.module_from_spec(spec)
    assert spec.loader is not None, "Module spec has no loader"
    spec.loader.exec_module(pymol_translator)

from api.utils import cache, security

# ---------------------------------------------------------------------------
# PyMOL integration
# ---------------------------------------------------------------------------
# We want the router to use the real PyMOL library when it is installed but still
# be import-able in environments (e.g. CI) where the native shared library is
# missing.  The existing integration tests patch `sys.modules["pymol"]` with a
# stub *before* importing this router, so the block below also behaves nicely in
# those cases.

try:
    # Set PyMOL environment variables for headless mode
    os.environ["PYMOL_QUIET"] = "1"
    os.environ["PYMOL_HEADLESS"] = "1"

    import pymol  # type: ignore

    # Launch PyMOL in quiet, headless mode exactly once.  This call is a no-op
    # on the stub object injected by the tests because the attribute will be
    # absent.
    if hasattr(pymol, "finish_launching"):
        pymol.finish_launching(["pymol", "-cq"])  # c = command-line, q = quiet

    pymol_cmd = pymol.cmd  # type: ignore[attr-defined]

except Exception:  # pragma: no cover – fallback when PyMOL truly unavailable
    # Create a minimal stub so the module still loads; the integration tests
    # will monkey-patch the attributes they need at runtime.
    from types import SimpleNamespace

    pymol_cmd = SimpleNamespace()  # type: ignore

    # Define stub methods for CI environments
    STUB_METHODS = [
        "reinitialize",
        "do",
        "png",
        "save",
        "get_view",
        "centerofmass",
        "get_extent",
    ]

    def noop(*args, **kwargs):
        return None

    # Set all stub methods
    for method in STUB_METHODS:
        setattr(pymol_cmd, method, noop)

router = APIRouter(prefix="/render", tags=["Render"])

# Only mutate PyMOL from one request at a time
lock = threading.Lock()


class RenderRequest(BaseModel):
    description: str
    format: Literal["image", "model", "animation"] = "image"
    # NEW: Advanced rendering options
    transparent_background: bool = False
    ray_trace: bool = True
    resolution: tuple[int, int] = (1920, 1080)
    dpi: int = 300
    ray_trace_mode: Literal["default", "cartoon_outline", "bw", "poster"] = "default"
    antialias: bool = True
    ray_shadow: bool = True
    depth_cue: bool = True
    background_color: str = "white"


def _output_path(key: str, fmt: str) -> Path:
    ext = {"image": "png", "model": "pdb", "animation": "mp4"}[fmt]
    return cache.CACHE_DIR / f"{key}.{ext}"


def _metadata() -> Dict[str, Any]:
    """Collect useful scene metadata from PyMOL (best-effort)."""
    data: Dict[str, Any] = {}
    try:
        data["camera"] = pymol_cmd.get_view()
    except Exception:
        pass
    try:
        data["center"] = pymol_cmd.centerofmass()
    except Exception:
        pass
    try:
        data["bbox"] = pymol_cmd.get_extent()
    except Exception:
        pass
    return data


@router.post("")
async def render(req: RenderRequest):
    key = sha256(f"{req.description}_{req.format}".encode()).hexdigest()

    # Serve from cache when possible
    cached = cache.get(key)
    if cached:
        path, meta = cached
        meta["cached"] = True
        return _build_response(path, req.format, meta)

    # Convert description to PyMOL commands via translator
    try:
        commands, llm_rendering_opts = pymol_translator.translate_with_options(req.description)  # type: ignore[attr-defined]
        # Merge LLM rendering options with request options (request takes precedence)
        merged_opts = {**llm_rendering_opts}
        for key, value in req.model_dump().items():
            if key in ["description", "format"]:
                continue
            if hasattr(req, key):
                merged_opts[key] = value
    except Exception:
        commands = []
        merged_opts = req.model_dump()
        merged_opts.pop("description", None)
        merged_opts.pop("format", None)

    # Fallback: if LLM could not translate description (e.g. no API key),
    # render a simple structure so that the PNG is not blank.
    if not commands:
        commands = [
            "cmd.fetch('1ubq')",  # Ubiquitin
            "cmd.hide('everything')",
            "cmd.show('cartoon')",
            "cmd.color('red', 'chain A')",
            "cmd.color('cyan', 'chain B')",
            "cmd.orient()",
            "cmd.label('name ca', 'resi')",  # label residues
        ]
    try:
        security.validate_commands(commands)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    out_path = _output_path(key, req.format)

    # Generate output under lock so concurrent requests don't clash
    with lock:
        pymol_cmd.reinitialize()
        for c in commands:
            # Supports both "cmd.fetch('1abc')" and "fetch 1abc" styles
            if c.startswith("cmd."):
                eval(c, {"cmd": pymol_cmd})
            else:
                pymol_cmd.do(c)

        # Apply advanced rendering settings using merged options
        if merged_opts.get("transparent_background", False):
            pymol_cmd.set("ray_opaque_background", 0)

        # Set background color
        pymol_cmd.do(f"bg_color {merged_opts.get('background_color', 'white')}")

        if merged_opts.get("ray_trace", True):
            # Configure ray-tracing quality
            pymol_cmd.set("antialias", 1 if merged_opts.get("antialias", True) else 0)
            pymol_cmd.set("ray_shadow", 1 if merged_opts.get("ray_shadow", True) else 0)
            pymol_cmd.set("depth_cue", 1 if merged_opts.get("depth_cue", True) else 0)

            # Set ray-trace mode
            ray_mode = merged_opts.get("ray_trace_mode", "default")
            if ray_mode == "cartoon_outline":
                pymol_cmd.set("ray_trace_mode", 3)
            elif ray_mode == "bw":
                pymol_cmd.set("ray_trace_mode", 2)
            elif ray_mode == "poster":
                pymol_cmd.set("ray_trace_mode", 1)

            # Ray-trace at specified resolution
            resolution = merged_opts.get("resolution", (1920, 1080))
            pymol_cmd.ray(resolution[0], resolution[1])

        if req.format == "image":
            # Use enhanced PNG export with DPI and ray-tracing
            dpi = merged_opts.get("dpi", 300)
            ray_trace = merged_opts.get("ray_trace", True)
            pymol_cmd.png(str(out_path), dpi=dpi, ray=1 if ray_trace else 0)
        elif req.format == "model":
            pymol_cmd.save(str(out_path))
        else:  # animation
            tmp_dir = tempfile.mkdtemp()
            frame = Path(tmp_dir) / "frame.png"
            # Apply same rendering settings for animation frames
            if merged_opts.get("ray_trace", True):
                resolution = merged_opts.get("resolution", (1920, 1080))
                pymol_cmd.ray(resolution[0], resolution[1])
            dpi = merged_opts.get("dpi", 300)
            ray_trace = merged_opts.get("ray_trace", True)
            pymol_cmd.png(str(frame), dpi=dpi, ray=1 if ray_trace else 0)
            subprocess.run(
                [
                    "ffmpeg",
                    "-y",
                    "-loop",
                    "1",
                    "-i",
                    str(frame),
                    "-t",
                    "2",
                    str(out_path),
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )

    meta = _metadata()
    meta.update({"file_path": str(out_path), "format": req.format})
    cache.set(key, meta)

    return _build_response(str(out_path), req.format, meta)


def _build_response(path: str, fmt: str, meta: Dict[str, Any]):
    """Return a FileResponse or, for large payloads, a JSON pointer to static
    storage."""
    media_types = {
        "image": "image/png",
        "model": "chemical/x-pdb",
        "animation": "video/mp4",
    }

    file_size = Path(path).stat().st_size
    if file_size > 25 * 1024 * 1024:  # 25 MB
        static_dir = Path("api/static")
        static_dir.mkdir(parents=True, exist_ok=True)
        static_path = static_dir / Path(path).name
        Path(path).rename(static_path)
        return JSONResponse({"url": f"/static/{static_path.name}", "metadata": meta})

    return FileResponse(
        path, media_type=media_types[fmt], headers={"X-Metadata": json.dumps(meta)}
    )
