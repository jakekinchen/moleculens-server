from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel
from typing import Literal, Dict, Any
from pathlib import Path
from hashlib import sha256
import tempfile
import subprocess
import json
import threading

from api.utils import llm, cache, security

from pymol import cmd

router = APIRouter(prefix="/render", tags=["Render"])

# Only mutate PyMOL from one request at a time
lock = threading.Lock()


class RenderRequest(BaseModel):
    description: str
    format: Literal["image", "model", "animation"] = "image"


def _output_path(key: str, fmt: str) -> Path:
    ext = {"image": "png", "model": "pdb", "animation": "mp4"}[fmt]
    return cache.CACHE_DIR / f"{key}.{ext}"


def _metadata() -> Dict[str, Any]:
    """Collect useful scene metadata from PyMOL (best-effort)."""
    data: Dict[str, Any] = {}
    try:
        data["camera"] = cmd.get_view()
    except Exception:
        pass
    try:
        data["center"] = cmd.centerofmass()
    except Exception:
        pass
    try:
        data["bbox"] = cmd.get_extent()
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

    # Convert description to PyMOL commands
    commands = llm.description_to_commands(req.description)
    # Fallback: if LLM could not translate description (e.g. no API key),
    # render a simple alanine fragment so that the PNG is not blank.
    if not commands:
        commands = [
            "cmd.fragment('ala')",
            "cmd.show('sticks')",
            "cmd.orient()",
        ]
    try:
        security.validate_commands(commands)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    out_path = _output_path(key, req.format)

    # Generate output under lock so concurrent requests don't clash
    with lock:
        cmd.reinitialize()
        for c in commands:
            # Supports both "cmd.fetch('1abc')" and "fetch 1abc" styles
            if c.startswith("cmd."):
                eval(c)
            else:
                cmd.do(c)

        if req.format == "image":
            cmd.png(str(out_path))
        elif req.format == "model":
            cmd.save(str(out_path))
        else:  # animation
            tmp_dir = tempfile.mkdtemp()
            frame = Path(tmp_dir) / "frame.png"
            cmd.png(str(frame))
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
    """Return a FileResponse or, for large payloads, a JSON pointer to static storage."""
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

    return FileResponse(path, media_type=media_types[fmt], headers={"X-Metadata": json.dumps(meta)})
