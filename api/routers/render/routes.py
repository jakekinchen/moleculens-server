from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel
from typing import Literal
from pathlib import Path
import tempfile
import json
from hashlib import sha256
import subprocess
from pymol import cmd
import threading

from api.utils import llm, cache, security

router = APIRouter(prefix="/render", tags=["Render"])

lock = threading.Lock()

class RenderRequest(BaseModel):
    description: str
    format: Literal["image", "model", "animation"] = "image"


def _output_path(key: str, fmt: str) -> Path:
    ext = {"image": "png", "model": "pdb", "animation": "mp4"}[fmt]
    return cache.CACHE_DIR / f"{key}.{ext}"


def _metadata() -> dict:
    return {
        "camera": cmd.get_view(),
        "center": cmd.centerofmass(),
        "bbox": cmd.get_extent()
    }


@router.post("")
async def render(req: RenderRequest):
    key = sha256(f"{req.description}_{req.format}".encode()).hexdigest()
    cached = cache.get(key)
    if cached:
        metadata = json.dumps({"cached": True})
        return FileResponse(cached, media_type="application/octet-stream", headers={"X-Metadata": metadata})

    commands = llm.description_to_commands(req.description)
    try:
        security.validate_commands(commands)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    out_path = _output_path(key, req.format)
    with lock:
        cmd.reinitialize()
        for c in commands:
            cmd.do(c)
        if req.format == "image":
            cmd.png(str(out_path))
        elif req.format == "model":
            cmd.save(str(out_path))
        else:  # animation
            tmp_dir = tempfile.mkdtemp()
            frame = Path(tmp_dir) / "frame.png"
            cmd.png(str(frame))
            subprocess.run(["ffmpeg", "-y", "-loop", "1", "-i", str(frame), "-t", "2", str(out_path)], check=True)

    cache.set(key, str(out_path))
    meta = json.dumps(_metadata())
    if out_path.stat().st_size > 25 * 1024 * 1024:
        url_path = f"/static/{out_path.name}"
        out_static = Path("api/static") / out_path.name
        out_path.rename(out_static)
        return JSONResponse({"url": url_path, "metadata": json.loads(meta)})

    return FileResponse(out_path, media_type="application/octet-stream", headers={"X-Metadata": meta})
