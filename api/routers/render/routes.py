from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import Optional, Dict, Any
from pathlib import Path
import json
import hashlib
import subprocess
import tempfile
import threading

from api.utils import llm, cache, security

router = APIRouter(prefix="/render", tags=["Render"])

lock = threading.Lock()

class RenderRequest(BaseModel):
    description: str
    format: Optional[str] = "image"


def _output_path(digest: str, fmt: str) -> Path:
    ext = {"image": "png", "model": "pdb", "animation": "mp4"}[fmt]
    return cache.CACHE_DIR / f"{digest}.{ext}"


def _get_metadata(cmd_module) -> Dict[str, Any]:
    meta: Dict[str, Any] = {}
    try:
        ext = cmd_module.get_extent()
        meta["bbox"] = ext
        meta["center"] = [(ext[0][i] + ext[1][i]) / 2.0 for i in range(3)]
    except Exception:
        pass
    try:
        meta["camera"] = cmd_module.get_view()
    except Exception:
        pass
    return meta


@router.post("/render")
async def render(req: RenderRequest):
    fmt = (req.format or "image").lower()
    if fmt not in {"image", "model", "animation"}:
        raise HTTPException(status_code=400, detail="Invalid format")
    digest = hashlib.sha256((req.description + fmt).encode()).hexdigest()
    cached = cache.get(digest)
    if cached:
        path, meta = cached
        return _build_response(path, fmt, meta)

    commands = await llm.description_to_commands(req.description)
    try:
        security.validate_commands(commands)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    from pymol import cmd

    out_path = _output_path(digest, fmt)
    with lock:
        cmd.reinitialize()
        for c in commands:
            eval(f"cmd.{c}")
        if fmt == "image":
            cmd.png(str(out_path))
        elif fmt == "model":
            cmd.save(str(out_path))
        else:
            tmp = tempfile.mkdtemp()
            frame = Path(tmp) / "frame.png"
            cmd.png(str(frame))
            subprocess.run([
                "ffmpeg",
                "-y",
                "-loop",
                "1",
                "-i",
                str(frame),
                "-t",
                "1",
                str(out_path),
            ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        meta = _get_metadata(cmd)

    meta.update({"file_path": str(out_path), "format": fmt})
    cache.set(digest, meta)
    return _build_response(str(out_path), fmt, meta)


def _build_response(path: str, fmt: str, meta: Dict[str, Any]) -> StreamingResponse:
    media = {
        "image": "image/png",
        "model": "chemical/x-pdb",
        "animation": "video/mp4",
    }[fmt]
    return StreamingResponse(open(path, "rb"), media_type=media, headers={"X-Meta": json.dumps(meta)})
