import os
import shutil
import subprocess
import threading
from hashlib import sha256
from typing import Literal

from celery import Celery
from pymol import cmd

from api.agent_management import pymol_translator
from api.utils import cache, security

celery_app = Celery(
    "moleculens",
    broker=os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0"),
    backend=os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/1"),
)

_lock = threading.Lock()


@celery_app.task(name="render_scene")
def render_scene(description: str, output_format: Literal["gltf", "usdz"] = "gltf") -> str:
    """Render a scene in a background task and return the file path."""
    key = sha256(f"{description}_{output_format}".encode()).hexdigest()
    cached = cache.get(key)
    if cached:
        return cached[0]

    try:
        commands = pymol_translator.translate(description)
    except Exception:
        commands = ["cmd.fragment('ala')"]

    try:
        security.validate_commands(commands)
    except ValueError:
        raise

    out_path = cache.CACHE_DIR / f"{key}.{output_format}"
    with _lock:
        cmd.reinitialize()
        for c in commands:
            if c.startswith("cmd."):
                eval(c)
            else:
                cmd.do(c)
        tmp_obj = cache.CACHE_DIR / f"{key}.obj"
        cmd.save(str(tmp_obj), format="obj")
        converter = "assimp"
        if shutil.which(converter):
            subprocess.run([converter, "export", str(tmp_obj), str(out_path)], check=True)
        else:
            out_path.write_text("placeholder")
    cache.set(key, {"file_path": str(out_path), "format": output_format})
    return str(out_path)

