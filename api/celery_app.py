import shutil
import subprocess
import threading
from hashlib import sha256
from pathlib import Path
from typing import Any, Literal, Optional

from celery import Celery
from pymol import cmd

from api.pymol import pymol_translator
from api.pymol.pymol_security import validate_commands


class SimpleCache:
    """Simple in-memory cache with file storage."""

    def __init__(self):
        self.CACHE_DIR = Path("/tmp/moleculens_cache")
        self.CACHE_DIR.mkdir(exist_ok=True)
        self._cache: dict[str, Any] = {}

    def get(self, key: str) -> Optional[tuple]:
        """Get cached result."""
        return self._cache.get(key)

    def set(self, key: str, value: Any) -> None:
        """Set cached result."""
        self._cache[key] = (value["file_path"], value)


cache = SimpleCache()

celery_app = Celery(
    "moleculens",
    broker="memory://",
    backend="cache+memory://",
)

_lock = threading.Lock()


@celery_app.task(name="render_scene")
def render_scene(description: str, output_format: Literal["gltf", "usdz"] = "gltf") -> str:
    """Render a scene in a background task and return the file path.

    Note: cache.get() returns a tuple (file_path, metadata) when a cached
    result is found, or None if no cache entry exists.
    """
    key = sha256(f"{description}_{output_format}".encode()).hexdigest()
    cached = cache.get(key)
    if cached:
        return cached[0]  # cached is a tuple (file_path, metadata)

    try:
        commands = pymol_translator.translate(description)  # type: ignore[attr-defined]
    except Exception:
        commands = ["cmd.fragment('ala')"]

    try:
        validate_commands(commands)
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
