import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

try:
    import redis  # type: ignore
except ImportError:  # pragma: no cover – fallback stub for test environments

    class _DummyRedis:
        _store: dict = {}

        def get(self, key):
            value = self._store.get(key)
            return value.encode() if isinstance(value, str) else value

        def set(self, key, value):
            self._store[key] = value

    redis = type("redis", (), {"Redis": lambda *args, **kwargs: _DummyRedis()})  # type: ignore

REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)

# Determine writable cache directory
_DEFAULT_CACHE_DIR = Path("/srv/cache")
try:
    _DEFAULT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
except (OSError, PermissionError):
    import tempfile

    _DEFAULT_CACHE_DIR = Path(tempfile.gettempdir()) / "sci_vis_cache"
    _DEFAULT_CACHE_DIR.mkdir(parents=True, exist_ok=True)

CACHE_DIR = _DEFAULT_CACHE_DIR


def get(key: str) -> Optional[Tuple[str, Dict[str, Any]]]:
    """Retrieve cached (file_path, metadata) tuple if it exists and the file is still present."""
    raw = redis_client.get(key)
    if not raw:
        return None

    try:
        meta: Dict[str, Any] = json.loads(raw)
    except (json.JSONDecodeError, TypeError):
        return None

    path = Path(meta.get("file_path", ""))
    if not path.exists():
        return None

    return str(path), meta


def set(key: str, meta: Dict[str, Any]) -> None:
    """Cache metadata under the provided key."""
    redis_client.set(key, json.dumps(meta))
