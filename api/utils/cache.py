import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Protocol, Tuple, Union, runtime_checkable


@runtime_checkable
class RedisLike(Protocol):
    """Protocol for Redis-like cache backends."""

    def get(self, key: str) -> Optional[Union[str, bytes]]:
        """Get a value by key."""
        ...

    def set(self, key: str, value: str) -> None:
        """Set a value by key."""
        ...

    def ping(self) -> bool:
        """Test connection."""
        ...


# Simple in-memory cache implementation
class _DummyRedis:
    _store: dict = {}

    def get(self, key):
        value = self._store.get(key)
        return value.encode() if isinstance(value, str) else value

    def set(self, key, value):
        self._store[key] = value

    def ping(self) -> bool:
        return True


# Only try to use real Redis if explicitly configured
REDIS_HOST = os.getenv("REDIS_HOST")
REDIS_PORT = int(os.getenv("REDIS_PORT", "6379"))

if REDIS_HOST:
    try:
        import redis

        redis_client: RedisLike = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)  # type: ignore[assignment]
        # Test the connection
        redis_client.ping()
    except (ImportError, Exception):
        redis_client = _DummyRedis()
else:
    redis_client = _DummyRedis()

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
    """Retrieve cached (file_path, metadata) tuple if it exists and the file is
    still present."""
    raw = redis_client.get(key)
    if not raw:
        return None

    try:
        # Handle both bytes and str return types from Redis
        if isinstance(raw, bytes):
            raw_str = raw.decode()
        elif isinstance(raw, str):
            raw_str = raw
        else:
            return None
        meta: Dict[str, Any] = json.loads(raw_str)
    except (json.JSONDecodeError, TypeError, AttributeError):
        return None

    path = Path(meta.get("file_path", ""))
    if not path.exists():
        return None

    return str(path), meta


def set(key: str, meta: Dict[str, Any]) -> None:
    """Cache metadata under the provided key."""
    redis_client.set(key, json.dumps(meta))
