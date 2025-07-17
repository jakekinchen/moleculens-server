import os
import json
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import redis

REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
CACHE_DIR = Path("/srv/cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)


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

