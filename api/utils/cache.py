import os
from pathlib import Path
from typing import Optional
import redis

REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

redis_client = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
CACHE_DIR = Path("/srv/cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def get(key: str) -> Optional[str]:
    value = redis_client.get(key)
    if value:
        path = Path(value.decode())
        if path.exists():
            return str(path)
    return None


def set(key: str, file_path: str):
    redis_client.set(key, file_path)

