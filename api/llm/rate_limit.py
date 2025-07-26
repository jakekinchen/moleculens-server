import time
from collections import defaultdict

from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import JSONResponse


class RateLimitMiddleware(BaseHTTPMiddleware):
    """Simple in-memory rate limiting middleware."""

    def __init__(self, app, per_ip: int = 5, global_limit: int = 30):
        super().__init__(app)
        self.per_ip = per_ip
        self.global_limit = global_limit
        self._ip_counts: dict[str, tuple[int, float]] = defaultdict(lambda: (0, 0.0))
        self._global_count = (0, 0.0)

    def _get_count(self, key: str, is_global: bool = False) -> int:
        """Get current count, cleaning expired entries."""
        current_time = time.time()

        if is_global:
            count, timestamp = self._global_count
            if current_time - timestamp > 60:  # 60 second window
                self._global_count = (1, current_time)
                return 1
            else:
                self._global_count = (count + 1, timestamp)
                return count + 1
        else:
            count, timestamp = self._ip_counts[key]
            if current_time - timestamp > 60:  # 60 second window
                self._ip_counts[key] = (1, current_time)
                return 1
            else:
                self._ip_counts[key] = (count + 1, timestamp)
                return count + 1

    async def dispatch(self, request: Request, call_next):
        ip = request.client.host if request.client else "unknown"

        ip_count = self._get_count(ip)
        global_count = self._get_count("global", is_global=True)

        if ip_count > self.per_ip or global_count > self.global_limit:
            return JSONResponse(status_code=429, content={"detail": "Rate limit exceeded"})

        response = await call_next(request)
        return response
