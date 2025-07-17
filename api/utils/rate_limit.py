import redis
from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import JSONResponse
import os

REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))

class RateLimitMiddleware(BaseHTTPMiddleware):
    def __init__(self, app, per_ip: int = 5, global_limit: int = 30):
        super().__init__(app)
        self.redis = redis.Redis(host=REDIS_HOST, port=REDIS_PORT)
        self.per_ip = per_ip
        self.global_limit = global_limit

    async def dispatch(self, request: Request, call_next):
        ip = request.client.host
        ip_key = f"rl:ip:{ip}"
        global_key = "rl:global"
        pipe = self.redis.pipeline()
        pipe.incr(ip_key, 1)
        pipe.expire(ip_key, 60)
        pipe.incr(global_key, 1)
        pipe.expire(global_key, 60)
        ip_count, _, global_count, _ = pipe.execute()
        if int(ip_count) > self.per_ip or int(global_count) > self.global_limit:
            return JSONResponse(status_code=429, content={"detail": "Rate limit exceeded"})
        response = await call_next(request)
        return response
