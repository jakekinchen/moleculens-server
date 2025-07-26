# Stub fastapi module for tests

from .testclient import TestClient


class FastAPI:
    def __init__(self, *args, **kwargs):
        pass

    def include_router(self, router, *args, **kwargs):
        # Stub method to include routers
        pass

    def post(self, *args, **kwargs):
        # Stub decorator for POST endpoints
        def decorator(func):
            return func

        return decorator

    def get(self, *args, **kwargs):
        # Stub decorator for GET endpoints
        def decorator(func):
            return func

        return decorator


# Stub APIRouter with decorator methods
class APIRouter:
    def __init__(self, *args, **kwargs):
        pass

    def get(self, *args, **kwargs):
        def decorator(func):
            return func

        return decorator

    def post(self, *args, **kwargs):
        def decorator(func):
            return func

        return decorator


# Stub HTTPException
class HTTPException(Exception):  # noqa: N818
    def __init__(self, status_code: int, detail: str):
        super().__init__(detail)
        self.status_code = status_code
        self.detail = detail


# Expose TestClient for import convenience
__all__ = ["FastAPI", "APIRouter", "HTTPException", "TestClient"]
