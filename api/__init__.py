"""API package root.

This file sets up import aliases so code/tests that expect top-level
`agent_management` and `llm_service` modules can still work even though the
actual implementation lives under `api.agent_management`.
"""

import importlib
import sys

# Expose `agent_management` package at top-level
_agent_mgmt_pkg = importlib.import_module("api.agent_management")
sys.modules.setdefault("agent_management", _agent_mgmt_pkg)

# Also expose `llm_service` (frequently imported directly) at top-level
sys.modules.setdefault(
    "llm_service", importlib.import_module("api.agent_management.llm_service")
)

# Backwards-compatibility single-file modules referenced in some tests

sys.modules.setdefault(
    "geometry_agent", importlib.import_module("api.agent_management.agents.geometry_agent")
)

# Alias the `api.dependencies` package so imports like `dependencies.use_llm` work
sys.modules.setdefault("dependencies", importlib.import_module("api.dependencies"))

# Provide a lightweight stub for `pymol` (used by render routes) when the real
# library is not available (e.g., during CI or unit testing where PyMOL isn't
# installed). The stub implements the minimal API surface needed by the code.

import types as _types

try:
    import pymol  # Attempt to import the real PyMOL if available
except ImportError:  # Fallback to a lightweight stub during CI / local dev
    if "pymol" not in sys.modules:
        _pymol = _types.ModuleType("pymol")

        class _CmdStub:
            def finish_launching(self, *args, **kwargs):
                pass

            def reinitialize(self):
                pass

            def get_view(self):
                return []

            def centerofmass(self):
                return []

            def get_extent(self):
                return ((), ())

            def png(self, *args, **kwargs):
                pass

            def save(self, *args, **kwargs):
                pass

            def do(self, *args, **kwargs):
                pass

        _pymol.cmd = _CmdStub()
        sys.modules["pymol"] = _pymol

# Provide a stub for the `redis` module (used by rate-limiting and caching)
# when Redis is not installed in the environment.

if "redis" not in sys.modules:
    class _RedisStubClient(dict):
        def pipeline(self):
            return self

        # Redis pipeline mimics
        def incr(self, key, amount=1):
            self[key] = int(self.get(key, 0)) + amount
            return self[key]

        def expire(self, key, ttl):
            return True

        def execute(self):
            # Return dummy values for ip_count, _, global_count, _
            return [1, None, 1, None]

        # Support same signature as redis-py's get(key) without default as
        # well as get(key, default).
        def get(self, key, default=None):  # type: ignore[override]
            value = super().get(key, default)
            return str(value).encode() if value is not None else None

        def set(self, key, value):
            super().__setitem__(key, value)

        # Pipeline context manager support
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            pass

    class _RedisStubModule:
        def Redis(self, *args, **kwargs):
            return _RedisStubClient()

    sys.modules["redis"] = _RedisStubModule()

# Alias routers package so `import routers` works
sys.modules.setdefault("routers", importlib.import_module("api.routers"))
