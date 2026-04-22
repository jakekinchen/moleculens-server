"""In-process cache event counters for operational visibility."""

from collections import Counter, defaultdict
from datetime import UTC, datetime
from threading import Lock
from typing import Any


class CacheMetrics:
    """Thread-safe cache metrics registry for the current process."""

    def __init__(self) -> None:
        self._lock = Lock()
        self._reset_unlocked()

    def _reset_unlocked(self) -> None:
        self._started_at = datetime.now(UTC)
        self._events: Counter[str] = Counter()
        self._job_types: dict[str, Counter[str]] = defaultdict(Counter)
        self._identity_sources: dict[str, Counter[str]] = defaultdict(Counter)

    def reset(self) -> None:
        with self._lock:
            self._reset_unlocked()

    def record(
        self,
        *,
        event: str,
        job_type: str | None = None,
        cache_identity_source: str | None = None,
    ) -> None:
        with self._lock:
            self._events[event] += 1
            if job_type is not None:
                self._job_types[job_type][event] += 1
            if cache_identity_source is not None:
                self._identity_sources[cache_identity_source][event] += 1

    def snapshot(self) -> dict[str, Any]:
        with self._lock:
            recorded_at = datetime.now(UTC)
            uptime_seconds = max(
                0.0,
                (recorded_at - self._started_at).total_seconds(),
            )
            return {
                "started_at": self._started_at.isoformat(),
                "recorded_at": recorded_at.isoformat(),
                "uptime_seconds": round(uptime_seconds, 3),
                "total_events": sum(self._events.values()),
                "events": dict(sorted(self._events.items())),
                "job_types": {
                    job_type: dict(sorted(counter.items()))
                    for job_type, counter in sorted(self._job_types.items())
                },
                "identity_sources": {
                    source: dict(sorted(counter.items()))
                    for source, counter in sorted(self._identity_sources.items())
                },
            }


cache_metrics = CacheMetrics()

