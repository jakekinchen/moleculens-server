"""SQLAlchemy database models."""

import enum
from datetime import datetime
from typing import Any

from sqlalchemy import JSON, DateTime, Enum, Float, Index, String, Text, func
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column


class Base(DeclarativeBase):
    """Base class for all models."""

    pass


class JobStatus(enum.Enum):
    """Job processing status."""

    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"


class Job(Base):
    """Job queue entry for orbital computation requests."""

    __tablename__ = "jobs"

    # Primary key
    id: Mapped[str] = mapped_column(String(64), primary_key=True)

    # Cache key for deduplication and artifact storage
    cache_key: Mapped[str] = mapped_column(String(64), nullable=False, index=True)

    # Job status
    status: Mapped[JobStatus] = mapped_column(
        Enum(JobStatus), default=JobStatus.PENDING, nullable=False
    )

    # Request parameters (stored as JSON)
    request_params: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)

    # Result metadata (filled on completion)
    result_meta: Mapped[dict[str, Any] | None] = mapped_column(JSON, nullable=True)

    # Error message (if status == ERROR)
    error_message: Mapped[str | None] = mapped_column(Text, nullable=True)

    # Timing
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )
    started_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    compute_time_ms: Mapped[float | None] = mapped_column(Float, nullable=True)

    # Worker info
    worker_id: Mapped[str | None] = mapped_column(String(64), nullable=True)

    # Retry tracking
    attempt_count: Mapped[int] = mapped_column(default=0)

    # Indexes for efficient queries
    __table_args__ = (
        Index("ix_jobs_status_created", "status", "created_at"),
        Index("ix_jobs_cache_key_status", "cache_key", "status"),
    )

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for API response."""
        return {
            "jobId": self.id,
            "cacheKey": self.cache_key,
            "status": self.status.value,
            "createdAt": self.created_at.isoformat() if self.created_at else None,
            "startedAt": self.started_at.isoformat() if self.started_at else None,
            "completedAt": self.completed_at.isoformat() if self.completed_at else None,
            "computeTimeMs": self.compute_time_ms,
            "attemptCount": self.attempt_count,
            "errorMessage": self.error_message,
        }


class CacheEntry(Base):
    """Cache metadata for computed artifacts."""

    __tablename__ = "cache_entries"

    # Cache key is the primary key
    cache_key: Mapped[str] = mapped_column(String(64), primary_key=True)

    # Input parameters hash (for verification)
    params_hash: Mapped[str] = mapped_column(String(64), nullable=False)

    # Result metadata
    result_meta: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)

    # Artifact paths (relative to cache dir)
    artifact_paths: Mapped[list[str]] = mapped_column(JSON, nullable=False)

    # Timing
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), server_default=func.now(), nullable=False
    )

    # Hit count for cache analytics
    hit_count: Mapped[int] = mapped_column(default=0)
    last_accessed_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True), nullable=True
    )
