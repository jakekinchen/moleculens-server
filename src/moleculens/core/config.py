"""Application configuration using pydantic-settings."""

from functools import lru_cache
from pathlib import Path
from typing import Literal

from pydantic import field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        case_sensitive=False,
        extra="ignore",
    )

    # Database
    database_url: str = "postgresql://moleculens:moleculens@localhost:5432/moleculens"

    # Filesystem paths
    molecule_cache_dir: Path = Path("/var/molecule-cache")
    psi_scratch: Path = Path("/var/psi4_scratch")
    log_dir: Path = Path("/var/log/molecule")

    # Psi4 configuration
    psi4_num_threads: int = 3
    psi4_memory_mb: int = 4096

    # Worker configuration
    worker_poll_seconds: float = 1.0
    worker_job_timeout_seconds: int = 600  # 10 minutes max per job

    # API configuration
    allowed_origins: str = "*"
    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = "INFO"

    # Job queue configuration
    max_retries: int = 3
    job_stale_seconds: int = 3600  # Consider job stale after 1 hour

    @field_validator("allowed_origins", mode="before")
    @classmethod
    def parse_origins(cls, v: str | list[str]) -> str:
        if isinstance(v, list):
            return ",".join(v)
        return v

    @property
    def cors_origins(self) -> list[str]:
        """Parse ALLOWED_ORIGINS into a list."""
        if self.allowed_origins == "*":
            return ["*"]
        return [origin.strip() for origin in self.allowed_origins.split(",") if origin.strip()]

    def ensure_directories(self) -> None:
        """Create required directories if they don't exist."""
        self.molecule_cache_dir.mkdir(parents=True, exist_ok=True)
        self.psi_scratch.mkdir(parents=True, exist_ok=True)
        self.log_dir.mkdir(parents=True, exist_ok=True)


@lru_cache
def get_settings() -> Settings:
    """Get cached settings instance."""
    return Settings()


settings = get_settings()
