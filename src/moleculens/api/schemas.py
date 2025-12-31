"""Pydantic schemas for API request/response validation."""

from typing import Literal

from pydantic import BaseModel, Field, field_validator


class ComputeRequest(BaseModel):
    """Request schema for orbital computation."""

    sdf_content: str = Field(
        ...,
        description="3D SDF file content with atom coordinates",
        min_length=10,
    )
    method: str = Field(
        default="scf",
        description="Computational method (scf, b3lyp, hf, etc.)",
        pattern=r"^[a-zA-Z0-9_-]+$",
    )
    basis: str = Field(
        default="sto-3g",
        description="Basis set (sto-3g, 6-31g*, cc-pvdz, etc.)",
    )
    orbitals: list[str] = Field(
        default=["homo", "lumo"],
        description="Orbital types to compute",
    )
    grid_spacing: float = Field(
        default=0.25,
        description="Grid spacing in Angstrom",
        ge=0.1,
        le=1.0,
        alias="gridSpacing",
    )
    isovalue: float = Field(
        default=0.05,
        description="Isosurface value for mesh generation",
        ge=0.001,
        le=0.5,
    )
    inchi_key: str | None = Field(
        default=None,
        description="Optional InChIKey for cache deduplication",
        alias="inchiKey",
    )

    model_config = {"populate_by_name": True}

    @field_validator("orbitals")
    @classmethod
    def validate_orbitals(cls, v: list[str]) -> list[str]:
        allowed = {"homo", "lumo", "density", "homo-1", "lumo+1"}
        for orb in v:
            if orb.lower() not in allowed:
                raise ValueError(f"Invalid orbital type: {orb}. Allowed: {allowed}")
        return [orb.lower() for orb in v]

    @field_validator("basis")
    @classmethod
    def normalize_basis(cls, v: str) -> str:
        # Normalize common basis set names
        return v.lower().replace("*", "star")


class MeshDataResponse(BaseModel):
    """Response schema for a single mesh."""

    vertices: str = Field(..., description="Base64 gzipped Float32Array of vertices")
    normals: str = Field(..., description="Base64 gzipped Float32Array of normals")
    indices: str = Field(..., description="Base64 gzipped Uint32Array of triangle indices")
    vertexCount: int = Field(..., ge=0)
    triangleCount: int = Field(..., ge=0)


class OrbitalResponse(BaseModel):
    """Response schema for a single orbital."""

    positive: MeshDataResponse | None = None
    negative: MeshDataResponse | None = None
    energy_eV: float | None = Field(default=None, alias="energyEv")
    isovalue: float

    model_config = {"populate_by_name": True}


class ComputeResultMeta(BaseModel):
    """Metadata about the computation."""

    method: str
    basis: str
    grid_spacing_angstrom: float = Field(alias="gridSpacingAngstrom")
    psi4_num_threads: int = Field(alias="psi4NumThreads")
    compute_time_ms: float | None = Field(default=None, alias="computeTimeMs")

    model_config = {"populate_by_name": True}


class ComputeResult(BaseModel):
    """Full computation result."""

    orbitals: dict[str, OrbitalResponse] = Field(default_factory=dict)
    density: MeshDataResponse | None = None
    meta: ComputeResultMeta


class JobResponse(BaseModel):
    """Response schema for job status/result."""

    cached: bool = False
    job_id: str = Field(alias="jobId")
    status: Literal["pending", "queued", "running", "done", "error"]
    cache_key: str = Field(alias="cacheKey")
    result: ComputeResult | None = None
    error_message: str | None = Field(default=None, alias="errorMessage")
    created_at: str | None = Field(default=None, alias="createdAt")
    compute_time_ms: float | None = Field(default=None, alias="computeTimeMs")

    model_config = {"populate_by_name": True}


class HealthResponse(BaseModel):
    """Health check response."""

    status: Literal["healthy", "degraded", "unhealthy"]
    version: str
    database: bool
    cache_dir: bool = Field(alias="cacheDir")

    model_config = {"populate_by_name": True}
