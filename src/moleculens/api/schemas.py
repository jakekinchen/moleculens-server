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


# =============================================================================
# Electrostatics Schemas
# =============================================================================


class SurfaceParams(BaseModel):
    """Surface generation parameters."""

    grid_spacing: float = Field(default=0.25, ge=0.1, le=1.0, alias="gridSpacing")
    probe_radius: float = Field(default=1.4, ge=0.0, le=5.0, alias="probeRadius")
    padding: float = Field(default=3.0, ge=1.0, le=10.0)

    model_config = {"populate_by_name": True}


class PotentialParams(BaseModel):
    """Potential calculation parameters."""

    units: str = Field(default="kcal/mol/e")
    softening_epsilon: float = Field(default=0.3, ge=0.1, le=2.0, alias="softeningEpsilon")
    clamp_percentiles: list[float] = Field(default=[5.0, 95.0], alias="clampPercentiles")

    model_config = {"populate_by_name": True}

    @field_validator("clamp_percentiles")
    @classmethod
    def validate_percentiles(cls, v: list[float]) -> list[float]:
        if len(v) != 2:
            raise ValueError("clamp_percentiles must have exactly 2 values")
        if not (0 <= v[0] < v[1] <= 100):
            raise ValueError("clamp_percentiles must be [low, high] with 0 <= low < high <= 100")
        return v


class ElectrostaticsRequest(BaseModel):
    """Request schema for electrostatics computation."""

    sdf_content: str = Field(
        ...,
        description="3D SDF file content with atom coordinates",
        min_length=10,
        alias="sdfContent",
    )
    smiles: str | None = Field(default=None, description="Optional SMILES string")
    inchi_key: str | None = Field(
        default=None,
        description="Optional InChIKey for cache deduplication",
        alias="inchiKey",
    )
    client_cache_key: str | None = Field(
        default=None,
        description="Optional client-provided cache key",
        alias="clientCacheKey",
    )
    charge: int = Field(default=0, ge=-10, le=10, description="Total molecular charge")
    multiplicity: int = Field(default=1, ge=1, le=10, description="Spin multiplicity")
    method: str = Field(
        default="xtb-gfn2",
        description="Charge method: 'xtb-gfn2' or 'gasteiger'",
        pattern=r"^(xtb-gfn2|xtb|gasteiger)$",
    )
    surface: SurfaceParams = Field(default_factory=SurfaceParams)
    potential: PotentialParams = Field(default_factory=PotentialParams)

    model_config = {"populate_by_name": True}


class ESPMeshResponse(BaseModel):
    """Response schema for ESP surface mesh."""

    vertices: str = Field(..., description="Base64 gzipped Float32Array of vertices [x0,y0,z0,...]")
    normals: str = Field(..., description="Base64 gzipped Float32Array of normals")
    indices: str = Field(..., description="Base64 gzipped Uint32Array of triangle indices")
    potential: str = Field(..., description="Base64 gzipped Float32Array of per-vertex potential")
    vertexCount: int = Field(..., ge=0)
    triangleCount: int = Field(..., ge=0)


class PotentialRangeResponse(BaseModel):
    """Response schema for potential range statistics."""

    min: float
    max: float
    p05: float | None = None
    p95: float | None = None


class UnitsResponse(BaseModel):
    """Response schema for units specification."""

    coords: str = "angstrom"
    potential: str = "kcal/mol/e"


class SurfaceResponse(BaseModel):
    """Response schema for ESP surface."""

    mesh: ESPMeshResponse
    range: PotentialRangeResponse
    units: UnitsResponse


class DipoleResponse(BaseModel):
    """Response schema for dipole moment."""

    origin_angstrom: list[float] = Field(alias="originAngstrom")
    vector_debye: list[float] = Field(alias="vectorDebye")
    magnitude_debye: float = Field(alias="magnitudeDebye")
    convention: Literal["physics", "chemistry"]

    model_config = {"populate_by_name": True}


class ElectrostaticsTimings(BaseModel):
    """Timing information for electrostatics computation."""

    total: float
    surface: float | None = None
    charges: float | None = None
    potential: float | None = None


class ElectrostaticsMeta(BaseModel):
    """Metadata for electrostatics computation."""

    charge_model: str = Field(alias="chargeModel")
    method: str
    timings_ms: ElectrostaticsTimings = Field(alias="timingsMs")
    sdf_hash: str = Field(alias="sdfHash")
    molecular_charge: int = Field(default=0, alias="molecularCharge")
    multiplicity: int = Field(default=1)

    model_config = {"populate_by_name": True}


class ElectrostaticsResult(BaseModel):
    """Full electrostatics computation result."""

    surface: SurfaceResponse
    dipole: DipoleResponse
    meta: ElectrostaticsMeta


class ElectrostaticsJobResponse(BaseModel):
    """Response schema for electrostatics job status/result."""

    cached: bool = False
    job_id: str = Field(alias="jobId")
    status: Literal["pending", "queued", "running", "done", "error"]
    cache_key: str = Field(alias="cacheKey")
    progress: float | None = Field(default=None, ge=0, le=100)
    step: str | None = None
    result: ElectrostaticsResult | None = None
    error_message: str | None = Field(default=None, alias="errorMessage")
    created_at: str | None = Field(default=None, alias="createdAt")
    compute_time_ms: float | None = Field(default=None, alias="computeTimeMs")

    model_config = {"populate_by_name": True}


# =============================================================================
# Conformer Schemas
# =============================================================================


class ConformerParams(BaseModel):
    """Parameters for conformer generation."""

    method: str = Field(
        default="etkdg_v3",
        description="Embedding method (etkdg_v3, etkdg_v2, etkdg)",
        pattern=r"^[a-zA-Z0-9_-]+$",
    )
    opt: str = Field(
        default="uff",
        description="Optimization method (uff, mmff, none)",
        pattern=r"^[a-zA-Z0-9_-]+$",
    )
    max_attempts: int = Field(default=10, ge=1, le=100, alias="maxAttempts")
    max_opt_iters: int = Field(default=200, ge=0, le=2000, alias="maxOptIters")
    add_hs: bool = Field(default=False, alias="addHs")

    model_config = {"populate_by_name": True}

    @field_validator("method")
    @classmethod
    def normalize_method(cls, v: str) -> str:
        return v.lower()

    @field_validator("opt")
    @classmethod
    def normalize_opt(cls, v: str) -> str:
        v = v.lower()
        if v in ("off", "false", "no"):
            return "none"
        return v


class ConformerRequest(BaseModel):
    """Request schema for conformer computation."""

    molblock2d: str = Field(
        ...,
        description="2D molblock/SDF content with atom ordering",
        min_length=10,
    )
    params: ConformerParams = Field(default_factory=ConformerParams)
    wait_ms: int = Field(
        default=0,
        ge=0,
        le=10000,
        description="Max time to wait for job completion",
        alias="waitMs",
    )
    geom_version: str = Field(
        ...,
        min_length=1,
        description="Geometry cache version",
        alias="geomVersion",
    )

    model_config = {"populate_by_name": True}


class ConformerTimingMs(BaseModel):
    """Timing information for conformer computation."""

    embed: float | None = None
    opt: float | None = None
    total: float | None = None


class ConformerMeta(BaseModel):
    """Metadata for conformer computation."""

    method: str
    opt: str
    max_attempts: int = Field(alias="maxAttempts")
    max_opt_iters: int = Field(alias="maxOptIters")
    add_hs: bool = Field(alias="addHs")
    geom_version: str = Field(alias="geomVersion")
    quality: str
    has_metal: bool = Field(alias="hasMetal")
    embed_status: str = Field(alias="embedStatus")
    opt_status: str = Field(alias="optStatus")

    model_config = {"populate_by_name": True}


class ConformerJobResponse(BaseModel):
    """Response schema for conformer job status/result."""

    status: Literal["pending", "done", "failed"]
    cache_key: str = Field(alias="cacheKey")
    job_id: str | None = Field(default=None, alias="jobId")
    sdf3d: str | None = None
    meta: ConformerMeta | None = None
    timing_ms: ConformerTimingMs | None = Field(default=None, alias="timingMs")
    error_message: str | None = Field(default=None, alias="errorMessage")

    model_config = {"populate_by_name": True}
