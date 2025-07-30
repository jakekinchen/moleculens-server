from __future__ import annotations

from typing import Any, List, Literal, Optional

from pydantic import BaseModel, Field


class RenderingOptions(BaseModel):
    """Advanced rendering options for PyMOL output."""

    transparent_background: bool = False
    ray_trace: bool = True
    resolution: List[int] = Field(default_factory=lambda: [1920, 1080])
    dpi: int = 300
    ray_trace_mode: Literal["default", "cartoon_outline", "bw", "poster"] = "default"
    antialias: bool = True
    ray_shadow: bool = True
    depth_cue: bool = True
    ray_trace_fog: bool = False
    background_color: str = "white"
    field_of_view: float = 20.0
    orthoscopic: bool = False


class SceneSpec(BaseModel):
    """Validated payload for every PyMOL request."""

    op: Literal[
        "overview",
        "binding_site",
        "mutation",
        "mutation_focus",
        "raw",
        # NEW: Advanced rendering operations
        "transparent_molecule",
        "publication_quality",
        "annotated_molecule",
        "transparent_binding_site",
    ]
    structure_id: str
    selection: Optional[str] = Field(
        default=None,
        description="PyMOL atom-selection targeting a subset of atoms.",
    )
    opts: dict[str, Any] = Field(
        default_factory=dict,
        description="Operation-specific keyword arguments.",
    )
    raw_cmds: Optional[List[str]] = Field(
        default=None,
        description="Explicit PyMOL commands when op == 'raw'.",
    )
    # NEW: Rendering specifications
    rendering: RenderingOptions = Field(
        default_factory=RenderingOptions,
        description="Advanced rendering options for high-quality output.",
    )
