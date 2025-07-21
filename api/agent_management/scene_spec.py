from __future__ import annotations

from typing import Literal, Optional, List, Dict, Any
from pydantic import BaseModel, Field


class SceneSpec(BaseModel):
    """Validated payload for every PyMOL request."""

    op: Literal["overview", "binding_site", "mutation", "raw"]
    structure_id: str
    selection: Optional[str] = Field(
        default=None,
        description="PyMOL atom-selection targeting a subset of atoms.",
    )
    opts: Dict[str, Any] = Field(
        default_factory=dict,
        description="Operation-specific keyword arguments.",
    )
    raw_cmds: Optional[List[str]] = Field(
        default=None,
        description="Explicit PyMOL commands when op == 'raw'.",
    )
