"""Pydantic models for molecular diagram generation."""

from typing import List, Literal, Optional, Tuple

from pydantic import BaseModel, Field


class MoleculePosition(BaseModel):
    """Position and name of a molecule in a diagram."""

    name: str = Field(description="Name of the molecule")
    position: Tuple[float, float] = Field(
        description="X, Y coordinates of the molecule"
    )


class Arrow(BaseModel):
    """Arrow connecting two points in a diagram."""

    start: Tuple[float, float] = Field(description="Starting coordinates")
    end: Tuple[float, float] = Field(description="Ending coordinates")
    text: Optional[str] = Field(description="Optional text to display along the arrow")


class DiagramPlan(BaseModel):
    """Plan for rendering a molecular diagram."""

    plan: str = Field(description="Description of the diagram plan")
    molecule_list: List["MoleculePlacement"] = Field(
        description="List of molecules and their positions"
    )
    arrows: List[Arrow] = Field(description="List of arrows connecting molecules")
    canvas_width: int = Field(default=800, description="Width of the canvas in pixels")
    canvas_height: int = Field(
        default=600, description="Height of the canvas in pixels"
    )


class MoleculePlacement(BaseModel):
    """Placement information for a molecule in a diagram."""

    molecule: str
    x: float
    y: float
    width: Optional[float] = None
    height: Optional[float] = None
    label: Optional[str] = None
    label_position: Literal["above", "below", "left", "right"] = "below"

    @property
    def position(self) -> Tuple[float, float]:
        """Get the position as a tuple."""
        return (self.x, self.y)
