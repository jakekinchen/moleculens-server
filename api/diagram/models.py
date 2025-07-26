"""Pydantic models for molecular diagram generation."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


class Arrow(BaseModel):
    """Arrow connecting two points in a diagram."""

    start: tuple[float, float] = Field(description="Starting coordinates")
    end: tuple[float, float] = Field(description="Ending coordinates")
    text: str | None = Field(description="Optional text to display along the arrow")


class MoleculePlacement(BaseModel):
    """Placement information for a molecule in a diagram."""

    molecule: str
    x: float
    y: float
    width: float | None = None
    height: float | None = None
    label: str | None = None
    label_position: Literal["above", "below", "left", "right"] = "below"

    @property
    def position(self) -> tuple[float, float]:
        """Get the position as a tuple."""
        return (self.x, self.y)


class DiagramPlan(BaseModel):
    """Plan for rendering a molecular diagram."""

    plan: str = Field(description="Description of the diagram plan")
    molecule_list: list[MoleculePlacement] = Field(description="List of molecules and their positions")
    arrows: list[Arrow] = Field(description="List of arrows connecting molecules")
    canvas_width: int = Field(default=800, description="Width of the canvas in pixels")
    canvas_height: int = Field(default=600, description="Height of the canvas in pixels")
