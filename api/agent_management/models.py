"""
Pydantic models for Three.js structured output generation.
"""

from typing import List, Dict, Optional
from pydantic import BaseModel, Field

class Vector3(BaseModel):
    """Three.js Vector3 representation"""
    x: float = Field(default=0.0, description="X coordinate")
    y: float = Field(default=0.0, description="Y coordinate")
    z: float = Field(default=0.0, description="Z coordinate")

class Material(BaseModel):
    """Three.js Material properties"""
    type: str = Field(
        description="Material type (MeshBasic, MeshPhong, etc.)",
        pattern="^Mesh[A-Za-z]+Material$"
    )
    color: int = Field(description="Hex color value (e.g., 0xff0000 for red)")
    opacity: float = Field(default=1.0, ge=0.0, le=1.0, description="Material opacity")
    transparent: bool = Field(default=False, description="Whether material is transparent")
    shininess: Optional[float] = Field(default=None, description="Shininess for PhongMaterial")

class Geometry(BaseModel):
    """Three.js Geometry definition"""
    type: str = Field(
        description="Geometry type (Sphere, Cylinder, etc.)",
        pattern="^[A-Za-z]+Geometry$"
    )
    parameters: Dict[str, float] = Field(
        description="Geometry parameters (radius, height, etc.)"
    )

class Mesh(BaseModel):
    """Three.js Mesh combining geometry and material"""
    name: str = Field(description="Unique identifier for the mesh")
    geometry: Geometry
    material: Material
    position: Vector3 = Field(default_factory=Vector3)
    rotation: Vector3 = Field(default_factory=Vector3)
    scale: Vector3 = Field(default_factory=lambda: Vector3(x=1.0, y=1.0, z=1.0))

class ThreeGroup(BaseModel):
    """Three.js Group containing multiple meshes"""
    name: str = Field(description="Unique identifier for the group")
    position: Vector3 = Field(default_factory=Vector3)
    rotation: Vector3 = Field(default_factory=Vector3)
    scale: Vector3 = Field(default_factory=lambda: Vector3(x=1.0, y=1.0, z=1.0))
    children: List[Mesh] = Field(default_factory=list, description="List of meshes in the group") 