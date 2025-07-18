"""
Pydantic models for Three.js structured output generation.
"""

from typing import List, Dict, Optional, Any, Union, Type, TypeVar, Callable, Tuple, Literal
from pydantic import BaseModel, Field


# Global Pydantic model config to avoid the "model_name" warning
class BaseModelWithConfig(BaseModel):
    """Base model with config that disables protected namespaces"""

    model_config = {
        "protected_namespaces": ()  # Remove protection for the 'model_' namespace
    }


# Define a type variable for models
T = TypeVar("T")


class ModelRegistry:
    """Registry for managing model classes and factory functions."""
    
    _registry: Dict[str, Tuple[Type[BaseModel], Callable[..., Any]]] = {}
    
    @classmethod
    def register(cls, model_name: str, model_cls: Type[BaseModel], factory_func: Callable[..., Any]):
        """
        Register a model class and its factory function.
        
        Args:
            model_name: Name to register the model under
            model_cls: The model class
            factory_func: Factory function to create model instances
        """
        cls._registry[model_name] = (model_cls, factory_func)
    
    @classmethod
    def get_model(cls, model_name: str) -> Type[BaseModel]:
        """
        Get the registered model class.
        
        Args:
            model_name: Name of the registered model
            
        Returns:
            The model class
            
        Raises:
            ValueError: If the model is not registered
        """
        if model_name not in cls._registry:
            raise ValueError(f"Model '{model_name}' not registered.")
        return cls._registry[model_name][0]
    
    @classmethod
    def create_instance(cls, model_name: str, *args, **kwargs):
        """
        Create an instance of a registered model.
        
        Args:
            model_name: Name of the registered model
            *args, **kwargs: Arguments to pass to the factory function
            
        Returns:
            An instance of the model
            
        Raises:
            ValueError: If the model is not registered
        """
        if model_name not in cls._registry:
            raise ValueError(f"Model '{model_name}' not registered.")
        _, factory_func = cls._registry[model_name]
        return factory_func(*args, **kwargs)
    
    @classmethod
    def list_models(cls) -> List[str]:
        """Get a list of all registered model names."""
        return list(cls._registry.keys())


class Vector3(BaseModel):
    """Three.js Vector3 representation"""

    x: float = Field(default=0.0, description="X coordinate")
    y: float = Field(default=0.0, description="Y coordinate")
    z: float = Field(default=0.0, description="Z coordinate")


class Material(BaseModel):
    """Three.js Material properties"""

    type: str = Field(
        description="Material type (MeshBasic, MeshPhong, etc.)",
        pattern="^Mesh[A-Za-z]+Material$",
    )
    color: int = Field(description="Hex color value (e.g., 0xff0000 for red)")
    opacity: float = Field(default=1.0, ge=0.0, le=1.0, description="Material opacity")
    transparent: bool = Field(
        default=False, description="Whether material is transparent"
    )
    shininess: Optional[float] = Field(
        default=None, description="Shininess for PhongMaterial"
    )


class Geometry(BaseModel):
    """Three.js Geometry definition"""

    type: str = Field(
        description="Geometry type (Sphere, Cylinder, etc.)",
        pattern="^[A-Za-z]+Geometry$",
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
    children: List[Mesh] = Field(
        default_factory=list, description="List of meshes in the group"
    )


class BooleanResponse(BaseModel):
    """Model for boolean validation responses"""

    is_true: bool = Field(description="Whether the statement is true")


class MolecularStructure(BaseModel):
    """Model for molecular validation responses"""

    molecular_structure: str = Field(description="Molecular structure in SMILES format")


class ScriptTimePoint(BaseModel):
    """Model for a single time point in an animation script"""

    timecode: str = Field(description="Timecode in MM:SS format (e.g., '00:15')")
    atoms: List[str] = Field(description="List of atoms to highlight in the scene")
    caption: str = Field(description="Brief caption text to be displayed at this time")

    model_config = {
        "protected_namespaces": (),  # Remove protection for model_ namespace
        "json_schema_extra": {
            "examples": [
                {
                    "timecode": "00:15",
                    "atoms": ["C1", "C2", "H1"],
                    "caption": "Carbon-carbon bonds form the backbone of organic molecules",
                }
            ]
        },
    }


class SceneScript(BaseModel):
    """Model for a complete scene script (not to be confused with animation)"""

    title: str = Field(description="Title of the scene")
    content: List[ScriptTimePoint] = Field(
        description="Sequence of time points in the scene"
    )


class SceneObject(BaseModel):
    """Model for a discrete 3D object needed in the scene"""

    name: str = Field(description="Unique identifier for the object")
    description: Union[str, List[str]] = Field(
        description="Detailed description of what this object represents"
    )
    properties: Dict[str, Any] = Field(
        description="Key properties of the object (size, color, etc.)"
    )
    appears_at: str = Field(description="Timecode when this object first appears")
    relationships: List[str] = Field(
        default_factory=list, description="Relationships to other objects"
    )


class OrchestrationPlan(BaseModel):
    """Model for the complete orchestration plan of objects needed for the scene"""

    scene_title: str = Field(description="Title of the scene")
    objects: List[SceneObject] = Field(
        description="List of discrete objects needed for the scene"
    )
    scene_overview: Optional[str] = Field(
        default=None, description="High-level overview of the overall scene"
    )


class AnimationKeyframe(BaseModel):
    """Model for an animation keyframe"""

    timecode: str = Field(
        description="Timecode in MM:SS format when this keyframe occurs"
    )
    actions: List[str] = Field(
        description="Animation actions to perform at this keyframe"
    )


class AnimationCode(BaseModelWithConfig):
    """Model for the complete animation code output"""

    code: str = Field(
        description="Animation code for the scene (content of the animate function)"
    )
    keyframes: List[AnimationKeyframe] = Field(
        description="List of keyframes in the animation"
    )


class FinalScenePackage(BaseModelWithConfig):
    """Model for the complete packaged Three.js scene"""

    html: str = Field(description="Complete HTML with embedded JavaScript")
    js: str = Field(description="Complete JavaScript code for the scene (standalone)")
    minimal_js: str = Field(
        description="Minimal JavaScript code without boilerplate for embedding"
    )
    title: str = Field(description="Title of the visualization")
    timecode_markers: List[str] = Field(
        description="List of timecode markers in the animation"
    )
    total_elements: int = Field(description="Total number of 3D elements in the scene")


class PubChemCompound(BaseModel):
    """Model representing a PubChem compound with its properties"""

    name: str = Field(description="The name used to query the compound")
    cid: int = Field(description="PubChem Compound ID")
    molecular_formula: str = Field(description="Molecular formula of the compound")
    molecular_weight: float = Field(description="Molecular weight of the compound")
    iupac_name: Optional[str] = Field(
        description="IUPAC name of the compound", default=None
    )
    sdf: Optional[str] = Field(description="SDF format structure data", default=None)
    canonical_smiles: Optional[str] = Field(
        description="Canonical SMILES representation of the compound", default=None
    )
    isomeric_smiles: Optional[str] = Field(
        description="Isomeric SMILES representation of the compound", default=None
    )
    elements: Optional[List[str]] = Field(
        description="List of elements in the compound", default=None
    )
    atoms: Optional[List[int]] = Field(
        description="List of atom indices in the compound", default=None
    )
    bonds: Optional[List[Tuple[int, int, int]]] = Field(
        description="List of bonds in the compound", default=None
    )
    charge: Optional[int] = Field(description="Charge of the compound", default=None)
    synonyms: Optional[List[str]] = Field(
        description="List of synonyms for the compound", default=None
    )


class PubChemSearchResult(BaseModel):
    """Model representing the results of a PubChem search"""

    query: str = Field(description="Original user query")
    interpreted_query: str = Field(description="Query interpreted by LLM")
    results: List[PubChemCompound] = Field(
        description="List of matching compounds", default_factory=list
    )


class Box2D(BaseModel):
    """2-D bounding box for molecule placement."""

    x: float
    y: float
    width: float
    height: float


class MoleculeBoxRequest(BaseModel):
    """Query for a molecule along with its layout box."""

    query: str
    box: Box2D


class MoleculeLayoutRequest(BaseModel):
    """Request containing multiple molecules and their boxes."""

    molecules: List[MoleculeBoxRequest]


class MoleculePosition(BaseModel):
    """Position and name of a molecule in a diagram"""
    name: str = Field(description="Name of the molecule")
    position: Tuple[float, float] = Field(description="X, Y coordinates of the molecule")


class Arrow(BaseModel):
    """Arrow connecting two points in a diagram"""
    start: Tuple[float, float] = Field(description="Starting coordinates")
    end: Tuple[float, float] = Field(description="Ending coordinates")
    text: Optional[str] = Field(description="Optional text to display along the arrow")


class DiagramPlan(BaseModel):
    """Plan for rendering a molecular diagram"""
    plan: str = Field(description="Description of the diagram plan")
    molecule_list: List[MoleculePosition] = Field(description="List of molecules and their positions")
    arrows: List[Arrow] = Field(description="List of arrows connecting molecules")
    canvas_width: int = Field(default=800, description="Width of the canvas in pixels")
    canvas_height: int = Field(default=600, description="Height of the canvas in pixels")


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
