"""
Pydantic models for Three.js structured output generation.
"""

from typing import List, Dict, Optional, Any, Union, Type, TypeVar, Callable
from pydantic import BaseModel, Field

# Global Pydantic model config to avoid the "model_name" warning
class BaseModelWithConfig(BaseModel):
    """Base model with config that disables protected namespaces"""
    model_config = {
        "protected_namespaces": ()  # Remove protection for the 'model_' namespace
    }

# Define a type variable for models
T = TypeVar('T')

class ModelRegistry:
    """
    Central registry for managing model classes and their instantiation.
    Implements the Registry design pattern to decouple model references.
    """
    _registry = {}
    
    @classmethod
    def register(cls, name: str, model_cls, factory: Optional[Callable] = None):
        """
        Register a model class with an optional factory function.
        
        Args:
            name: Unique identifier for the model
            model_cls: The model class to register
            factory: Optional factory function to instantiate the model
        """
        if factory is None:
            # Create a default factory that accepts kwargs
            factory = lambda *args, **kwargs: model_cls(**kwargs)
            
        cls._registry[name] = {
            "class": model_cls,
            "factory": factory
        }
        return model_cls  # Return for decorator use
    
    @classmethod
    def get_model(cls, name: str):
        """
        Get a model class by name.
        
        Args:
            name: Name of the registered model
            
        Returns:
            The model class
            
        Raises:
            ValueError: If the model is not registered
        """
        model_data = cls._registry.get(name)
        if model_data is None:
            raise ValueError(f"Model '{name}' not registered.")
        return model_data["class"]
    
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
        model_data = cls._registry.get(model_name)
        if model_data is None:
            raise ValueError(f"Model '{model_name}' not registered.")
        
        factory = model_data["factory"]
        return factory(*args, **kwargs)
    
    @classmethod
    def list_models(cls):
        """List all registered model names"""
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

class BooleanResponse(BaseModel):
    """Model for boolean validation responses"""
    is_true: bool = Field(description="Whether the statement is true")

class MolecularStructure(BaseModel):
    """Model for molecular validation responses"""
    molecular_structure: str = Field(description="Molecular structure in SMILES format")

class ScriptTimePoint(BaseModel):
    """Model for a single time point in an animation script"""
    timecode: str = Field(description="Timecode in MM:SS format (e.g., '00:15')")
    description: str = Field(description="Detailed description of what is occurring in the scene at this time")
    caption: str = Field(description="Brief caption text to be displayed at this time")

class SceneScript(BaseModel):
    """Model for a complete scene script (not to be confused with animation)"""
    title: str = Field(description="Title of the scene")
    content: List[ScriptTimePoint] = Field(description="Sequence of time points in the scene")

class SceneObject(BaseModel):
    """Model for a discrete 3D object needed in the scene"""
    name: str = Field(description="Unique identifier for the object")
    description: Union[str, List[str]] = Field(description="Detailed description of what this object represents")
    properties: Dict[str, Any] = Field(description="Key properties of the object (size, color, etc.)")
    appears_at: str = Field(description="Timecode when this object first appears")
    relationships: List[str] = Field(default_factory=list, description="Relationships to other objects")

class OrchestrationPlan(BaseModel):
    """Model for the complete orchestration plan of objects needed for the scene"""
    scene_title: str = Field(description="Title of the scene")
    objects: List[SceneObject] = Field(description="List of discrete objects needed for the scene")
    scene_overview: Optional[str] = Field(default=None, description="High-level overview of the overall scene")

class AnimationKeyframe(BaseModel):
    """Model for an animation keyframe"""
    timecode: str = Field(description="Timecode in MM:SS format when this keyframe occurs")
    actions: List[str] = Field(description="Animation actions to perform at this keyframe")
    
class AnimationCode(BaseModelWithConfig):
    """Model for the complete animation code output"""
    code: str = Field(description="Animation code for the scene (content of the animate function)")
    keyframes: List[AnimationKeyframe] = Field(description="List of keyframes in the animation")

class FinalScenePackage(BaseModelWithConfig):
    """Model for the complete packaged Three.js scene"""
    html: str = Field(description="Complete HTML with embedded JavaScript")
    js: str = Field(description="Complete JavaScript code for the scene (standalone)")
    minimal_js: str = Field(description="Minimal JavaScript code without boilerplate for embedding")
    title: str = Field(description="Title of the visualization")
    timecode_markers: List[str] = Field(description="List of timecode markers in the animation")
    total_elements: int = Field(description="Total number of 3D elements in the scene")

class PubChemCompound(BaseModel):
    """Model representing a PubChem compound with its properties"""
    name: str = Field(description="The name used to query the compound")
    cid: int = Field(description="PubChem Compound ID")
    molecular_formula: str = Field(description="Molecular formula of the compound")
    molecular_weight: float = Field(description="Molecular weight of the compound")
    iupac_name: Optional[str] = Field(description="IUPAC name of the compound", default=None)
    sdf: Optional[str] = Field(description="SDF format structure data", default=None)

class PubChemSearchResult(BaseModel):
    """Model representing the results of a PubChem search"""
    query: str = Field(description="Original user query")
    interpreted_query: str = Field(description="Query interpreted by LLM")
    results: List[PubChemCompound] = Field(description="List of matching compounds", default_factory=list)