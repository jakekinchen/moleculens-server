"""Routes for handling prompts and generating visualizations."""

import json
import logging
import os
import traceback
from datetime import datetime
from typing import Any, Dict, List, Literal, Optional, Tuple
from uuid import uuid4

from pydantic import BaseModel, Field

from api.agent_management.agent_factory import AgentFactory
from api.agent_management.agents.pubchem_agent import PubChemAgent, _sdf_to_pdb_block
from api.agent_management.model_config import ModelCategory
from api.agent_management.models import (
    BaseModelWithConfig,
    OrchestrationPlan,
    SceneScript,
)
from api.agent_management.providers.openai_provider import generate_structured
from api.agent_management.scene_packager import ScenePackager
from api.dependencies.use_llm import use_llm
from api.diagram.diagram_renderer import render_diagram
from api.diagram.models import DiagramPlan
from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException

# Diagram Pydantic models are imported from api.agent_management.models to avoid duplication

# Instantiate the router AFTER all Pydantic models and necessary imports are defined
router = APIRouter(
    prefix="/prompt",
    tags=["Prompt"],
    responses={404: {"description": "Not found"}},
)

# Set up logger
logger = logging.getLogger(__name__)


# Other Pydantic models used by other endpoints can be defined here or remain in their original positions
class PromptRequest(BaseModel):
    """Request model for prompts, optionally with preferred model type."""

    prompt: str
    model: Optional[str] = None
    preferred_model_category: Optional[ModelCategory] = None


class PromptWorking(BaseModel):
    prompt: str
    working: bool
    script: Optional[SceneScript] = None
    orchestration_plan: Optional[OrchestrationPlan] = None


class CompletePipelineResponse(BaseModelWithConfig):
    """Response model for the complete prompt-to-visualization pipeline."""

    html: str
    js: str
    minimal_js: str
    title: str
    timecode_markers: List[str]
    total_elements: int


class DiagramPromptRequest(BaseModel):
    prompt: str
    canvas_width: int = 960
    canvas_height: int = 640
    model: Optional[str] = None
    preferred_model_category: Optional[ModelCategory] = None


class DiagramResponse(BaseModel):
    diagram_image: str  # base-64 PNG or SVG string
    diagram_plan: DiagramPlan
    status: Literal["completed", "failed", "processing"] = "completed"
    job_id: Optional[str] = None
    error: Optional[str] = None


class GeometryResponse(BaseModel):
    result: str
    is_molecular: bool = True
    validation_message: Optional[str] = None


class GeometryRequest(BaseModel):
    """Request model for geometry generation, optionally with preferred model
    type."""

    prompt: str
    model: Optional[str] = None  # Model override for geometry agent
    preferred_model_category: Optional[ModelCategory] = None


# Module-level structured models
class ValidationResponse(BaseModel):
    is_molecular: bool


class MoleculePackage(BaseModel):
    molecule_names: List[str]


class Molecule2DData(BaseModel):
    molecule: str
    x: float
    y: float
    width: float
    height: float
    label: Optional[str] = None
    label_position: Optional[str] = None


class Molecule2DDataList(BaseModel):
    __root__: List[Molecule2DData]


class VisualizationData(BaseModel):
    html: str
    pdb_data: str
    title: str
    timecode_markers: List[str]
    total_elements: int


class JobResponse(BaseModel):
    job_id: str
    status: str
    message: str
    result: Optional[str] = None
    progress: Optional[float] = None
    visualization: Optional[VisualizationData] = None
    error: Optional[str] = None


@router.post("/generate-from-pubchem/", response_model=dict)
async def generate_from_pubchem(request: PromptRequest):
    """Generate a molecule name or names from a prompt.

    Returns:
        - molecule_names: List[str]
    """
    try:
        logger.info(f"Processing PubChem request: {request.prompt}")

        # Use module-level ValidationResponse

        validation_result = generate_structured(
            user_prompt=request.prompt,
            response_model=ValidationResponse,
            system_prompt="You are a helpful assistant that validates if a prompt is scientific in nature.",
        )

        if not validation_result.is_molecular:
            raise HTTPException(status_code=400, detail="Non-molecular prompt rejected")

        # Use module-level MoleculePackage

        logger.info(f"Generating molecule package for: {request.prompt}")
        result = generate_structured(
            user_prompt=request.prompt,
            response_model=MoleculePackage,
            system_prompt="You are a helpful assistant that generates a molecule name or names for a given prompt.",
        )
        logger.info(f"Successfully generated molecule package: {result.molecule_names}")

        return {
            "molecule_names": result.molecule_names,
        }
    except Exception as e:
        error_message = str(e)
        logger.error(f"Error generating from PubChem: {error_message}")
        logger.error(traceback.format_exc())
        raise HTTPException(status_code=500, detail=error_message)


@router.post("/", response_model=JobResponse)
async def submit_prompt(request: PromptRequest, background_tasks: BackgroundTasks):
    """End-to-end endpoint to process a scientific prompt into a complete 3D
    visualization. This endpoint now launches a background task and returns
    immediately with a job ID.

    The client should poll the /prompt/process/{job_id} endpoint to check the status
    and get the final result.

    Pipeline steps:
    1. Validates that the prompt is scientific
    2. Generates an animation script with timecodes, descriptions, and captions
    3. Creates an orchestration plan with discrete objects needed for the animation
    4. Generates Three.js geometry for each object in the plan
    5. Creates animation code based on the script and objects
    6. Packages everything into a complete scene
    """
    try:
        # Use module-level ValidationResponse

        validation_result = generate_structured(
            user_prompt=request.prompt,
            response_model=ValidationResponse,
            system_prompt="You are a helpful assistant that validates if a prompt is scientific in nature.",
        )

        # If not scientific, reject the prompt immediately
        if not validation_result.is_molecular:
            return {
                "job_id": "rejected",
                "status": "failed",
                "message": "Non-molecular prompt rejected",
                "error": "The prompt does not contain molecular content",
            }

        job_id = str(uuid4())

        # Schedule background pipeline processing (stub)
        async def _process_pipeline(job_id: str, prompt: str):
            logger.info(f"Background processing for job {job_id}: {prompt}")

        background_tasks.add_task(_process_pipeline, job_id, request.prompt)
        return {
            "job_id": job_id,
            "status": "processing",
            "message": "Processing started in the background. Poll /prompt/process/{job_id} for updates.",
        }
    except Exception as e:
        print(f"Unhandled exception in submit_prompt: {str(e)}")
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Unhandled error: {str(e)}")


@router.post("/validate-scientific/", response_model=ValidationResponse)
async def validate_scientific(request: PromptRequest):
    """Endpoint to validate if a prompt is scientific in nature."""
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500, detail="OPENAI_API_KEY environment variable is not set"
            )

        # Create domain validator with appropriate model
        domain_validator = AgentFactory.create_domain_validator(request.model)

        validation_result = domain_validator.is_molecular(request.prompt)

        return {"is_molecular": validation_result.is_true}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# System Prompt for Diagram Generation - Now a template
DIAGRAM_SYSTEM_PROMPT_TEMPLATE = """
You are a helpful AI assistant that generates plans for 2D chemical diagrams based on user prompts.
Your goal is to return a JSON object that strictly conforms to the following Pydantic schema. Do NOT include any other text, explanations, or markdown formatting around the JSON.

The target canvas for the diagram has a width of {{canvas_width}} pixels and a height of {{canvas_height}} pixels.
All coordinates (x, y) and dimensions (width, height) you provide for molecules and arrows MUST be chosen to fit entirely within these canvas boundaries. Assume (0,0) is the top-left corner.

Expected JSON Schema:
```json
{{
  "plan": "string (A brief explanation of the diagram layout and chemical process)",
  "molecule_list": [
    {{
      "molecule": "string (MUST be a valid common chemical name, IUPAC name, or chemical formula, e.g., H2O, CH4, water, methane. Do NOT use placeholders like 'N/A' or generic terms if a specific chemical is implied by the user prompt.)",
      "x": "float (x-coordinate for top-left of molecule's bounding box, must be >= 0 and < {{canvas_width}})",
      "y": "float (y-coordinate for top-left of molecule's bounding box, must be >= 0 and < {{canvas_height}})",
      "width": "float (optional, width of bounding box, e.g., 100.0; x + width must be <= {{canvas_width}})",
      "height": "float (optional, height of bounding box, e.g., 100.0; y + height must be <= {{canvas_height}})",
      "label": "string (optional, label for the molecule, e.g., 2H₂O. If omitted, the 'molecule' value will be used as the label.)",
      "label_position": "string (optional, one of: above, below, left, right; defaults to below)"
    }}
    // ... more molecules
  ],
  "arrows": [
    {{
      "start": "[float, float] (x,y coordinates for arrow start, within canvas bounds)",
      "end": "[float, float] (x,y coordinates for arrow end, within canvas bounds)",
      "style": "string (optional, one of: straight, curved; defaults to straight)",
      "text": "string (optional, text label for the arrow)"
    }}
    // ... more arrows
  ]
}}
```

Example of a valid JSON response for a prompt like "Show water H2O decomposing into Hydrogen H2 and Oxygen O2" on a 600x400 canvas:
```json
{{
  "plan": "Water (H2O) on the left decomposes into Hydrogen (H2) and Oxygen (O2) on the right, separated by a reaction arrow, fitting within a 600x400 canvas.",
  "molecule_list": [
    {{"molecule": "H2O", "x": 50.0, "y": 150.0, "width": 100.0, "height": 100.0, "label": "H₂O"}},
    {{"molecule": "H2", "x": 350.0, "y": 100.0, "width": 100.0, "height": 100.0, "label": "H₂"}},
    {{"molecule": "O2", "x": 350.0, "y": 200.0, "width": 100.0, "height": 100.0, "label": "O₂"}}
  ],
  "arrows": [
    {{"start": [180.0, 175.0], "end": [320.0, 175.0], "text": "reaction"}}
  ]
}}
```

Ensure all field names and types match exactly. For the 'molecule' field, provide specific, recognizable chemical identifiers. Floats should be numbers (e.g., 150.0). Optional fields can be omitted if not applicable.
User prompt will follow.
"""


@router.post("/generate-molecule-diagram/", response_model=DiagramResponse)
async def generate_molecule_diagram(request: DiagramPromptRequest) -> DiagramResponse:
    """Generates a 2-D molecule diagram from a text prompt.

    1. Uses an LLM to create a plan (molecules, positions, arrows).
    2. Fetches 2D molecule data from PubChem.
    3. Renders an SVG diagram.
    """
    # Add validation for empty prompt
    if not request.prompt.strip():
        raise HTTPException(status_code=422, detail="Prompt cannot be empty.")

    try:

        # Format the system prompt with actual canvas dimensions
        formatted_system_prompt = DIAGRAM_SYSTEM_PROMPT_TEMPLATE.replace(
            "{{canvas_width}}", str(request.canvas_width)
        ).replace("{{canvas_height}}", str(request.canvas_height))
        # Escape curly braces for JSON schema examples within the f-string formatted prompt if any remain
        # For this template, direct replacement is fine since schema examples use double curlies.

        diagram_plan_from_llm = generate_structured(
            user_prompt=request.prompt,
            response_model=DiagramPlan,
            system_prompt=formatted_system_prompt,
        )

        if not diagram_plan_from_llm or not diagram_plan_from_llm.molecule_list:
            raise ValueError(
                "LLM failed to generate a valid diagram plan or molecule list."
            )

        # Augment the plan with canvas dimensions from the request
        # We create a new DiagramPlan instance or update the existing one if mutable
        # For simplicity, let's assume we'll construct the final plan dict for the response
        # or update the pydantic model if it's easy.

        # Create a dictionary from the LLM's plan and add/overwrite canvas dimensions
        final_diagram_plan_dict = diagram_plan_from_llm.model_dump()
        final_diagram_plan_dict["canvas_width"] = request.canvas_width
        final_diagram_plan_dict["canvas_height"] = request.canvas_height

        # Validate this potentially augmented dict back into a DiagramPlan object if needed for type consistency downstream
        # or ensure DiagramResponse can take this dict for its diagram_plan field.
        # If DiagramResponse.diagram_plan expects a DiagramPlan object:
        final_diagram_plan_obj = DiagramPlan(**final_diagram_plan_dict)

        layout_requests = []
        for mp in final_diagram_plan_obj.molecule_list:  # Use the augmented plan
            layout_requests.append(
                {
                    "query": mp.molecule,
                    "box": {
                        "x": mp.x,
                        "y": mp.y,
                        "width": mp.width or 200,
                        "height": mp.height or 200,
                    },
                }
            )

        # Prepare layout prompt and fetch bounding-box data via structured LLM
        layout_prompt = (
            "Given the following layout requests, return an array of objects "
            "conforming to the schema Molecule2DDataList:\n"
            + json.dumps(layout_requests)
        )
        fetched = generate_structured(
            user_prompt=layout_prompt,
            response_model=Molecule2DDataList,
            system_prompt="You are a helpful assistant that fills in bounding-box data for molecules.",
        )
        molecules = fetched.__root__

        if len(molecules) != len(final_diagram_plan_obj.molecule_list):
            raise ValueError(
                "Mismatch between planned molecules and fetched molecule data."
            )

        # Create renderable molecules list
        renderable_molecules = {}
        # Use LLM-provided bounding-box data to build renderable molecules
        for mol in molecules:
            renderable_molecules[mol.molecule] = {
                "x": mol.x,
                "y": mol.y,
                "width": mol.width,
                "height": mol.height,
                "label": getattr(mol, "label", mol.molecule) or mol.molecule,
                "label_position": getattr(mol, "label_position", "below"),
            }

        svg_image = render_diagram(
            plan=final_diagram_plan_obj,
            molecule_data=renderable_molecules,
            canvas_width=request.canvas_width,
            canvas_height=request.canvas_height,
        )

        return DiagramResponse(
            diagram_image=svg_image,
            diagram_plan=final_diagram_plan_obj,
            status="completed",
        )

    except HTTPException:
        raise
    except ValueError as ve:
        logger.error(
            f"ValueError in generate_molecule_diagram: {str(ve)}\n{traceback.format_exc()}"
        )
        # Include canvas dimensions in error plan as well
        error_plan = DiagramPlan(
            plan=f"ValueError: {str(ve)}",
            molecule_list=[],
            arrows=[],
            canvas_width=request.canvas_width,  # Add here
            canvas_height=request.canvas_height,  # Add here
        )
        return DiagramResponse(
            diagram_image="",
            diagram_plan=error_plan,
            status="failed",
            error=str(ve),
        )
    except Exception as e:
        logger.error(
            f"Unexpected error in generate_molecule_diagram: {str(e)}\n{traceback.format_exc()}"
        )
        raise HTTPException(status_code=500, detail=str(e))


class FetchMoleculeRequest(BaseModel):
    query: str


class Box2D(BaseModel):
    """2-D bounding box for placing a molecule."""

    x: float
    y: float
    width: float
    height: float


class MoleculeBoxRequest(BaseModel):
    """Pair a molecule query with a placement box."""

    query: str
    box: Box2D


class MoleculeLayoutRequest(BaseModel):
    """Request multiple molecules with layout information."""

    molecules: List[MoleculeBoxRequest]


class GenerateHTMLRequest(BaseModel):
    molecule_data: Dict[str, Any]


class SDFToPDBRequest(BaseModel):
    """Input SDF text to convert to PDB."""

    sdf: str


class SDFToPDBResponse(BaseModel):
    pdb_data: str


@router.post("/fetch-molecule-data/")
async def fetch_molecule_data(request: FetchMoleculeRequest):
    """
    Step A: Return molecule info + PDB block (no HTML).
    Frontend can store or pass it later to generate HTML.
    """
    try:
        # Create or get the PubChemAgent
        pubchem_agent = AgentFactory.create_pubchem_agent(
            script_model=None, use_element_labels=True, convert_back_to_indices=True
        )
        data = pubchem_agent.get_molecule_data(request.query)
        return data
    except Exception as e:
        logger.error(f"Error in fetch_molecule_data: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/fetch-molecule-2d/")
async def fetch_molecule_2d_data(request: FetchMoleculeRequest):
    """Return 2D coordinate information for a molecule."""
    try:
        pubchem_agent = AgentFactory.create_pubchem_agent(
            script_model=None,
            use_element_labels=True,
            convert_back_to_indices=True,
        )
        data = pubchem_agent.get_molecule_2d_info(request.query)
        return data
    except Exception as e:
        logger.error(f"Error in fetch_molecule_2d_data: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/fetch-molecule-layout/")
async def fetch_molecule_layout(request: MoleculeLayoutRequest):
    """Return 2D info for multiple molecules with layout boxes."""
    try:
        pubchem_agent = AgentFactory.create_pubchem_agent(
            script_model=None,
            use_element_labels=True,
            convert_back_to_indices=True,
        )
        queries = [
            {"query": m.query, "box": m.box.model_dump()} for m in request.molecules
        ]
        data = pubchem_agent.get_molecules_2d_layout(queries)
        return {"molecules": data}
    except Exception as e:
        logger.error(f"Error in fetch_molecule_layout: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/sdf-to-pdb/", response_model=SDFToPDBResponse)
async def convert_sdf_to_pdb(request: SDFToPDBRequest):
    """Convert SDF text to PDB format using RDKit."""
    try:
        pdb_data = _sdf_to_pdb_block(request.sdf)
        if not pdb_data:
            raise ValueError("Failed to convert SDF to PDB")
        return {"pdb_data": pdb_data}
    except Exception as e:
        logger.error(f"Error in convert_sdf_to_pdb: {str(e)}")
        raise HTTPException(status_code=400, detail=str(e))
