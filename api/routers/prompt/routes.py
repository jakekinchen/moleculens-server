"""Routes for handling prompts and generating visualizations."""

import logging
import os
import traceback
from datetime import datetime
from typing import Any, Dict, List, Literal, Optional, Tuple

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException
from pydantic import BaseModel, Field

from api.agent_management.agent_factory import AgentFactory
from api.agent_management.agents.pubchem_agent import PubChemAgent, _sdf_to_pdb_block
from api.agent_management.llm_service import (
    LLMModelConfig,
    LLMService,
    StructuredLLMRequest,
)
from api.agent_management.model_config import ModelCategory, get_default_model
from api.agent_management.model_registry import ModelRegistry
from api.agent_management.models import (
    AnimationCode,
    BaseModelWithConfig,
    OrchestrationPlan,
    SceneScript,
)
from api.agent_management.scene_packager import ScenePackager
from api.dependencies.use_llm import use_llm
from api.diagram.diagram_renderer import render_diagram
from api.diagram.models import DiagramPlan

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
    preferred_model_category: Optional[str] = None


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


class ValidationResponse(BaseModel):
    is_molecular: bool


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

    # Add configuration settings for better validation handling
    model_config = {
        "protected_namespaces": (),
        "arbitrary_types_allowed": True,
        "json_schema_extra": {
            "examples": [
                {
                    "job_id": "b96a7408-d3b0-4022-9a98-79e71a798be9",
                    "status": "processing",
                    "message": "Processing in progress",
                    "progress": 0.5,
                    "result": "",
                }
            ]
        },
    }


# Import ModelRegistry and related functions
from api.agent_management.model_config import get_default_model

# Simple in-memory job stores
# In a production app
geometry_jobs: Dict[str, Dict[str, Any]] = {}
pipeline_jobs: Dict[str, Dict[str, Any]] = {}


# Background task function for processing prompts
async def process_prompt_pipeline_task(
    job_id: str, prompt: str, override_model: Optional[str] = None
):
    """Background task to process a prompt through the entire pipeline. Updates
    the pipeline_jobs dictionary with the result when complete.

    Args:
        job_id: Unique identifier for this job
        prompt: The user's prompt to process
        override_model: Optional model name to override all agents
    """
    try:
        print(f"[Job {job_id}] Starting processing for prompt: {prompt[:50]}...")

        # Initialize all agents with appropriate models
        try:
            # Create agents using the factory with optimal model selection
            # Validation already happened in the main endpoint
            script_agent = AgentFactory.create_script_agent(override_model)
            orchestration_agent = AgentFactory.create_orchestration_agent(
                override_model
            )
            animation_agent = AgentFactory.create_animation_agent(override_model)

            print(f"[Job {job_id}] All agents initialized successfully")

            # Update progress after initialization
            pipeline_jobs[job_id]["progress"] = 0.1
        except Exception as e:
            print(f"[Job {job_id}] Error initializing agents: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error initializing agents: {str(e)}",
            }
            return

        # Step 1: Generate an animation script (validation already done)
        try:
            print(f"[Job {job_id}] Generating animation script...")
            animation_script: SceneScript = script_agent.generate_script(prompt)
            print(
                f"[Job {job_id}] Script generated with {len(animation_script.content)} time points"
            )

            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.3
        except Exception as e:
            print(f"[Job {job_id}] Error generating script: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating script: {str(e)}",
            }
            return

        # Step 2: Generate an orchestration plan based on the script
        try:
            print(f"[Job {job_id}] Generating orchestration plan...")
            orchestration_plan = orchestration_agent.generate_orchestration_plan(
                animation_script
            )
            print(
                f"[Job {job_id}] Orchestration plan generated with {len(orchestration_plan.objects)} objects"
            )

            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.4
        except Exception as e:
            print(f"[Job {job_id}] Error generating orchestration plan: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating orchestration plan: {str(e)}",
            }
            return

        # Step 3: Generate geometry for all objects in the plan
        try:
            print(f"[Job {job_id}] Generating geometry for all objects...")
            object_geometries = await orchestration_agent.generate_geometry_from_plan(
                orchestration_plan
            )

            # Check if there are any successful geometries
            successful_count = sum(
                1
                for obj_data in object_geometries.values()
                if obj_data.get("status") == "success" and obj_data != "_summary"
            )
            print(f"[Job {job_id}] Generated {successful_count} successful geometries")

            if successful_count == 0:
                print(
                    f"[Job {job_id}] Warning: No successful geometries were generated"
                )

            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.7
        except Exception as e:
            print(f"[Job {job_id}] Error generating geometries: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating geometries: {str(e)}",
            }
            return

        # Step 4: Generate animation code
        try:
            print(f"[Job {job_id}] Generating animation code...")
            animation_code = animation_agent.generate_animation_code(
                script=animation_script,
                object_geometries=object_geometries,
                orchestration_plan=orchestration_plan,
            )
            print(
                f"[Job {job_id}] Animation code generated with {len(animation_code.keyframes)} keyframes"
            )

            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.8
        except Exception as e:
            print(f"[Job {job_id}] Error generating animation code: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating animation code: {str(e)}",
            }
            return

        # Step 5: Package everything into a complete scene
        try:
            print(f"[Job {job_id}] Packaging scene...")
            scene_package = ScenePackager.create_scene_package(
                script=animation_script,
                orchestration_plan=orchestration_plan,
                object_geometries=object_geometries,
                animation_code=animation_code,
            )
            print(f"[Job {job_id}] Scene packaged successfully")

            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.9
        except Exception as e:
            print(f"[Job {job_id}] Error packaging scene: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error packaging scene: {str(e)}",
            }
            return

        # Save the JS and HTML files to the static directory
        try:
            print(f"[Job {job_id}] Saving JS and HTML files to static directory...")
            static_dir = os.path.join(
                os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "static"
            )
            os.makedirs(static_dir, exist_ok=True)

            # Generate filenames based on job ID
            js_filename = f"scene_{job_id}.js"
            html_filename = f"scene_{job_id}.html"

            # Write the JavaScript to a file (useful for debugging/development)
            with open(os.path.join(static_dir, js_filename), "w") as f:
                f.write(scene_package.js)

            # Write the HTML file with embedded JavaScript
            with open(os.path.join(static_dir, html_filename), "w") as f:
                f.write(scene_package.html)

            print(
                f"[Job {job_id}] Files saved successfully: JS ({js_filename}), HTML with embedded JS ({html_filename})"
            )

            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.95
        except Exception as e:
            print(f"[Job {job_id}] Error saving JS/HTML files: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error saving JS/HTML files: {str(e)}",
            }
            return

        # Update job with the completed result
        print(f"[Job {job_id}] Processing complete, storing results")

        result = {
            "html": scene_package.html,
            "js": scene_package.js,
            "minimal_js": scene_package.minimal_js,
            "title": scene_package.title,
            "timecode_markers": scene_package.timecode_markers,
            "total_elements": scene_package.total_elements,
            "js_filename": js_filename,
        }

        pipeline_jobs[job_id] = {
            **pipeline_jobs[job_id],
            "status": "completed",
            "progress": 1.0,
            "result": result,
        }

    except Exception as e:
        print(
            f"[Job {job_id}] Unhandled exception in process_prompt_pipeline_task: {str(e)}"
        )
        traceback.print_exc()
        pipeline_jobs[job_id] = {
            **pipeline_jobs[job_id],
            "status": "error",
            "progress": pipeline_jobs[job_id].get("progress", 0),
            "error": f"Unhandled error: {str(e)}",
        }


@router.post("/process/", response_model=dict)
async def submit_prompt_background(
    request: PromptRequest, background_tasks: BackgroundTasks
):
    """Starts the scientific visualization pipeline as a background task and
    returns immediately. This endpoint begins processing the prompt through all
    pipeline steps but doesn't wait for completion.

    Returns a job ID that can be used to check the status of processing.

    Query parameters:
    - model: Optional specific model to use for all agents
    - preferred_model_category: Optional model category preference
    """
    # Generate a unique job ID
    import uuid

    job_id = str(uuid.uuid4())

    # Start the background task with the selected model (if specified)
    background_tasks.add_task(
        process_prompt_pipeline_task,
        job_id=job_id,
        prompt=request.prompt,
        override_model=request.model,
    )

    # Create job entry
    pipeline_jobs[job_id] = {
        "status": "processing",
        "prompt": request.prompt,
        "created_at": datetime.now().isoformat(),
        "result": None,
        "progress": 0.0,  # Initial progress
    }

    # Return immediately with the job ID
    return {
        "job_id": job_id,
        "status": "processing",
        "message": "Processing started in the background",
        "progress": 0.0,
    }


@router.get("/process/{job_id}")
async def check_process_status(job_id: str):
    """Check the status of a background processing job. Returns the current
    status and result if processing is complete.

    Response format:
    - job_id: The unique job ID
    - status: 'processing', 'completed', or 'failed'
    - progress: A float between 0 and 1 indicating progress
    - message: A human-readable status message
    - visualization: (When completed) The full visualization data
    - error: (When failed) The error message

    The visualization object contains:
    - html: The complete HTML for the visualization with embedded JavaScript
    - title: The title of the scene
    - timecode_markers: List of timecodes for the animation
    - total_elements: Total number of elements in the scene
    """
    if job_id not in pipeline_jobs:
        raise HTTPException(status_code=404, detail=f"Job ID {job_id} not found")

    job = pipeline_jobs[job_id]

    # If the job is completed, return the full result
    if job["status"] == "completed":
        # Convert the result to the VisualizationData format
        visualization = None
        if job["result"]:
            visualization = {
                "html": job["result"].get("html", ""),
                "title": job["result"].get("title", "Scientific Visualization"),
                "timecode_markers": job["result"].get("timecode_markers", []),
                "total_elements": job["result"].get("total_elements", 0),
            }

        return {
            "job_id": job_id,
            "status": "completed",
            "progress": 1.0,
            "message": "Visualization processing completed successfully",
            "visualization": visualization,
            # Legacy fields for backward compatibility - convert to string for compatibility
            "result": str(job["result"].get("js", "")) if job["result"] else "",
            "geometry_result": job["result"].get("js", "") if job["result"] else "",
        }
    # If there was an error, include the error message
    elif job["status"] == "error":
        return {
            "job_id": job_id,
            "status": "failed",
            "progress": job.get("progress", 0.0),
            "message": "Processing failed",
            "error": job.get("error", "Unknown error occurred"),
            # Add empty result field for compatibility
            "result": "",
        }
    # Otherwise just return the status info
    else:
        return {
            "job_id": job_id,
            "status": "processing",
            "progress": job.get("progress", 0.0),
            "message": "Processing in progress",
            # Add empty result field for compatibility
            "result": "",
        }


@router.post("/generate-from-pubchem/", response_model=dict)
async def generate_from_pubchem(request: PromptRequest):
    """Generate a 3D visualization from a query to PubChem.

    Returns:
    - pdb_data: PDB data for the molecule
    - result_html: HTML content
    - title: Molecule title
    """
    try:
        logger.info(f"Processing PubChem request: {request.prompt}")

        # Create domain validator agent using the factory
        domain_validator = AgentFactory.create_domain_validator(request.model)

        # Initial validation can be done synchronously to reject bad prompts immediately
        validation_result = domain_validator.is_molecular(request.prompt)

        # If not molecular, reject the prompt immediately
        if not validation_result.is_true:
            logger.warning(f"Non-molecular prompt rejected: {request.prompt}")
            return {
                "job_id": "rejected",
                "status": "failed",
                "message": "Non-molecular prompt rejected",
                "error": "The prompt does not contain molecular content",
            }

        # Create PubChem agent with script agent override - use the specified model for script generation
        pubchem_agent = AgentFactory.create_pubchem_agent(
            script_model=request.model,
            use_element_labels=True,  # Use element-based labels for better LLM understanding during script generation
            convert_back_to_indices=True,  # ALWAYS convert element-based labels back to numeric indices for visualization
        )

        # Generate the geometry directly for immediate response
        logger.info(f"Generating molecule package for: {request.prompt}")
        result = pubchem_agent.get_molecule_package(request.prompt)
        logger.info(f"Successfully generated molecule package: {result.title}")

        return {
            "pdb_data": result.pdb_data,
            "result_html": result.html,
            "title": result.title,
        }
    except Exception as e:
        error_message = str(e)
        logger.error(f"Error generating from PubChem: {error_message}")
        logger.error(traceback.format_exc())

        # Provide a more user-friendly error message
        if (
            "Anthropic API error: " in error_message
            and "overloaded" in error_message.lower()
        ):
            error_message = "The service is currently experiencing high demand. Please try again in a few moments."
        elif "No molecules found" in error_message:
            error_message = f"Could not find a molecule matching your query. Please try a different molecule name or formula."

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
        # Create domain validator agent using the factory
        domain_validator = AgentFactory.create_domain_validator(request.model)

        # Initial validation can be done synchronously to reject bad prompts immediately
        validation_result = domain_validator.is_molecular(request.prompt)

        # If not scientific, reject the prompt immediately
        if not validation_result.is_true:
            return {
                "job_id": "rejected",
                "status": "failed",
                "message": "Non-molecular prompt rejected",
                "error": "The prompt does not contain molecular content",
            }

        # Generate a unique job ID
        import uuid

        job_id = str(uuid.uuid4())

        # Start the background task with the selected model (if specified)
        background_tasks.add_task(
            process_prompt_pipeline_task,
            job_id=job_id,
            prompt=request.prompt,
            override_model=request.model,
        )

        # Create job entry
        pipeline_jobs[job_id] = {
            "status": "processing",
            "prompt": request.prompt,
            "created_at": datetime.now().isoformat(),
            "result": None,
            "progress": 0.0,  # Initial progress
        }

        print(
            f"Started background processing job {job_id} for prompt: {request.prompt[:50]}..."
        )

        # Return immediately with the job ID and status
        return {
            "job_id": job_id,
            "status": "processing",
            "message": "Processing started in the background. Poll /prompt/process/{job_id} for updates.",
            "progress": 0.0,
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


def _create_pubchem_agent(model: Optional[str] = None) -> PubChemAgent:
    return AgentFactory.create_pubchem_agent(
        override_model=model,
        script_model=None,
        use_element_labels=True,
        convert_back_to_indices=True,
    )


@router.post("/generate-molecule-diagram/", response_model=DiagramResponse)
async def generate_molecule_diagram(
    request: DiagramPromptRequest, llm_service: LLMService = Depends(use_llm)
) -> DiagramResponse:
    """Generates a 2-D molecule diagram from a text prompt.

    1. Uses an LLM to create a plan (molecules, positions, arrows).
    2. Fetches 2D molecule data from PubChem.
    3. Renders an SVG diagram.
    """
    # Add validation for empty prompt
    if not request.prompt.strip():
        raise HTTPException(status_code=422, detail="Prompt cannot be empty.")

    try:
        model_info_instance = ModelRegistry.create_instance(
            request.model or get_default_model()
        )

        llm_model_config = LLMModelConfig(
            provider=model_info_instance.provider,
            model_name=model_info_instance.name,
        )

        # Format the system prompt with actual canvas dimensions
        formatted_system_prompt = DIAGRAM_SYSTEM_PROMPT_TEMPLATE.replace(
            "{{canvas_width}}", str(request.canvas_width)
        ).replace("{{canvas_height}}", str(request.canvas_height))
        # Escape curly braces for JSON schema examples within the f-string formatted prompt if any remain
        # For this template, direct replacement is fine since schema examples use double curlies.

        structured_llm_request = StructuredLLMRequest[DiagramPlan](
            user_prompt=request.prompt,
            system_prompt=formatted_system_prompt,  # Use the formatted prompt
            response_model=DiagramPlan,
            llm_config=llm_model_config,
        )

        diagram_plan_from_llm = llm_service.generate_structured(structured_llm_request)

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

        pubchem_agent = _create_pubchem_agent(request.model)
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
        fetched_molecules_data = pubchem_agent.get_molecules_2d_layout(layout_requests)

        if len(fetched_molecules_data) != len(final_diagram_plan_obj.molecule_list):
            raise ValueError(
                "Mismatch between planned molecules and fetched molecule data."
            )

        # Create renderable molecules list
        renderable_molecules = {}
        for plan_data_from_list in final_diagram_plan_obj.molecule_list:
            name = plan_data_from_list.molecule
            fetched_data = await fetch_molecule_2d_data(
                FetchMoleculeRequest(query=name)
            )
            renderable_molecules[name] = {
                **fetched_data,
                "label": plan_data_from_list.label or plan_data_from_list.molecule,
                "label_position": plan_data_from_list.label_position,
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
        # When raising HTTPException for a 500, the detail is the error message.
        # The DiagramPlan is not directly part of the 500 response body in this case.
        raise HTTPException(
            status_code=500, detail=f"Unexpected server error: {str(e)}"
        )


@router.post("/generate-geometry/", response_model=GeometryResponse)
async def generate_geometry(
    request: GeometryRequest, llm_service: LLMService = Depends(use_llm)
):
    """Endpoint to generate Three.js geometry based on user prompt. Only
    generates geometry for scientific content.

    This endpoint processes the request immediately and returns the result.
    For longer processing, use the /prompt/process/ endpoint.

    Query parameters:
    - model: Optional global model override for all agents
    - preferred_model_category: Optional model category preference
    - use_case: Set to 'geometry' to use the geometry-specific default model
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500, detail="OPENAI_API_KEY environment variable is not set"
            )

        # Get override models from request if specified
        global_override_model = request.model

        # Create domain validator agent using the factory with potential override
        domain_validator = AgentFactory.create_domain_validator(global_override_model)

        # Validate the prompt is scientific
        validation_result = domain_validator.is_molecular(request.prompt)

        if not validation_result.is_true:
            # Instead of raising an error, return a response indicating non-molecular content
            return {
                "result": "",
                "is_molecular": False,
                "validation_message": "The prompt does not contain molecular content",
            }

        # Create geometry agent with appropriate model - use specific override if provided
        geometry_agent = AgentFactory.create_geometry_agent(global_override_model)

        # Generate the geometry directly for immediate response
        generated_code = geometry_agent.get_geometry_snippet(request.prompt)

        return {
            "result": generated_code,
            "is_molecular": True,
            "validation_message": None,
        }
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/models/", response_model=List[Dict[str, Any]])
async def get_models():
    """Endpoint to get information about all available models.

    Returns a list of registered models with their capabilities.
    """
    models = []
    for model_name in ModelRegistry.list_models():
        model_info = ModelRegistry.create_instance(model_name)
        models.append(
            {
                "name": model_name,
                "display_name": model_info.display_name,
                "provider": model_info.provider,
                "categories": [category for category in model_info.categories],
                "context_length": model_info.context_length,
                "is_default": model_info.is_default,
                "default_for": model_info.default_for,
            }
        )
    return models


@router.get("/agent-models/", response_model=List[Dict[str, Any]])
async def get_agent_model_configs():
    """Endpoint to get information about all agent-model configurations.

    Returns a list of agents with their preferred models.
    """
    from api.agent_management.agent_model_config import AGENT_MODEL_MAP, AgentType

    configs = []
    for agent_type in AgentType:
        if agent_type in AGENT_MODEL_MAP:
            config = AGENT_MODEL_MAP[agent_type]
            configs.append(
                {
                    "agent_type": agent_type,
                    "preferred_model": config.preferred_model,
                    "fallback_models": config.fallback_models,
                    "required_categories": [
                        category for category in config.required_categories
                    ],
                    "description": config.description,
                }
            )

    return configs


class ScriptRequest(BaseModel):
    script: SceneScript


class OrchestrationRequest(BaseModel):
    plan: OrchestrationPlan


class GeometryGenerationResponse(BaseModel):
    job_id: str
    status: str
    message: str


class GeometryResultResponse(BaseModel):
    job_id: str
    status: str
    completed: int
    total: int
    results: Optional[Dict[str, Any]] = None


class AnimationRequest(BaseModel):
    script: SceneScript
    orchestration_plan: OrchestrationPlan
    object_geometries: Dict[str, Any]


class AnimationResponse(BaseModel):
    code: str
    keyframes: List[Dict[str, Any]]


class PackagedSceneRequest(BaseModelWithConfig):
    script: SceneScript
    orchestration_plan: OrchestrationPlan
    object_geometries: Dict[str, Any]
    animation_code: AnimationCode


class PackagedSceneResponse(BaseModelWithConfig):
    html: str
    js: str
    minimal_js: str
    title: str
    timecode_markers: List[str]
    total_elements: int


@router.post("/generate-orchestration/", response_model=OrchestrationPlan)
async def generate_orchestration(
    request: ScriptRequest,
    model: Optional[str] = None,  # Add query parameter for model override
):
    """Endpoint to generate an orchestration plan from a scene script. Breaks
    down the script into discrete objects needed for the visualization.

    Query parameters:
    - model: Optional specific model to use for orchestration agent
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500, detail="OPENAI_API_KEY environment variable is not set"
            )

        # Create orchestration agent with appropriate model (potentially overriding with query param)
        orchestration_agent = AgentFactory.create_orchestration_agent(model)

        # Generate the orchestration plan from the script
        orchestration_plan = orchestration_agent.generate_orchestration_plan(
            request.script
        )
        return orchestration_plan
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# Background task to generate geometry for all objects in the plan
async def generate_geometry_task(
    job_id: str, plan: OrchestrationPlan, override_model: Optional[str] = None
):
    """Background task to generate geometry for all objects in the plan.

    Args:
        job_id: The unique job identifier
        plan: The orchestration plan to process
        override_model: Optional model name to override the default for the geometry agent
    """
    try:
        # Initialize job status
        geometry_jobs[job_id] = {
            "status": "processing",
            "total": len(plan.objects),
            "completed": 0,
            "results": None,
        }

        # Create orchestration agent with appropriate model
        orchestration_agent = AgentFactory.create_orchestration_agent(override_model)

        # Generate geometry for all objects
        results = await orchestration_agent.generate_geometry_from_plan(plan)

        # Update job status
        geometry_jobs[job_id]["status"] = "completed"
        geometry_jobs[job_id]["completed"] = len(plan.objects)
        geometry_jobs[job_id]["results"] = results

    except Exception as e:
        # Update job status in case of error
        geometry_jobs[job_id]["status"] = "failed"
        geometry_jobs[job_id]["error"] = str(e)


@router.post("/generate-geometry-for-plan/", response_model=GeometryGenerationResponse)
async def generate_geometry_for_plan(
    request: OrchestrationRequest,
    background_tasks: BackgroundTasks,
    model: Optional[str] = None,  # Add query parameter for model override
):
    """Endpoint to generate Three.js geometry for all objects in an
    orchestration plan. This starts an asynchronous job that processes each
    object sequentially.

    Query parameters:
    - model: Optional specific model to use for geometry generation
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500, detail="OPENAI_API_KEY environment variable is not set"
            )

        # Generate a unique job ID
        import uuid

        job_id = str(uuid.uuid4())

        # Start the background task to generate geometry with potential model override
        background_tasks.add_task(generate_geometry_task, job_id, request.plan, model)

        # Return the job ID so the client can check status
        return {
            "job_id": job_id,
            "status": "processing",
            "message": f"Started processing {len(request.plan.objects)} objects",
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/geometry-job-status/{job_id}", response_model=GeometryResultResponse)
async def get_geometry_job_status(job_id: str):
    """Check the status of a geometry generation job.

    Returns the results if the job is completed.
    """
    if job_id not in geometry_jobs:
        raise HTTPException(status_code=404, detail=f"Job ID {job_id} not found")

    job = geometry_jobs[job_id]

    return {
        "job_id": job_id,
        "status": job["status"],
        "completed": job.get("completed", 0),
        "total": job.get("total", 0),
        "results": job.get("results", None),
    }


@router.post("/generate-animation/", response_model=AnimationResponse)
async def generate_animation(
    request: AnimationRequest,
    model: Optional[str] = None,  # Add query parameter for model override
):
    """Generate animation code for a scene based on the script, orchestration
    plan, and generated geometries. This endpoint takes the outputs from
    previous steps in the pipeline.

    Query parameters:
    - model: Optional specific model to use for animation code generation
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500, detail="OPENAI_API_KEY environment variable is not set"
            )

        # Create animation agent with appropriate model (potentially overriding with query param)
        animation_agent = AgentFactory.create_animation_agent(model)

        # Generate animation code
        animation = animation_agent.generate_animation_code(
            script=request.script,
            object_geometries=request.object_geometries,
            orchestration_plan=request.orchestration_plan,
        )

        # Convert keyframes to serializable format
        keyframes_dict = [
            {"timecode": keyframe.timecode, "actions": keyframe.actions}
            for keyframe in animation.keyframes
        ]

        return {"code": animation.code, "keyframes": keyframes_dict}
    except Exception as e:
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


@router.post("/generate-molecule-html/")
async def generate_molecule_html(request: GenerateHTMLRequest):
    """
    Step B: Given existing molecule data, generate script + HTML,
    without re-downloading from PubChem or re-calling earlier steps.
    """
    try:
        pubchem_agent = AgentFactory.create_pubchem_agent(
            script_model=None, use_element_labels=True, convert_back_to_indices=True
        )
        html = pubchem_agent.generate_visualization(request.molecule_data)
        return {"html": html}
    except Exception as e:
        logger.error(f"Error generating HTML: {str(e)}")
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


@router.post("/package-scene/", response_model=PackagedSceneResponse)
async def package_scene(request: PackagedSceneRequest):
    """Create a complete, packaged Three.js scene from the generated
    components.

    This endpoint combines all outputs from previous steps into a ready-
    to-use scene.
    """
    try:
        # Use the ScenePackager to create a complete scene
        scene_package = ScenePackager.create_scene_package(
            script=request.script,
            orchestration_plan=request.orchestration_plan,
            object_geometries=request.object_geometries,
            animation_code=request.animation_code,
        )

        # Save the JS and HTML files to the static directory
        import os

        static_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "static"
        )
        os.makedirs(static_dir, exist_ok=True)

        # Generate filenames with timestamp to avoid conflicts
        import datetime

        timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        js_filename = f"scene_{timestamp}.js"
        html_filename = f"scene_{timestamp}.html"

        # Write the JavaScript to a file
        with open(os.path.join(static_dir, js_filename), "w") as f:
            f.write(scene_package.js)

        # Also save as scene.js for backward compatibility
        with open(os.path.join(static_dir, "scene.js"), "w") as f:
            f.write(scene_package.js)

        # Write the HTML file with embedded JavaScript
        with open(os.path.join(static_dir, html_filename), "w") as f:
            f.write(scene_package.html)

        return {
            "html": scene_package.html,
            "js": scene_package.js,
            "minimal_js": scene_package.minimal_js,
            "title": scene_package.title,
            "timecode_markers": scene_package.timecode_markers,
            "total_elements": scene_package.total_elements,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
