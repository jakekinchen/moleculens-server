from fastapi import APIRouter, HTTPException, Depends, BackgroundTasks
from pydantic import BaseModel
from typing import Optional, Dict, Any, List
from dependencies.use_llm import use_llm 
from agent_management.agents.geometry_agent import GeometryAgent
from agent_management.agents.domain_bool_agent import DomainValidator
from agent_management.agents.script_agent import ScriptAgent
from agent_management.agents.orchestration_agent import OrchestrationAgent
from agent_management.agents.animation_agent import AnimationAgent
from agent_management.models import SceneScript, ScriptTimePoint, OrchestrationPlan, AnimationCode, FinalScenePackage, BaseModelWithConfig
from agent_management.scene_packager import ScenePackager
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
import os
import asyncio
import traceback
from datetime import datetime


router = APIRouter(
    prefix="/prompt",
    tags=["Prompt"],
    responses={404: {"description": "Not found"}},
)


class PromptRequest(BaseModel):
    prompt: str

class PromptWorking(BaseModel):
    prompt: str
    working: bool
    script: Optional[SceneScript] = None
    orchestration_plan: Optional[OrchestrationPlan] = None
    
class CompletePipelineResponse(BaseModelWithConfig):
    """Response model for the complete prompt-to-visualization pipeline"""
    html: str
    js: str
    minimal_js: str
    title: str
    timecode_markers: List[str]
    total_elements: int

class GeometryResponse(BaseModel):
    result: str

class GeometryRequest(BaseModel):
    prompt: str

class ValidationResponse(BaseModel):
    is_scientific: bool
    confidence: float
    reasoning: str
    
class VisualizationData(BaseModel):
    html: str
    js: str
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
                    "result": ""
                }
            ]
        }
    }


# Initialize LLM config
llm_config = LLMModelConfig(
    provider=ProviderType.OPENAI,
    model_name="o3-mini",
    api_key=os.getenv("OPENAI_API_KEY")
)

# Simple in-memory job stores
# In a production app, these would use Redis or another persistent store
geometry_jobs = {}
pipeline_jobs = {}

# Background task function for processing prompts
async def process_prompt_pipeline_task(job_id: str, prompt: str, llm_service: LLMService):
    """
    Background task to process a prompt through the entire pipeline.
    Updates the pipeline_jobs dictionary with the result when complete.
    """
    try:
        print(f"[Job {job_id}] Starting processing for prompt: {prompt[:50]}...")
        
        # Initialize all agents
        try:
            # Validation already happened in the main endpoint
            # domain_validator = DomainValidator(llm_service)
            script_agent = ScriptAgent(llm_service)
            orchestration_agent = OrchestrationAgent(llm_service)
            animation_agent = AnimationAgent(llm_service)
            print(f"[Job {job_id}] All agents initialized successfully")
            
            # Update progress after initialization
            pipeline_jobs[job_id]["progress"] = 0.1
        except Exception as e:
            print(f"[Job {job_id}] Error initializing agents: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error initializing agents: {str(e)}"
            }
            return
        
        # Step 1: Generate an animation script (validation already done)
        try:
            print(f"[Job {job_id}] Generating animation script...")
            animation_script = script_agent.generate_script(prompt)
            print(f"[Job {job_id}] Script generated with {len(animation_script.content)} time points")
            
            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.3
        except Exception as e:
            print(f"[Job {job_id}] Error generating script: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating script: {str(e)}"
            }
            return
        
        # Step 2: Generate an orchestration plan based on the script
        try:
            print(f"[Job {job_id}] Generating orchestration plan...")
            orchestration_plan = orchestration_agent.generate_orchestration_plan(animation_script)
            print(f"[Job {job_id}] Orchestration plan generated with {len(orchestration_plan.objects)} objects")
            
            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.4
        except Exception as e:
            print(f"[Job {job_id}] Error generating orchestration plan: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating orchestration plan: {str(e)}"
            }
            return
        
        # Step 3: Generate geometry for all objects in the plan
        try:
            print(f"[Job {job_id}] Generating geometry for all objects...")
            object_geometries = await orchestration_agent.generate_geometry_from_plan(orchestration_plan)
            
            # Check if there are any successful geometries
            successful_count = sum(1 for obj_data in object_geometries.values() 
                               if obj_data.get("status") == "success" and obj_data != "_summary")
            print(f"[Job {job_id}] Generated {successful_count} successful geometries")
            
            if successful_count == 0:
                print(f"[Job {job_id}] Warning: No successful geometries were generated")
                
            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.7
        except Exception as e:
            print(f"[Job {job_id}] Error generating geometries: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating geometries: {str(e)}"
            }
            return
        
        # Step 4: Generate animation code
        try:
            print(f"[Job {job_id}] Generating animation code...")
            animation_code = animation_agent.generate_animation_code(
                script=animation_script,
                object_geometries=object_geometries,
                orchestration_plan=orchestration_plan
            )
            print(f"[Job {job_id}] Animation code generated with {len(animation_code.keyframes)} keyframes")
            
            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.8
        except Exception as e:
            print(f"[Job {job_id}] Error generating animation code: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error generating animation code: {str(e)}"
            }
            return
        
        # Step 5: Package everything into a complete scene
        try:
            print(f"[Job {job_id}] Packaging scene...")
            scene_package = ScenePackager.create_scene_package(
                script=animation_script,
                orchestration_plan=orchestration_plan,
                object_geometries=object_geometries,
                animation_code=animation_code
            )
            print(f"[Job {job_id}] Scene packaged successfully")
            
            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.9
        except Exception as e:
            print(f"[Job {job_id}] Error packaging scene: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error packaging scene: {str(e)}"
            }
            return
        
        # Save the JS file to the static directory
        try:
            print(f"[Job {job_id}] Saving JS file to static directory...")
            static_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "static")
            os.makedirs(static_dir, exist_ok=True)
            
            # Generate a JS filename based on job ID
            js_filename = f"scene_{job_id}.js"
            
            # Write the JavaScript to a file
            with open(os.path.join(static_dir, js_filename), "w") as f:
                f.write(scene_package.js)
                
            # Update the HTML to reference the correct JS file
            html = scene_package.html.replace('/static/scene.js', f'/static/{js_filename}')
            print(f"[Job {job_id}] JS file saved successfully as {js_filename}")
            
            # Update progress
            pipeline_jobs[job_id]["progress"] = 0.95
        except Exception as e:
            print(f"[Job {job_id}] Error saving JS file: {str(e)}")
            pipeline_jobs[job_id] = {
                **pipeline_jobs[job_id],
                "status": "error",
                "error": f"Error saving JS file: {str(e)}"
            }
            return
        
        # Update job with the completed result
        print(f"[Job {job_id}] Processing complete, storing results")
        
        result = {
            "html": html,
            "js": scene_package.js,
            "minimal_js": scene_package.minimal_js,
            "title": scene_package.title,
            "timecode_markers": scene_package.timecode_markers,
            "total_elements": scene_package.total_elements,
            "js_filename": js_filename
        }
        
        pipeline_jobs[job_id] = {
            **pipeline_jobs[job_id],
            "status": "completed",
            "progress": 1.0,
            "result": result
        }
        
    except Exception as e:
        print(f"[Job {job_id}] Unhandled exception in process_prompt_pipeline_task: {str(e)}")
        traceback.print_exc()
        pipeline_jobs[job_id] = {
            **pipeline_jobs[job_id],
            "status": "error",
            "progress": pipeline_jobs[job_id].get("progress", 0),
            "error": f"Unhandled error: {str(e)}"
        }



@router.post("/process/", response_model=dict)
async def submit_prompt_background(
    request: PromptRequest,
    background_tasks: BackgroundTasks,
    llm_service=Depends(use_llm)
    ):
    """
    Starts the scientific visualization pipeline as a background task and returns immediately.
    This endpoint begins processing the prompt through all pipeline steps but doesn't wait for completion.
    
    Returns a job ID that can be used to check the status of processing.
    """
    # Generate a unique job ID
    import uuid
    job_id = str(uuid.uuid4())
    
    # Start the background task
    background_tasks.add_task(
        process_prompt_pipeline_task,
        job_id=job_id,
        prompt=request.prompt,
        llm_service=llm_service
    )
    
    # Create job entry
    pipeline_jobs[job_id] = {
        "status": "processing",
        "prompt": request.prompt,
        "created_at": datetime.now().isoformat(),
        "result": None
    }
    
    # Return immediately with the job ID
    return {
        "job_id": job_id,
        "status": "processing",
        "message": "Processing started in the background"
    }

@router.get("/process/{job_id}")
async def check_process_status(job_id: str):
    """
    Check the status of a background processing job.
    Returns the current status and result if processing is complete.
    
    Response format:
    - job_id: The unique job ID
    - status: 'processing', 'completed', or 'failed'
    - progress: A float between 0 and 1 indicating progress
    - message: A human-readable status message
    - visualization: (When completed) The full visualization data
    - error: (When failed) The error message
    
    The visualization object contains:
    - html: The complete HTML for the visualization
    - js: The full JavaScript code
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
                "js": job["result"].get("js", ""),
                "title": job["result"].get("title", "Scientific Visualization"),
                "timecode_markers": job["result"].get("timecode_markers", []),
                "total_elements": job["result"].get("total_elements", 0)
            }
        
        return {
            "job_id": job_id,
            "status": "completed",
            "progress": 1.0,
            "message": "Visualization processing completed successfully",
            "visualization": visualization,
            # Legacy fields for backward compatibility - convert to string for compatibility
            "result": str(job["result"].get("js", "")) if job["result"] else "",
            "geometry_result": job["result"].get("js", "") if job["result"] else ""
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
            "result": ""
        }
    # Otherwise just return the status info
    else:
        return {
            "job_id": job_id,
            "status": "processing",
            "progress": job.get("progress", 0.0),
            "message": "Processing in progress",
            # Add empty result field for compatibility
            "result": ""
        }

@router.post("/", response_model=JobResponse)
async def submit_prompt(
    request: PromptRequest,
    background_tasks: BackgroundTasks,
    llm_service=Depends(use_llm)
    ):
    """
    End-to-end endpoint to process a scientific prompt into a complete 3D visualization.
    This endpoint now launches a background task and returns immediately with a job ID.
    
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
        # Check API key early
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        
        # Initial validation can be done synchronously to reject bad prompts immediately
        domain_validator = DomainValidator(llm_service)
        validation_result = domain_validator.is_scientific(request.prompt)
        
        # If not scientific, reject the prompt immediately
        if not validation_result.is_true:
            return {
                "job_id": "rejected",
                "status": "failed",
                "message": "Non-scientific prompt rejected",
                "error": validation_result.reasoning
            }
        
        # Generate a unique job ID
        import uuid
        job_id = str(uuid.uuid4())
        
        # Start the background task
        background_tasks.add_task(
            process_prompt_pipeline_task,
            job_id=job_id,
            prompt=request.prompt,
            llm_service=llm_service
        )
        
        # Create job entry
        pipeline_jobs[job_id] = {
            "status": "processing",
            "prompt": request.prompt,
            "created_at": datetime.now().isoformat(),
            "result": None,
            "progress": 0.0  # Initial progress
        }
        
        print(f"Started background processing job {job_id} for prompt: {request.prompt[:50]}...")
        
        # Return immediately with the job ID and status
        return {
            "job_id": job_id,
            "status": "processing",
            "message": "Processing started in the background. Poll /prompt/process/{job_id} for updates.",
            "progress": 0.0
        }
    except Exception as e:
        print(f"Unhandled exception in submit_prompt: {str(e)}")
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Unhandled error: {str(e)}")

@router.post("/validate-scientific/", response_model=ValidationResponse)
async def validate_scientific(
    request: PromptRequest,
    llm_service=Depends(use_llm)
    ):
    """
    Endpoint to validate if a prompt is scientific in nature.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        domain_validator = DomainValidator(llm_service)

        validation_result = domain_validator.is_scientific(request.prompt)
        
        return {
            "is_scientific": validation_result.is_true,
            "confidence": validation_result.confidence,
            "reasoning": validation_result.reasoning or ""
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/generate-geometry/", response_model=GeometryResponse)
async def generate_geometry(
    request: GeometryRequest,
    llm_service=Depends(use_llm)
    ):
    """
    Endpoint to generate Three.js geometry based on user prompt.
    Only generates geometry for scientific content.
    
    This endpoint processes the request immediately and returns the result.
    For longer processing, use the /prompt/process/ endpoint.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        
        # Validate the prompt is scientific
        domain_validator = DomainValidator(llm_service)
        validation_result = domain_validator.is_scientific(request.prompt)
        
        if not validation_result.is_true:
            raise HTTPException(
                status_code=400,
                detail=f"Non-scientific prompt rejected: {validation_result.reasoning}"
            )
            
        # Generate the geometry directly for immediate response
        # This is what the client expects
        geometry_agent = GeometryAgent(llm_service)
        generated_code = geometry_agent.get_geometry_snippet(request.prompt)
        
        return {"result": generated_code}
    except Exception as e:
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/generate-script/", response_model=SceneScript)
async def generate_script(
    request: PromptRequest,
    llm_service=Depends(use_llm)
    ):
    """
    Endpoint to generate a structured scene script based on user prompt.
    Only generates scripts for scientific content.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        
        # First validate that the prompt is scientific
        validation_result = domain_validator.is_scientific(request.prompt)
        
        # If not scientific, reject the prompt
        if not validation_result.is_true:
            raise HTTPException(
                status_code=400,
                detail=f"Non-scientific prompt rejected: {validation_result.reasoning}"
            )
        script_agent = ScriptAgent(llm_service)
        # If scientific, generate the script
        animation_script = script_agent.generate_script(request.prompt)
        return animation_script
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

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
    llm_service=Depends(use_llm)
    ):
    """
    Endpoint to generate an orchestration plan from a scene script.
    Breaks down the script into discrete objects needed for the visualization.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        orchestration_agent = OrchestrationAgent(llm_service)
        # Generate the orchestration plan from the script
        orchestration_plan = orchestration_agent.generate_orchestration_plan(request.script)
        return orchestration_plan
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# Background task to generate geometry for all objects in the plan
async def generate_geometry_task(job_id: str, plan: OrchestrationPlan):
    """Background task to generate geometry for all objects in the plan"""
    try:
        # Initialize job status
        geometry_jobs[job_id] = {
            "status": "processing",
            "total": len(plan.objects),
            "completed": 0,
            "results": None
        }
        
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
async def generate_geometry_for_plan(request: OrchestrationRequest, background_tasks: BackgroundTasks):
    """
    Endpoint to generate Three.js geometry for all objects in an orchestration plan.
    This starts an asynchronous job that processes each object sequentially.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        
        # Generate a unique job ID
        import uuid
        job_id = str(uuid.uuid4())
        
        # Start the background task to generate geometry
        background_tasks.add_task(generate_geometry_task, job_id, request.plan)
        
        # Return the job ID so the client can check status
        return {
            "job_id": job_id,
            "status": "processing",
            "message": f"Started processing {len(request.plan.objects)} objects"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/geometry-job-status/{job_id}", response_model=GeometryResultResponse)
async def get_geometry_job_status(job_id: str):
    """
    Check the status of a geometry generation job.
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
        "results": job.get("results", None)
    }

@router.post("/generate-animation/", response_model=AnimationResponse)
async def generate_animation(
    request: AnimationRequest,
    llm_service=Depends(use_llm)
    ):
    """
    Generate animation code for a scene based on the script, orchestration plan, and generated geometries.
    This endpoint takes the outputs from previous steps in the pipeline.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        animation_agent = AnimationAgent(llm_service)
        # Generate animation code
        animation = animation_agent.generate_animation_code(
            script=request.script,
            object_geometries=request.object_geometries,
            orchestration_plan=request.orchestration_plan
        )
        
        # Convert keyframes to serializable format
        keyframes_dict = [
            {
                "timecode": keyframe.timecode,
                "actions": keyframe.actions
            }
            for keyframe in animation.keyframes
        ]
        
        return {
            "code": animation.code,
            "keyframes": keyframes_dict
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/package-scene/", response_model=PackagedSceneResponse)
async def package_scene(request: PackagedSceneRequest):
    """
    Create a complete, packaged Three.js scene from the generated components.
    This endpoint combines all outputs from previous steps into a ready-to-use scene.
    """
    try:
        # Use the ScenePackager to create a complete scene
        scene_package = ScenePackager.create_scene_package(
            script=request.script,
            orchestration_plan=request.orchestration_plan,
            object_geometries=request.object_geometries,
            animation_code=request.animation_code
        )
        
        # Save the JS file to the static directory
        import os
        static_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "static")
        os.makedirs(static_dir, exist_ok=True)
        
        # Write the JavaScript to a file
        with open(os.path.join(static_dir, "scene.js"), "w") as f:
            f.write(scene_package.js)
        
        return {
            "html": scene_package.html,
            "js": scene_package.js,
            "minimal_js": scene_package.minimal_js,
            "title": scene_package.title,
            "timecode_markers": scene_package.timecode_markers,
            "total_elements": scene_package.total_elements
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))