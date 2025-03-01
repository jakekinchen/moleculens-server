from fastapi import APIRouter, HTTPException, Depends, BackgroundTasks
from pydantic import BaseModel
from typing import Optional, Dict, Any, List
from agent_management.agents.geometry_agent import GeometryAgent
from agent_management.agents.domain_bool_agent import DomainValidator
from agent_management.agents.script_agent import ScriptAgent
from agent_management.agents.orchestration_agent import OrchestrationAgent
from agent_management.agents.animation_agent import AnimationAgent
from agent_management.models import SceneScript, ScriptTimePoint, OrchestrationPlan, AnimationCode
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
import os
import asyncio


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

class GeometryResponse(BaseModel):
    result: str

class GeometryRequest(BaseModel):
    prompt: str

class ValidationResponse(BaseModel):
    is_scientific: bool
    confidence: float
    reasoning: str


# Initialize LLMService and agents
llm_config = LLMModelConfig(
    provider=ProviderType.OPENAI,
    model_name="o3-mini",
    api_key=os.getenv("OPENAI_API_KEY")
)
llm_service = LLMService(llm_config)
geometry_agent = GeometryAgent(llm_service)
domain_validator = DomainValidator(llm_service)
script_agent = ScriptAgent(llm_service)
orchestration_agent = OrchestrationAgent(llm_service)
animation_agent = AnimationAgent(llm_service)

# Simple in-memory job store for geometry generation jobs
# In a production app, this would use Redis or another persistent store
geometry_jobs = {}


@router.post("/", response_model=PromptWorking)
async def submit_prompt(request: PromptRequest):
    """
    Initial endpoint to validate and submit a prompt for processing.
    Validates if the prompt is scientific before accepting it, then generates:
    1. An animation script with timecodes, descriptions, and captions
    2. An orchestration plan with discrete objects needed for the animation
    
    This method kicks off the main pipeline for processing scientific visualizations.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        
        # Validate that the prompt is scientific
        validation_result = domain_validator.is_scientific(request.prompt)
        
        # If not scientific, reject the prompt
        if not validation_result.is_true:
            raise HTTPException(
                status_code=400,
                detail=f"Non-scientific prompt rejected: {validation_result.reasoning}"
            )
        
        # If scientific, generate an animation script
        animation_script = script_agent.generate_script(request.prompt)
        
        # Generate an orchestration plan based on the script
        orchestration_plan = orchestration_agent.generate_orchestration_plan(animation_script)
        
        # Return the prompt, script, and orchestration plan
        return {
            "prompt": request.prompt,
            "working": True,
            "script": animation_script,
            "orchestration_plan": orchestration_plan
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/validate-scientific/", response_model=ValidationResponse)
async def validate_scientific(request: PromptRequest):
    """
    Endpoint to validate if a prompt is scientific in nature.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        
        validation_result = domain_validator.is_scientific(request.prompt)
        
        return {
            "is_scientific": validation_result.is_true,
            "confidence": validation_result.confidence,
            "reasoning": validation_result.reasoning or ""
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/generate-geometry/", response_model=GeometryResponse)
async def generate_geometry(request: GeometryRequest):
    """
    Endpoint to generate Three.js geometry based on user prompt.
    Only generates geometry for scientific content.
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
        
        # If scientific, generate the geometry
        generated_code = geometry_agent.get_geometry_snippet(request.prompt)
        return {"result": generated_code}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/generate-script/", response_model=SceneScript)
async def generate_script(request: PromptRequest):
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

@router.post("/generate-orchestration/", response_model=OrchestrationPlan)
async def generate_orchestration(request: ScriptRequest):
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
async def generate_animation(request: AnimationRequest):
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