from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from typing import Optional
from agent_management.agents.geometry_agent import GeometryAgent
from agent_management.agents.domain_bool_agent import DomainValidator
from agent_management.agents.script_agent import ScriptAgent
from agent_management.agents.orchestration_agent import OrchestrationAgent
from agent_management.models import SceneScript, ScriptTimePoint, OrchestrationPlan
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType
import os


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