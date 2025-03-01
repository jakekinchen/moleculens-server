from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from agent_management.agents.geometry_agent import GeometryAgent
from agent_management.agents.domain_bool_agent import DomainValidator
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


@router.post("/", response_model=PromptWorking)
def submit_prompt(request: PromptRequest):
    """
    Initial endpoint to validate and submit a prompt for processing.
    Validates if the prompt is scientific before accepting it. This meethod is supposed to kick off a background task to process the prompt through our main pipeline that will return a final output (string) to the user.
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
        
        # If scientific, accept the prompt
        return {
            "prompt": request.prompt,
            "working": True
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