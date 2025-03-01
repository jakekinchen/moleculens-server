from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from agent_management.agents.geometry_agent import GeometryAgent
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


# Initialize LLMService and GeometryAgent
llm_config = LLMModelConfig(
    provider=ProviderType.OPENAI,
    model_name="o3-mini",
    api_key=os.getenv("OPENAI_API_KEY")
)
llm_service = LLMService(llm_config)
geometry_agent = GeometryAgent(llm_service)


@router.post("/", response_model=PromptWorking)
def submit_prompt(request: PromptRequest):
    return {
        "prompt": request.prompt,
        "working": True
    }

@router.post("/generate-geometry/", response_model=GeometryResponse)
async def generate_geometry(request: GeometryRequest):
    """
    Endpoint to generate Three.js geometry based on user prompt.
    """
    try:
        if not os.getenv("OPENAI_API_KEY"):
            raise HTTPException(
                status_code=500,
                detail="OPENAI_API_KEY environment variable is not set"
            )
        generated_code = geometry_agent.get_geometry_snippet(request.prompt)
        return {"result": generated_code}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))