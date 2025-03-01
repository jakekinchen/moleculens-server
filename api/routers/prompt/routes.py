from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
# from agent_management.agents.geometry_agent import GeometryAgent
# from agent_management.llm_service import LLMService


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
# llm_service = LLMService()
# geometry_agent = GeometryAgent(llm_service)


@router.post("/", response_model=PromptWorking)
def submit_prompt(request: PromptRequest):
    return {
        "prompt": request.prompt,
        "working": True
    }

# @router.post("/generate-geometry/")
# async def generate_geometry(request: GeometryRequest):
#     """
#     Endpoint to generate Three.js geometry based on user prompt.
#     """
#     generated_code = geometry_agent.get_geometry_snippet(request.prompt)
#     return {"javascript_code": generated_code}