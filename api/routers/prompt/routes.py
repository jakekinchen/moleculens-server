from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel

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


@router.post("/", response_model=PromptWorking)
def submit_prompt(request: PromptRequest):
    return {
        "prompt": request.prompt,
        "working": True
    }
