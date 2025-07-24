"""Planner LLM service for generating YAML specs from natural language briefs."""

import logging
from typing import Any, Dict, Optional

from pydantic import BaseModel

from api.utils.llm_direct import generate_structured

logger = logging.getLogger(__name__)

# Planner-LLM Prompt Template from MOLECULE_DIAGRAM_PLAN.md
PLANNER_SYSTEM_PROMPT = """You are Planner-LLM. Emit YAML conforming to Schema v1.1 only.

Generate a molecular diagram specification in YAML format that includes:
- meta: title and version
- canvas: width, height, and optional DPI
- cells: list of diagram elements with id, type, bbox, and type-specific fields

Supported cell types:
- TEXT: Text labels and annotations
- DIAGRAM: Molecular structures and arrows
- IMAGE_GEN: Placeholder for generative content
- COMPUTE: Calculated values or properties
- GROUP: Container for other elements

Each cell must have:
- id: snake_case identifier (≤32 chars)
- type: one of the supported types
- bbox: {x, y, w, h} bounding box coordinates

For DIAGRAM cells, include:
- layout: "manual" or "auto"
- nodes: molecular structures with id, label, shape, color
- edges: arrows/connections with src, dst, label

Use valid CSS colors and ensure all node references are valid."""

PLANNER_USER_TEMPLATE = """Brief: {brief}
Context: {context}
Theme: {theme}
Canvas: {width}x{height}
Sections: {sections}
Notes: {notes}

Generate a complete YAML specification for this molecular diagram."""


class YAMLResponse(BaseModel):
    content: str


def plan(
    brief: str,
    context: str = "Molecular visualization for educational purposes",
    theme: str = "Clean scientific style with blue/white palette",
    width: int = 960,
    height: int = 640,
    sections: Optional[str] = None,
    notes: Optional[str] = None,
    model_name: str = "gpt-4o-mini",
) -> str:
    """Generate YAML spec from natural language brief.

    Args:
        brief: Natural language description of the desired diagram
        context: Why this graphic exists (default: educational)
        theme: Visual style and color preferences
        width: Canvas width in pixels
        height: Canvas height in pixels
        sections: Optional section breakdown
        notes: Additional requirements or constraints
        model_name: LLM model to use for planning

    Returns:
        YAML string conforming to Schema v1.1

    Raises:
        ValueError: If planning fails or invalid parameters provided
    """
    if not brief.strip():
        raise ValueError("Brief cannot be empty")

    if width <= 0 or height <= 0:
        raise ValueError("Canvas dimensions must be positive")

    if width * height > 8192 * 8192:
        raise ValueError("Canvas too large (max 8192x8192)")

    # Prepare sections and notes
    sections_text = (
        sections or "- id: main_diagram\n  purpose: Primary molecular visualization"
    )
    notes_text = (
        notes
        or "• Use clear molecular representations\n• Ensure proper spacing and alignment"
    )

    # Format the user prompt
    user_prompt = PLANNER_USER_TEMPLATE.format(
        brief=brief,
        context=context,
        theme=theme,
        width=width,
        height=height,
        sections=sections_text,
        notes=notes_text,
    )

    try:
        # Generate the plan using direct OpenAI provider
        response = generate_structured(
            user_prompt=user_prompt,
            response_model=YAMLResponse,
            system_prompt=PLANNER_SYSTEM_PROMPT,
            model_name=model_name,
            max_tokens=2000,
            temperature=(
                None if model_name.startswith("o") else 0.3
            ),  # o* models don't support temperature
        )

        # Extract YAML content
        yaml_content = response.content

        if not yaml_content.strip():
            raise ValueError("Empty YAML response from LLM")

        logger.info(f"Generated YAML spec for brief: {brief[:50]}...")
        return yaml_content.strip()

    except Exception as e:
        logger.error(f"Failed to generate plan for brief '{brief[:50]}...': {str(e)}")
        raise ValueError(f"Planning failed: {str(e)}") from e


def plan_simple(brief: str, width: int = 960, height: int = 640) -> str:
    """Simplified planning interface with minimal parameters.

    Args:
        brief: Natural language description
        width: Canvas width (default 960)
        height: Canvas height (default 640)

    Returns:
        YAML specification string
    """
    return plan(
        brief=brief,
        width=width,
        height=height,
        context="Educational molecular diagram",
        theme="Clean scientific visualization",
    )
