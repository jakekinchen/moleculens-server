"""Routes for graphic generation using YAML specs."""

import base64
import io
import logging
import tempfile
import traceback
from pathlib import Path
from typing import Any, Literal, Optional

import httpx
import yaml
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from api.diagram.models import DiagramPlan
from api.diagram.planner_llm.service import plan, plan_simple
from api.diagram.renderer.diagram import render_diagram
from api.diagram.validator.core import validate

logger = logging.getLogger(__name__)

router = APIRouter(
    prefix="/graphic",
    tags=["Graphic"],
    responses={404: {"description": "Not found"}},
)


class GraphicPlanRequest(BaseModel):
    """Request model for generating YAML specs from natural language."""

    brief: str = Field(description="Natural language description of the desired diagram")
    context: str = Field(default="Educational molecular diagram", description="Context for the graphic")
    theme: str = Field(default="Clean scientific visualization", description="Visual theme")
    width: int = Field(default=960, description="Canvas width in pixels")
    height: int = Field(default=640, description="Canvas height in pixels")
    sections: Optional[str] = Field(default=None, description="Section breakdown")
    notes: Optional[str] = Field(default=None, description="Additional requirements")
    model_name: str = Field(default="o3-mini", description="LLM model to use")


class GraphicPlanResponse(BaseModel):
    """Response model for YAML spec generation."""

    yaml_spec: str = Field(description="Generated YAML specification")
    status: str = Field(default="completed", description="Generation status")
    error: Optional[str] = Field(default=None, description="Error message if failed")


class GraphicValidateRequest(BaseModel):
    """Request model for validating YAML specs."""

    yaml_spec: str = Field(description="YAML specification to validate")


class GraphicValidateResponse(BaseModel):
    """Response model for validation results."""

    valid: bool = Field(description="Whether the spec is valid")
    errors: list = Field(default_factory=list, description="Validation errors if any")
    status: str = Field(default="completed", description="Validation status")


class GraphicRenderRequest(BaseModel):
    """Request model for rendering graphics from YAML specs."""

    yaml_spec: str = Field(description="YAML specification to render")
    deterministic: bool = Field(default=True, description="Use deterministic rendering")
    output_format: Literal["svg", "png", "both"] = Field(default="svg", description="Output format")


class GraphicRenderResponse(BaseModel):
    """Response model for rendered graphics."""

    svg_content: str = Field(description="Rendered SVG content")
    png_base64: Optional[str] = Field(default=None, description="Base64 encoded PNG if requested")
    status: str = Field(default="completed", description="Rendering status")
    error: Optional[str] = Field(default=None, description="Error message if failed")


class GraphicMakeRequest(BaseModel):
    """Request model for full pipeline processing."""

    brief: str = Field(description="Natural language description")
    width: int = Field(default=960, description="Canvas width")
    height: int = Field(default=640, description="Canvas height")
    model_name: str = Field(default="o3-mini", description="LLM model to use")
    output_format: Literal["svg", "png", "both"] = Field(default="both", description="Output format")


class GraphicMakeResponse(BaseModel):
    """Response model for full pipeline results."""

    yaml_spec: str = Field(description="Generated YAML specification")
    svg_content: str = Field(description="Rendered SVG content")
    png_base64: Optional[str] = Field(default=None, description="Base64 encoded PNG if requested")
    job_id: Optional[str] = Field(default=None, description="Job ID for async processing")
    status: str = Field(default="completed", description="Processing status")
    error: Optional[str] = Field(default=None, description="Error message if failed")


async def _generate_molecule_images(plan: DiagramPlan) -> dict[str, str]:
    """Generate transparent PNG images for each molecule in the diagram plan.

    Args:
        plan: DiagramPlan containing molecules to render

    Returns:
        Dictionary mapping molecule names to base64-encoded PNG data
    """
    import asyncio

    molecule_images = {}
    unique_molecules = list(set(mol.molecule for mol in plan.molecule_list))

    async with httpx.AsyncClient(timeout=60.0) as client:
        for i, molecule_name in enumerate(unique_molecules):
            try:
                # Add delay between requests to avoid rate limiting
                if i > 0:
                    await asyncio.sleep(2)  # 2 second delay between requests

                # Create render request for transparent 2D molecule with enhanced descriptions
                enhanced_descriptions = {
                    "glucose": "D-glucose alpha-D-glucopyranose C6H12O6 sugar molecule",
                    "Glucose": "D-glucose alpha-D-glucopyranose C6H12O6 sugar molecule",
                    "C6H12O6": "D-glucose alpha-D-glucopyranose C6H12O6 sugar molecule",
                    "RuBP": "ribulose-1,5-bisphosphate RuBP C5H12O11P2 molecule",
                    "CO2": "carbon dioxide CO2 molecule",
                    "H2O": "water H2O molecule",
                    "O2": "oxygen O2 molecule",
                    "ATP": "adenosine triphosphate ATP energy molecule",
                    "ADP": "adenosine diphosphate ADP molecule",
                }

                enhanced_description = enhanced_descriptions.get(molecule_name, f"{molecule_name} molecule structure")

                render_request = {
                    "description": f"Show {enhanced_description} with transparent background for diagram overlay",
                    "format": "image",
                    "transparent_background": True,
                    "ray_trace": True,
                    "resolution": [300, 300],  # Smaller size for faster rendering
                    "dpi": 100,
                    "background_color": "white",
                }

                logger.info(f"Rendering molecule {i + 1}/{len(unique_molecules)}: {molecule_name}")

                # Call the render endpoint
                response = await client.post("http://localhost:8000/render", json=render_request)

                if response.status_code == 200:
                    # Check if response is JSON (large file redirect) or direct image
                    content_type = response.headers.get("content-type", "")

                    if "application/json" in content_type:
                        # Large file - get URL and fetch it
                        json_response = response.json()
                        image_url = json_response.get("url")
                        if image_url:
                            image_response = await client.get(f"http://localhost:8000{image_url}")
                            if image_response.status_code == 200:
                                image_data = base64.b64encode(image_response.content).decode("utf-8")
                                molecule_images[molecule_name] = image_data
                    else:
                        # Direct image response
                        image_data = base64.b64encode(response.content).decode("utf-8")
                        molecule_images[molecule_name] = image_data

                    logger.info(f"Successfully generated image for molecule: {molecule_name}")
                elif response.status_code == 429:
                    logger.warning(f"Rate limited for molecule {molecule_name}, waiting longer...")
                    await asyncio.sleep(5)  # Wait 5 seconds and retry once
                    # Retry once
                    retry_response = await client.post("http://localhost:8000/render", json=render_request)
                    if retry_response.status_code == 200:
                        image_data = base64.b64encode(retry_response.content).decode("utf-8")
                        molecule_images[molecule_name] = image_data
                        logger.info(f"Successfully generated image for molecule on retry: {molecule_name}")
                    else:
                        logger.warning(
                            f"Failed to render molecule {molecule_name} even after retry: {retry_response.status_code}"
                        )
                else:
                    logger.warning(f"Failed to render molecule {molecule_name}: {response.status_code}")

            except Exception as e:
                logger.error(f"Error rendering molecule {molecule_name}: {str(e)}")
                continue

    logger.info(f"Generated {len(molecule_images)} molecule images out of {len(unique_molecules)} requested")
    return molecule_images


def _render_diagram_with_molecules(
    plan: DiagramPlan,
    molecule_images: dict[str, str],
    canvas_width: int,
    canvas_height: int,
) -> str:
    """Render diagram SVG with embedded molecular images.

    Args:
        plan: DiagramPlan containing layout information
        molecule_images: Dictionary of molecule names to base64 PNG data
        canvas_width: Canvas width in pixels
        canvas_height: Canvas height in pixels

    Returns:
        SVG string with embedded molecular images and arrows
    """
    import math

    import svgwrite

    # Create SVG drawing
    dwg = svgwrite.Drawing(size=(canvas_width, canvas_height))

    # Add white background
    dwg.add(dwg.rect(insert=(0, 0), size=(canvas_width, canvas_height), fill="white"))

    # Add molecular images or fallback representations
    for molecule in plan.molecule_list:
        molecule_name = molecule.molecule

        if molecule_name in molecule_images:
            # Calculate image size and position
            image_size = min(molecule.width or 100, molecule.height or 100)
            x = molecule.x - image_size / 2
            y = molecule.y - image_size / 2

            # Embed base64 image
            image_data = molecule_images[molecule_name]
            dwg.add(
                dwg.image(
                    href=f"data:image/png;base64,{image_data}",
                    insert=(x, y),
                    size=(image_size, image_size),
                )
            )
        else:
            # Fallback: draw a simple colored circle with molecule name
            radius = min(molecule.width or 80, molecule.height or 80) / 2

            # Choose color based on molecule type
            color_map = {
                "CO2": "#90EE90",  # Light green for CO2
                "RuBP": "#FFB6C1",  # Light pink for RuBP
                "glucose": "#87CEEB",  # Sky blue for glucose
                "RuBisCO": "#DDA0DD",  # Plum for enzyme
            }
            fill_color = color_map.get(molecule_name, "#F0F0F0")  # Light gray default

            # Draw circle
            dwg.add(
                dwg.circle(
                    center=(molecule.x, molecule.y),
                    r=radius,
                    fill=fill_color,
                    stroke="black",
                    stroke_width=2,
                )
            )

            # Add molecule name in circle
            dwg.add(
                dwg.text(
                    molecule_name,
                    insert=(molecule.x, molecule.y),
                    text_anchor="middle",
                    dominant_baseline="middle",
                    font_size="12px",
                    font_family="Arial, sans-serif",
                    fill="black",
                    font_weight="bold",
                )
            )

        # Add molecule label below (ensure it stays within canvas)
        if molecule.label:
            label_y = min(molecule.y + (molecule.height or 100) / 2 + 20, canvas_height - 10)
            dwg.add(
                dwg.text(
                    molecule.label,
                    insert=(molecule.x, label_y),
                    text_anchor="middle",
                    font_size="14px",
                    font_family="Arial, sans-serif",
                    fill="black",
                )
            )

    # Add arrows
    for arrow in plan.arrows:
        start = arrow.start
        end = arrow.end

        # Calculate arrow direction
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = math.sqrt(dx * dx + dy * dy)

        if length > 0:
            # Unit direction vector
            udx = dx / length
            udy = dy / length

            # Shorten arrow to not overlap with molecules
            margin = 60  # Distance from molecule center
            start_adj = (start[0] + udx * margin, start[1] + udy * margin)
            end_adj = (end[0] - udx * margin, end[1] - udy * margin)

            # Arrow head parameters
            arrow_size = 10
            angle = math.pi / 6  # 30 degrees

            # Calculate arrow head points
            p1x = end_adj[0] - arrow_size * (udx * math.cos(angle) - udy * math.sin(angle))
            p1y = end_adj[1] - arrow_size * (udx * math.sin(angle) + udy * math.cos(angle))
            p2x = end_adj[0] - arrow_size * (udx * math.cos(angle) + udy * math.sin(angle))
            p2y = end_adj[1] - arrow_size * (-udx * math.sin(angle) + udy * math.cos(angle))

            # Draw arrow line
            dwg.add(
                dwg.line(
                    start=start_adj,
                    end=end_adj,
                    stroke="black",
                    stroke_width=2,
                    marker_end="url(#arrowhead)",
                )
            )

            # Draw arrow head
            dwg.add(dwg.polygon(points=[end_adj, (p1x, p1y), (p2x, p2y)], fill="black"))

            # Add arrow label if present
            if arrow.text:
                mid_x = (start_adj[0] + end_adj[0]) / 2
                mid_y = (start_adj[1] + end_adj[1]) / 2 - 10  # Slightly above the arrow
                dwg.add(
                    dwg.text(
                        arrow.text,
                        insert=(mid_x, mid_y),
                        text_anchor="middle",
                        font_size="12px",
                        font_family="Arial, sans-serif",
                        fill="black",
                    )
                )

    return dwg.tostring()


def _svg_to_png(svg_content: str, width: int = 960, height: int = 640) -> str:
    """Convert SVG content to PNG and return as base64 string.

    Args:
        svg_content: SVG content as string
        width: Output width in pixels
        height: Output height in pixels

    Returns:
        Base64 encoded PNG data
    """
    try:
        # Try using svglib + reportlab + PIL (most compatible)
        from PIL import Image
        from reportlab.graphics import renderPM
        from svglib.svglib import renderSVG

        # Parse SVG and create a temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg_content)
            svg_path = f.name

        try:
            # Convert SVG to ReportLab drawing
            drawing = renderSVG.svg2rlg(svg_path)

            # Scale the drawing to desired size
            if drawing:
                drawing.width = width
                drawing.height = height
                drawing.scale(
                    width / drawing.width if drawing.width else 1,
                    height / drawing.height if drawing.height else 1,
                )

                # Render to PNG bytes
                png_data = renderPM.drawToString(drawing, fmt="PNG")
                return base64.b64encode(png_data).decode("utf-8")
        finally:
            Path(svg_path).unlink(missing_ok=True)

    except ImportError:
        pass
    except Exception as e:
        logger.warning(f"svglib conversion failed: {e}")

    try:
        # Try using cairosvg (if Cairo is properly installed)
        import cairosvg

        png_data = cairosvg.svg2png(
            bytestring=svg_content.encode("utf-8"),
            output_width=width,
            output_height=height,
        )
        return base64.b64encode(png_data).decode("utf-8")
    except ImportError:
        pass
    except Exception as e:
        logger.warning(f"cairosvg conversion failed: {e}")

    try:
        # Fallback to wand (ImageMagick)
        from wand.color import Color
        from wand.image import Image as WandImage

        with WandImage() as img:
            img.format = "svg"
            img.read(blob=svg_content.encode("utf-8"))
            img.format = "png"
            img.resize(width, height)
            img.background_color = Color("white")
            return base64.b64encode(img.make_blob()).decode("utf-8")
    except ImportError:
        pass
    except Exception as e:
        logger.warning(f"wand conversion failed: {e}")

    try:
        # Fallback: Create a simple placeholder PNG using PIL
        from PIL import Image, ImageDraw, ImageFont

        # Create a white background image
        img = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(img)

        # Add a simple message
        try:
            font = ImageFont.truetype("Arial.ttf", 24)
        except OSError:
            font = ImageFont.load_default()

        text = "SVG Diagram\n(PNG conversion not available)"
        bbox = draw.textbbox((0, 0), text, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]

        x = (width - text_width) // 2
        y = (height - text_height) // 2

        draw.text((x, y), text, fill="black", font=font)

        # Convert to PNG bytes
        buffer = io.BytesIO()
        img.save(buffer, format="PNG")
        png_data = buffer.getvalue()

        return base64.b64encode(png_data).decode("utf-8")

    except ImportError:
        pass
    except Exception as e:
        logger.warning(f"PIL fallback failed: {e}")

    # If all conversion methods fail, return empty string
    logger.warning("No SVG to PNG conversion library available (svglib, cairosvg, wand, PIL)")
    return ""


@router.post("/plan", response_model=GraphicPlanResponse)
async def plan_graphic(request: GraphicPlanRequest) -> GraphicPlanResponse:
    """Generate YAML spec from natural language brief."""
    try:
        if not request.brief.strip():
            raise HTTPException(status_code=422, detail="Brief cannot be empty")

        yaml_spec = plan(
            brief=request.brief,
            context=request.context,
            theme=request.theme,
            width=request.width,
            height=request.height,
            sections=request.sections,
            notes=request.notes,
            model_name=request.model_name,
        )

        return GraphicPlanResponse(yaml_spec=yaml_spec, status="completed")

    except Exception as e:
        logger.error(f"Failed to generate plan: {str(e)}\n{traceback.format_exc()}")
        return GraphicPlanResponse(yaml_spec="", status="failed", error=str(e))


@router.post("/validate", response_model=GraphicValidateResponse)
async def validate_graphic(request: GraphicValidateRequest) -> GraphicValidateResponse:
    """Validate a YAML specification."""
    try:
        if not request.yaml_spec.strip():
            raise HTTPException(status_code=422, detail="YAML spec cannot be empty")

        # Parse YAML first
        try:
            yaml.safe_load(request.yaml_spec)
        except yaml.YAMLError as e:
            return GraphicValidateResponse(valid=False, errors=[f"YAML syntax error: {str(e)}"], status="completed")

        # Validate using our validator
        errors = validate(request.yaml_spec)

        return GraphicValidateResponse(valid=len(errors) == 0, errors=errors, status="completed")

    except Exception as e:
        logger.error(f"Failed to validate spec: {str(e)}\n{traceback.format_exc()}")
        return GraphicValidateResponse(valid=False, errors=[f"Validation error: {str(e)}"], status="failed")


@router.post("/render", response_model=GraphicRenderResponse)
async def render_graphic(request: GraphicRenderRequest) -> GraphicRenderResponse:
    """Render a graphic from YAML specification."""
    try:
        if not request.yaml_spec.strip():
            raise HTTPException(status_code=422, detail="YAML spec cannot be empty")

        # Parse YAML
        try:
            spec_dict = yaml.safe_load(request.yaml_spec)
        except yaml.YAMLError as e:
            raise HTTPException(status_code=422, detail=f"Invalid YAML: {str(e)}") from e

        # Validate first
        errors = validate(request.yaml_spec)
        if errors:
            raise HTTPException(status_code=422, detail=f"Validation errors: {errors}")

        # Convert YAML spec to DiagramPlan for compatibility with existing renderer
        diagram_plan = _yaml_to_diagram_plan(spec_dict)

        # Render using existing diagram renderer
        canvas_width = spec_dict.get("canvas", {}).get("w", 960)
        canvas_height = spec_dict.get("canvas", {}).get("h", 640)

        svg_content = render_diagram(
            plan=diagram_plan,
            molecule_data={},  # Will be populated by the renderer
            canvas_width=canvas_width,
            canvas_height=canvas_height,
        )

        # Convert to PNG if requested
        png_base64 = None
        if request.output_format in ["png", "both"]:
            png_base64 = _svg_to_png(svg_content, canvas_width, canvas_height)

        return GraphicRenderResponse(svg_content=svg_content, png_base64=png_base64, status="completed")

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to render graphic: {str(e)}\n{traceback.format_exc()}")
        return GraphicRenderResponse(svg_content="", status="failed", error=str(e))


@router.post("/make", response_model=GraphicMakeResponse)
async def make_graphic(request: GraphicMakeRequest) -> GraphicMakeResponse:
    """Full pipeline: plan, validate, and render a graphic."""
    try:
        if not request.brief.strip():
            raise HTTPException(status_code=422, detail="Brief cannot be empty")

        # Step 1: Generate YAML spec
        yaml_spec = plan_simple(brief=request.brief, width=request.width, height=request.height)

        # Step 2: Validate the spec
        errors = validate(yaml_spec)
        if errors:
            logger.warning(f"Generated spec has validation errors: {errors}")
            # Continue anyway for now, but log the issues

        # Step 3: Convert to DiagramPlan
        spec_dict = yaml.safe_load(yaml_spec)
        diagram_plan = _yaml_to_diagram_plan(spec_dict)

        # Step 4: Generate molecular images for each molecule
        molecule_images = await _generate_molecule_images(diagram_plan)

        # Step 5: Render the graphic with molecular images
        svg_content = _render_diagram_with_molecules(
            plan=diagram_plan,
            molecule_images=molecule_images,
            canvas_width=request.width,
            canvas_height=request.height,
        )

        # Step 6: Convert to PNG if requested
        png_base64 = None
        if request.output_format in ["png", "both"]:
            png_base64 = _svg_to_png(svg_content, request.width, request.height)

        return GraphicMakeResponse(
            yaml_spec=yaml_spec,
            svg_content=svg_content,
            png_base64=png_base64,
            status="completed",
        )

    except Exception as e:
        logger.error(f"Failed to make graphic: {str(e)}\n{traceback.format_exc()}")
        return GraphicMakeResponse(yaml_spec="", svg_content="", png_base64=None, status="failed", error=str(e))


def _yaml_to_diagram_plan(spec_dict: dict[str, Any]) -> DiagramPlan:
    """Convert YAML spec to DiagramPlan for compatibility with existing renderer."""
    # Extract basic info
    meta = spec_dict.get("meta", {})
    canvas = spec_dict.get("canvas", {})
    cells = spec_dict.get("cells", [])

    # Find DIAGRAM cells and convert to molecule list
    molecule_list = []
    arrows = []

    for cell in cells:
        if cell.get("type") == "DIAGRAM":
            nodes = cell.get("nodes", [])
            edges = cell.get("edges", [])
            # bbox = cell.get("bbox", {})  # Currently unused

            # Convert nodes to molecules with proper spacing
            node_count = len(nodes)
            canvas_width = canvas.get("w", 960)
            canvas_height = canvas.get("h", 640)

            # Leave margin for labels
            margin = 50
            usable_width = canvas_width - 2 * margin
            # usable_height = canvas_height - 2 * margin  # Currently unused

            for i, node in enumerate(nodes):
                # Distribute nodes horizontally across the usable canvas area
                x_pos = margin + (usable_width / (node_count + 1)) * (i + 1)
                y_pos = canvas_height / 2 - 30  # Slightly above center to leave room for labels

                molecule_list.append(
                    {
                        "molecule": node.get("label", node.get("id", "unknown")),
                        "x": x_pos,
                        "y": y_pos,
                        "width": 80,  # Smaller to fit better
                        "height": 80,
                        "label": node.get("label"),
                        "label_position": "below",  # Use valid literal value
                    }
                )

            # Convert edges to arrows with proper coordinates
            for edge in edges:
                # Find source and destination nodes
                src_id = edge.get("src")
                dst_id = edge.get("dst")

                # Find the positions of source and destination nodes
                src_pos = None
                dst_pos = None

                for i, node in enumerate(nodes):
                    if node.get("id") == src_id:
                        src_pos = (
                            (canvas_width / (node_count + 1)) * (i + 1),
                            canvas_height / 2,
                        )
                    if node.get("id") == dst_id:
                        dst_pos = (
                            (canvas_width / (node_count + 1)) * (i + 1),
                            canvas_height / 2,
                        )

                # Create arrow if both positions found
                if src_pos and dst_pos:
                    arrows.append(
                        {
                            "start": src_pos,
                            "end": dst_pos,
                            "text": edge.get("label", ""),
                        }
                    )

    return DiagramPlan(
        plan=meta.get("title", "Generated diagram"),
        molecule_list=molecule_list,
        arrows=arrows,
        canvas_width=canvas.get("w", 960),
        canvas_height=canvas.get("h", 640),
    )


@router.get("/job/{job_id}")
async def get_job_status(job_id: str):
    """Get status of an async job (placeholder for future async processing)."""
    # For now, return a simple response since we're doing synchronous processing
    return {
        "job_id": job_id,
        "status": "completed",
        "message": "Synchronous processing - job completed immediately",
    }
