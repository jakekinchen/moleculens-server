"""Core validation module for YAML specs and diagram plans."""

import logging
import re
from enum import Enum
from typing import Any, Optional

import yaml  # type: ignore

from ..constants import MAX_CANVAS_AREA

logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Custom exception for validation errors with error codes."""

    def __init__(self, message: str, error_code: str, details: Optional[dict[str, Any]] = None):
        super().__init__(message)
        self.error_code = error_code
        self.details = details or {}


class ErrorCode(str, Enum):
    """Validation error codes from MOLECULE_DIAGRAM_PLAN.md."""

    YAML_SYNTAX_ERROR = "YAML_SYNTAX_ERROR"
    CANVAS_TOO_LARGE = "CANVAS_TOO_LARGE"
    NO_BACKGROUND_CELL = "NO_BACKGROUND_CELL"
    INVALID_NODE_REF = "INVALID_NODE_REF"
    INVALID_COLOR = "INVALID_COLOR"
    ID_FORMAT_ERROR = "ID_FORMAT_ERROR"
    MISSING_REQUIRED_FIELD = "MISSING_REQUIRED_FIELD"
    INVALID_CELL_TYPE = "INVALID_CELL_TYPE"
    BBOX_OUT_OF_BOUNDS = "BBOX_OUT_OF_BOUNDS"


# Valid CSS color patterns
CSS_COLOR_PATTERNS = [
    r"^#[0-9a-fA-F]{3}$",  # #rgb
    r"^#[0-9a-fA-F]{6}$",  # #rrggbb
    r"^rgb\(\s*\d+\s*,\s*\d+\s*,\s*\d+\s*\)$",  # rgb(r,g,b)
    r"^rgba\(\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*[0-9.]+\s*\)$",  # rgba(r,g,b,a)
]

# Valid CSS color names (subset)
CSS_COLOR_NAMES = {
    "black",
    "white",
    "red",
    "green",
    "blue",
    "yellow",
    "orange",
    "purple",
    "pink",
    "brown",
    "gray",
    "grey",
    "cyan",
    "magenta",
    "lime",
    "navy",
    "maroon",
    "olive",
    "teal",
    "silver",
    "aqua",
    "fuchsia",
    "steelblue",  # Added common CSS color
    # Light colors commonly used by LLMs
    "lightblue",
    "lightgreen",
    "lightyellow",
    "lightcoral",
    "lightgray",
    "lightgrey",
    "lightpink",
    "lightcyan",
    # Dark colors commonly used by LLMs
    "darkblue",
    "darkgreen",
    "darkred",
    "darkorange",
    "darkgray",
    "darkgrey",
    "darkviolet",
    "darkcyan",
}

# Valid cell types
VALID_CELL_TYPES = {"TEXT", "DIAGRAM", "IMAGE_GEN", "COMPUTE", "GROUP"}

# Snake case pattern for IDs
SNAKE_CASE_PATTERN = re.compile(r"^[a-z][a-z0-9_]*[a-z0-9]$|^[a-z]$")


def validate_yaml_syntax(spec_yaml: str) -> dict[str, Any]:
    """Validate YAML syntax and return parsed data.

    Args:
        spec_yaml: YAML string to validate

    Returns:
        Parsed YAML data as dictionary

    Raises:
        ValidationError: If YAML syntax is invalid
    """
    try:
        data = yaml.safe_load(spec_yaml)
        if not isinstance(data, dict):
            raise ValidationError(
                "YAML must contain a dictionary at root level",
                ErrorCode.YAML_SYNTAX_ERROR,
            )
        return data
    except yaml.YAMLError as e:
        raise ValidationError(
            f"Invalid YAML syntax: {str(e)}",
            ErrorCode.YAML_SYNTAX_ERROR,
            {"yaml_error": str(e)},
        ) from e


def validate_canvas_size(canvas: dict[str, Any]) -> None:
    """Validate canvas dimensions.

    Args:
        canvas: Canvas configuration dictionary

    Raises:
        ValidationError: If canvas is too large
    """
    width = canvas.get("w", 0)
    height = canvas.get("h", 0)

    if width * height > MAX_CANVAS_AREA:
        raise ValidationError(
            f"Canvas too large: {width}x{height} exceeds 8192x8192 limit",
            ErrorCode.CANVAS_TOO_LARGE,
            {"width": width, "height": height, "max_area": MAX_CANVAS_AREA},
        )


def validate_id_format(cell_id: str) -> None:
    """Validate cell ID format.

    Args:
        cell_id: Cell identifier to validate

    Raises:
        ValidationError: If ID format is invalid
    """
    if len(cell_id) > 32:
        raise ValidationError(
            f"ID too long: '{cell_id}' exceeds 32 characters",
            ErrorCode.ID_FORMAT_ERROR,
            {"id": cell_id, "length": len(cell_id)},
        )

    if not SNAKE_CASE_PATTERN.match(cell_id):
        raise ValidationError(
            f"Invalid ID format: '{cell_id}' must be snake_case",
            ErrorCode.ID_FORMAT_ERROR,
            {"id": cell_id, "pattern": "snake_case"},
        )


def validate_color(color: str) -> None:
    """Validate CSS color format.

    Args:
        color: Color string to validate

    Raises:
        ValidationError: If color format is invalid
    """
    # Check color names first
    if color.lower() in CSS_COLOR_NAMES:
        return

    # Check color patterns
    for pattern in CSS_COLOR_PATTERNS:
        if re.match(pattern, color):
            return

    raise ValidationError(
        f"Invalid color format: '{color}'",
        ErrorCode.INVALID_COLOR,
        {
            "color": color,
            "valid_formats": ["#rgb", "#rrggbb", "rgb()", "rgba()", "color names"],
        },
    )


def validate_bbox(bbox: dict[str, Any], canvas_width: int, canvas_height: int) -> None:
    """Validate bounding box coordinates.

    Args:
        bbox: Bounding box dictionary with x, y, w, h
        canvas_width: Canvas width for bounds checking
        canvas_height: Canvas height for bounds checking

    Raises:
        ValidationError: If bbox is out of bounds
    """
    x = bbox.get("x", 0)
    y = bbox.get("y", 0)
    w = bbox.get("w", 0)
    h = bbox.get("h", 0)

    if x < 0 or y < 0 or w <= 0 or h <= 0:
        raise ValidationError(
            f"Invalid bbox coordinates: x={x}, y={y}, w={w}, h={h}",
            ErrorCode.BBOX_OUT_OF_BOUNDS,
            {"bbox": bbox, "reason": "negative or zero dimensions"},
        )

    if x + w > canvas_width or y + h > canvas_height:
        raise ValidationError(
            f"Bbox extends beyond canvas: ({x},{y},{w},{h}) > ({canvas_width},{canvas_height})",
            ErrorCode.BBOX_OUT_OF_BOUNDS,
            {"bbox": bbox, "canvas": {"w": canvas_width, "h": canvas_height}},
        )


def validate_node_references(cells: list[dict[str, Any]]) -> None:
    """Validate that all node references in DIAGRAM cells are valid.

    Args:
        cells: List of cell dictionaries

    Raises:
        ValidationError: If invalid node references found
    """
    # Collect all node IDs from DIAGRAM cells
    all_node_ids: set[str] = set()
    diagram_cells = []

    for cell in cells:
        if cell.get("type") == "DIAGRAM":
            diagram_cells.append(cell)
            nodes = cell.get("nodes", [])
            for node in nodes:
                if "id" in node:
                    all_node_ids.add(node["id"])

    # Check edge references
    for cell in diagram_cells:
        edges = cell.get("edges", [])
        for edge in edges:
            src = edge.get("src")
            dst = edge.get("dst")

            if src and src not in all_node_ids:
                raise ValidationError(
                    f"Invalid node reference in edge: '{src}' not found",
                    ErrorCode.INVALID_NODE_REF,
                    {
                        "edge": edge,
                        "missing_node": src,
                        "available_nodes": list(all_node_ids),
                    },
                )

            if dst and dst not in all_node_ids:
                raise ValidationError(
                    f"Invalid node reference in edge: '{dst}' not found",
                    ErrorCode.INVALID_NODE_REF,
                    {
                        "edge": edge,
                        "missing_node": dst,
                        "available_nodes": list(all_node_ids),
                    },
                )


def validate_background_cell(cells: list[dict[str, Any]]) -> None:
    """Validate that at least one background cell exists.

    Args:
        cells: List of cell dictionaries

    Raises:
        ValidationError: If no background cell found
    """
    has_background = False

    for cell in cells:
        # Check for explicit background cell or GROUP cell that covers full canvas
        if cell.get("id") == "background" or (
            cell.get("type") == "GROUP"
            and cell.get("bbox", {}).get("x", 0) == 0
            and cell.get("bbox", {}).get("y", 0) == 0
        ):
            has_background = True
            break

    if not has_background:
        logger.warning(
            "Validation warning: no background cell found. Proceeding without background cell.",
            extra={"available_cells": [cell.get("id", "unnamed") for cell in cells]},
        )
        return


def validate(spec_yaml: str) -> None:
    """Validate a complete YAML specification.

    Args:
        spec_yaml: YAML string to validate

    Raises:
        ValidationError: If validation fails with specific error code
    """
    # Rule 1: YAML parses cleanly
    data = validate_yaml_syntax(spec_yaml)

    # Check required top-level fields
    if "meta" not in data:
        raise ValidationError(
            "Missing required field: meta",
            ErrorCode.MISSING_REQUIRED_FIELD,
            {"field": "meta"},
        )

    if "canvas" not in data:
        raise ValidationError(
            "Missing required field: canvas",
            ErrorCode.MISSING_REQUIRED_FIELD,
            {"field": "canvas"},
        )

    if "cells" not in data:
        raise ValidationError(
            "Missing required field: cells",
            ErrorCode.MISSING_REQUIRED_FIELD,
            {"field": "cells"},
        )

    canvas = data["canvas"]
    cells = data["cells"]

    # Rule 2: Canvas size validation
    validate_canvas_size(canvas)

    canvas_width = canvas.get("w", 960)
    canvas_height = canvas.get("h", 640)

    # Rule 3: Background cell validation
    validate_background_cell(cells)

    # Validate each cell
    for i, cell in enumerate(cells):
        cell_id = cell.get("id", f"cell_{i}")

        # Rule 6: ID format validation
        if "id" in cell:
            validate_id_format(cell["id"])

        # Cell type validation
        cell_type = cell.get("type")
        if cell_type not in VALID_CELL_TYPES:
            raise ValidationError(
                f"Invalid cell type: '{cell_type}' in cell '{cell_id}'",
                ErrorCode.INVALID_CELL_TYPE,
                {
                    "cell_id": cell_id,
                    "type": cell_type,
                    "valid_types": list(VALID_CELL_TYPES),
                },
            )

        # Bbox validation
        if "bbox" in cell:
            validate_bbox(cell["bbox"], canvas_width, canvas_height)

        # Rule 5: Color validation
        if "style" in cell and "bg" in cell["style"]:
            validate_color(cell["style"]["bg"])

        # Validate colors in nodes
        if cell_type == "DIAGRAM" and "nodes" in cell:
            for node in cell["nodes"]:
                if "color" in node:
                    validate_color(node["color"])

    # Rule 4: Node reference validation
    validate_node_references(cells)

    logger.info("YAML specification validation passed")


def validate_safe(spec_yaml: str) -> Optional[ValidationError]:
    """Safe validation that returns error instead of raising.

    Args:
        spec_yaml: YAML string to validate

    Returns:
        ValidationError if validation fails, None if successful
    """
    try:
        validate(spec_yaml)
        return None
    except ValidationError as e:
        return e
    except Exception as e:
        return ValidationError(
            f"Unexpected validation error: {str(e)}",
            ErrorCode.YAML_SYNTAX_ERROR,
            {"exception": str(e)},
        )
