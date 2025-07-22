import math
from typing import Any, Dict, List, Optional

import svgwrite  # type: ignore

from ..agent_management.models import DiagramPlan

# CPK colors dictionary (common elements)
CPK_COLORS = {
    "H": "white",
    "C": "black",
    "N": "blue",
    "O": "red",
    "F": "lightgreen",
    "Cl": "green",
    "Br": "darkred",
    "I": "purple",
    "S": "yellow",
    "P": "orange",
    "B": "tan",
    "Si": "grey",
    # Add more elements and colors as needed
    "DEFAULT": "pink",  # Default color for unknown elements
}

TEXT_COLORS = {
    "white": "black",  # Text color for light backgrounds
    "black": "white",  # Text color for dark backgrounds
    "blue": "white",
    "red": "white",
    "lightgreen": "black",
    "green": "white",
    "darkred": "white",
    "purple": "white",
    "yellow": "black",
    "orange": "black",
    "tan": "black",
    "grey": "white",
    "DEFAULT": "black",
}


def _render_single_molecule(dwg: svgwrite.Drawing, data: Dict[str, Any]):
    atoms = data.get("atoms", [])
    bonds = data.get("bonds", [])
    box = data.get("box", {})
    label = data.get("label")
    label_position = data.get("label_position", "below")

    if not atoms:
        return

    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    # Add a small epsilon to prevent division by zero for single atoms or linear molecules with no height/width
    mol_width = (max_x - min_x) or 0.1
    mol_height = (max_y - min_y) or 0.1

    # Use a slightly larger radius for atoms to make them more visible
    atom_radius = 5
    font_size_pixels = atom_radius * 1.5  # Make font size relative to radius

    # Calculate scale and offset to fit molecule within its designated box
    # The box width/height are from the LLM's plan
    box_w = box.get("width", 100.0)  # Default box width if not provided by LLM
    box_h = box.get("height", 100.0)  # Default box height if not provided by LLM

    scale_w = box_w / mol_width
    scale_h = box_h / mol_height
    scale = min(scale_w, scale_h) * 0.8  # Use 80% of box to leave some padding

    # Center the molecule within its box
    offset_x = box.get("x", 0) + (box_w - (mol_width * scale)) / 2 - (min_x * scale)
    offset_y = box.get("y", 0) + (box_h - (mol_height * scale)) / 2 - (min_y * scale)

    for bond in bonds:
        a1 = atoms[bond["start"]]
        a2 = atoms[bond["end"]]
        x1 = offset_x + a1["x"] * scale
        y1 = offset_y + a1["y"] * scale
        x2 = offset_x + a2["x"] * scale
        y2 = offset_y + a2["y"] * scale
        dwg.add(
            dwg.line(start=(x1, y1), end=(x2, y2), stroke="black", stroke_width=1.5)
        )

    for atom_data in atoms:
        element = atom_data.get("element", "X")  # Default to X if no element
        cx = offset_x + atom_data["x"] * scale
        cy = offset_y + atom_data["y"] * scale

        atom_color = CPK_COLORS.get(element.capitalize(), CPK_COLORS["DEFAULT"])
        text_color = TEXT_COLORS.get(atom_color, TEXT_COLORS["DEFAULT"])

        dwg.add(
            dwg.circle(
                center=(cx, cy),
                r=atom_radius,
                fill=atom_color,
                stroke="black",
                stroke_width=0.5,
            )
        )
        dwg.add(
            dwg.text(
                element,
                insert=(cx, cy),
                font_size=f"{font_size_pixels}px",
                fill=text_color,
                text_anchor="middle",
                dominant_baseline="central",  # More reliable for vertical centering
            )
        )

    if label:
        label_font_size = 12
        padding = 5  # Padding between molecule box and label
        if label_position == "above":
            x = box.get("x", 0) + box_w / 2
            y = (
                box.get("y", 0) - padding - (label_font_size / 2)
            )  # Position above the box
            anchor = "middle"
        elif label_position == "below":
            x = box.get("x", 0) + box_w / 2
            y = (
                box.get("y", 0) + box_h + padding + label_font_size
            )  # Position below the box
            anchor = "middle"
        elif label_position == "left":
            x = box.get("x", 0) - padding
            y = box.get("y", 0) + box_h / 2
            anchor = "end"  # Anchor text to its end for left positioning
        elif label_position == "right":
            x = box.get("x", 0) + box_w + padding
            y = box.get("y", 0) + box_h / 2
            anchor = "start"  # Anchor text to its start for right positioning
        else:  # Default to below if invalid position
            x = box.get("x", 0) + box_w / 2
            y = box.get("y", 0) + box_h + padding + label_font_size
            anchor = "middle"

        dwg.add(
            dwg.text(
                label,
                insert=(x, y),
                text_anchor=anchor,
                font_size=f"{label_font_size}px",
            )
        )


def render_diagram(
    plan: DiagramPlan,
    molecule_data: Dict[
        str, Any
    ],  # Can be Dict[str, Any] or Dict[str, List[Dict[str, Any]]]
    canvas_width: Optional[int] = None,
    canvas_height: Optional[int] = None,
) -> str:
    """Render a molecular diagram as SVG.

    Args:
        plan: The diagram plan containing molecule positions and arrows
        molecule_data: Dictionary mapping molecule names to their 2D structure data
        canvas_width: Optional canvas width (default from plan)
        canvas_height: Optional canvas height (default from plan)

    Returns:
        SVG string representation of the diagram
    """
    # Use provided dimensions or fall back to plan defaults
    width = canvas_width or plan.canvas_width
    height = canvas_height or plan.canvas_height

    # Create SVG drawing
    dwg = svgwrite.Drawing(size=(width, height))

    # Add molecules
    for molecule in plan.molecule_list:
        name = molecule.molecule
        if name in molecule_data:
            data_list = molecule_data[name]
            # molecule_data contains a list, but we need a single dict for rendering
            if isinstance(data_list, list) and len(data_list) > 0:
                data = data_list[0].copy()  # Use first item and make a copy
            elif isinstance(data_list, dict):
                data = data_list.copy()  # Handle case where it's already a dict
            else:
                continue  # Skip if no valid data

            data["box"] = {"x": molecule.x, "y": molecule.y}
            _render_single_molecule(dwg, data)

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

            # Arrow head parameters
            arrow_size = 10
            angle = math.pi / 6  # 30 degrees

            # Calculate arrow head points
            p1x = end[0] - arrow_size * (udx * math.cos(angle) - udy * math.sin(angle))
            p1y = end[1] - arrow_size * (udx * math.sin(angle) + udy * math.cos(angle))
            p2x = end[0] - arrow_size * (udx * math.cos(angle) + udy * math.sin(angle))
            p2y = end[1] - arrow_size * (-udx * math.sin(angle) + udy * math.cos(angle))

            # Draw arrow line
            dwg.add(dwg.line(start=start, end=end, stroke="black", stroke_width=2))

            # Draw arrow head
            dwg.add(dwg.polygon(points=[end, (p1x, p1y), (p2x, p2y)], fill="black"))

    return dwg.tostring()
