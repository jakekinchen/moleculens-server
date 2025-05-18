import svgwrite
from typing import List, Dict, Any


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
    width = max_x - min_x or 1.0
    height = max_y - min_y or 1.0

    scale = min(box.get("width", 1.0) / width, box.get("height", 1.0) / height)
    offset_x = box.get("x", 0) + (box.get("width", 0) - width * scale) / 2 - min_x * scale
    offset_y = box.get("y", 0) + (box.get("height", 0) - height * scale) / 2 - min_y * scale

    for bond in bonds:
        a1 = atoms[bond["start"]]
        a2 = atoms[bond["end"]]
        x1 = offset_x + a1["x"] * scale
        y1 = offset_y + a1["y"] * scale
        x2 = offset_x + a2["x"] * scale
        y2 = offset_y + a2["y"] * scale
        dwg.add(dwg.line(start=(x1, y1), end=(x2, y2), stroke="black"))

    for atom in atoms:
        cx = offset_x + atom["x"] * scale
        cy = offset_y + atom["y"] * scale
        dwg.add(dwg.circle(center=(cx, cy), r=2, fill="white", stroke="black"))
        dwg.add(dwg.text(atom["element"], insert=(cx + 4, cy - 4), font_size="10px"))

    if label:
        if label_position in ("above", "below"):
            x = box.get("x", 0) + box.get("width", 0) / 2
            y = box.get("y", 0) - 5 if label_position == "above" else box.get("y", 0) + box.get("height", 0) + 12
            anchor = "middle"
        else:
            x = box.get("x", 0) - 5 if label_position == "left" else box.get("x", 0) + box.get("width", 0) + 5
            y = box.get("y", 0) + box.get("height", 0) / 2
            anchor = "start"
        dwg.add(dwg.text(label, insert=(x, y), text_anchor=anchor, font_size="12px"))


def render_diagram(molecules: List[Dict[str, Any]], arrows: List[Dict[str, Any]], width: int, height: int) -> str:
    dwg = svgwrite.Drawing(size=(width, height))

    for mol in molecules:
        _render_single_molecule(dwg, mol)

    for arrow in arrows or []:
        start = arrow.get("start")
        end = arrow.get("end")
        if start and end:
            dwg.add(dwg.line(start=start, end=end, stroke="black"))
            if arrow.get("text"):
                mid_x = (start[0] + end[0]) / 2
                mid_y = (start[1] + end[1]) / 2
                dwg.add(dwg.text(arrow["text"], insert=(mid_x, mid_y - 4), text_anchor="middle", font_size="12px"))

    return dwg.tostring()
