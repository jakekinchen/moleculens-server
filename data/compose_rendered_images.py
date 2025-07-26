import io
import os

import requests
from PIL import Image, ImageDraw, ImageFont

# API endpoint (adjust if needed)
API_URL = "http://localhost:8000/render"

# Four detailed render requests with transparent backgrounds and annotations
RENDER_REQUESTS = [
    {
        "description": "Caffeine molecule with distance and angle annotations",
        "transparent_background": True,
        "ray_trace": True,
        "dpi": 300,
        "annotations": [
            {
                "type": "distance",
                "name": "bond1",
                "atom1": "resi 1 and name CA",
                "atom2": "resi 2 and name CA",
            },
            {
                "type": "angle",
                "name": "angle1",
                "atom1": "resi 1 and name CA",
                "atom2": "resi 2 and name CA",
                "atom3": "resi 3 and name CA",
            },
            {"type": "label", "selection": "resi 10", "text": "Active Site"},
        ],
        "structure_id": "caffeine",
    },
    {
        "description": "Protein binding site with surface and label",
        "transparent_background": True,
        "ray_trace": True,
        "dpi": 300,
        "annotations": [
            {"type": "label", "selection": "resi 50-60", "text": "Binding Site"},
            {
                "type": "distance",
                "name": "site_dist",
                "atom1": "resi 50 and name CA",
                "atom2": "resi 60 and name CA",
            },
        ],
        "structure_id": "1ubq",
        "highlight_selection": "resi 50-60",
    },
    {
        "description": "Enzyme with publication quality and annotation",
        "transparent_background": True,
        "ray_trace": True,
        "dpi": 400,
        "ray_trace_mode": "poster",
        "annotations": [
            {"type": "label", "selection": "resi 100", "text": "Catalytic Residue"},
            {
                "type": "angle",
                "name": "cat_angle",
                "atom1": "resi 100 and name CA",
                "atom2": "resi 101 and name CA",
                "atom3": "resi 102 and name CA",
            },
        ],
        "structure_id": "3eiy",
    },
    {
        "description": "Annotated molecule with cartoon outline",
        "transparent_background": True,
        "ray_trace": True,
        "dpi": 200,
        "ray_trace_mode": "cartoon_outline",
        "annotations": [
            {"type": "label", "selection": "resi 5", "text": "N-Terminus"},
            {"type": "label", "selection": "resi 120", "text": "C-Terminus"},
        ],
        "structure_id": "1crn",
    },
]

LABELS = [
    "A: Caffeine annotated",
    "B: Binding site",
    "C: Enzyme quality",
    "D: Cartoon outline",
]

OUTPUT_DIR = "rendered_images"
COMPOSITE_PATH = "composition.png"

os.makedirs(OUTPUT_DIR, exist_ok=True)

images = []

for i, req in enumerate(RENDER_REQUESTS):
    print(f"Requesting image {i + 1}...")
    response = requests.post(API_URL, json=req)
    if response.status_code != 200:
        print(f"Failed to render image {i + 1}: {response.text}")
        continue
    img = Image.open(io.BytesIO(response.content)).convert("RGBA")
    img_path = os.path.join(OUTPUT_DIR, f"img_{i + 1}.png")
    img.save(img_path)
    images.append(img)
    print(f"Saved {img_path}")

# Determine grid size (assume all images same size)
if len(images) < 4:
    print("Not enough images to compose a 2x2 grid.")
    exit(1)

w, h = images[0].size
margin = 40
label_height = 40
bg_color = (240, 240, 255, 255)

composite = Image.new("RGBA", (2 * w + 3 * margin, 2 * h + 3 * margin + 2 * label_height), bg_color)


def draw_label(draw, text, x, y, w, h):
    try:
        font = ImageFont.truetype("arial.ttf", 28)
    except OSError:
        font = ImageFont.load_default()
    # Use textbbox if available (Pillow >=8.0), else fallback to textsize
    try:
        bbox = draw.textbbox((0, 0), text, font=font)
        text_w, text_h = bbox[2] - bbox[0], bbox[3] - bbox[1]
    except AttributeError:
        text_w, text_h = draw.textsize(text, font=font)
    draw.rectangle([x, y, x + w, y + label_height], fill=(220, 220, 240, 255))
    draw.text(
        (x + (w - text_w) // 2, y + (label_height - text_h) // 2),
        text,
        fill=(30, 30, 60, 255),
        font=font,
    )


positions = [
    (margin, margin + label_height),
    (w + 2 * margin, margin + label_height),
    (margin, h + 2 * margin + label_height),
    (w + 2 * margin, h + 2 * margin + label_height),
]

draw = ImageDraw.Draw(composite)

for _i, (img, label, pos) in enumerate(zip(images, LABELS, positions)):
    x, y = pos
    composite.paste(img, (x, y), img)
    draw_label(draw, label, x, y - label_height, w, label_height)

composite.save(COMPOSITE_PATH)
print(f"Saved composite image to {COMPOSITE_PATH}")
