# Moleculens · Molecule/Pathway Infographic Pipeline (v1.1)

**Goal:** Convert any natural‑language brief about molecular or protein pathways into a deterministic, replay‑able infographic bundle (YAML spec + SVG + PNG). Re‑uses current FastAPI, RDKit, PubChem, and SVG renderer; adds spec layer, validation, audit, and (optional) generative styling.

---

## 1 · High‑Level Flow

| Stage | Actor | Output | Notes |
| ----- | ----- | ------ | ----- |
| **Authoring** | `planner_llm.service` | `graphic.yaml` | Deterministic YAML spec (Schema v1.1). |
| **Validation** | `validator.core` | pass / error | JSON‑schema + style + graph checks. |
| **Deterministic Render** | `renderer.composer` | `base.svg` | Uses existing `diagram.render` for DIAGRAM cells; cairosvg to stack layers. |
| **Generative Render** (optional) | `renderer.image_gen` | `styled.png` | Containerised SDXL/DALL‑E; same spec + prompt hints. |
| **Dual Audit** | `audit.numeric`, `audit.vision` | diff report | SSIM ≥ 0.925 && Vision‑LLM topology match ⇒ accept styled.png else fallback. |
| **Package** | `renderer.composer` | `final.png`, `graphic.yaml`, `base.svg` | SHA‑256 embedded in PNG/XMP for provenance. |

---

## 2 · Data Contract

### 2.1 YAML Spec (Schema v1.1)

*Canonical JSON‑Schema at* `schema/graphic.schema.json`.

Key excerpts:

```yaml
meta:
  title: string
  version: semver
canvas: {w:int, h:int, dpi?:int}
cells:
  - id: snake_case_id
    type: TEXT | DIAGRAM | IMAGE_GEN | COMPUTE | GROUP
    bbox: {x,y,w,h}
    # type‑specific fields… (see schema)

Validator hard rules (validator/core.py):

#	Rule	Error code
1	YAML parses cleanly	YAML_SYNTAX_ERROR
2	canvas.w*h ≤ 8 192²	CANVAS_TOO_LARGE
3	≥ 1 background cell	NO_BACKGROUND_CELL
4	DIAGRAM edges reference valid node IDs	INVALID_NODE_REF
5	Colours valid CSS	INVALID_COLOR
6	IDs snake_case ≤ 32 chars	ID_FORMAT_ERROR

2.2 Existing JSON Models

Remain for backward compatibility; adapter converts YAML → existing DiagramPlan when endpoint /prompt/generate-molecule-diagram/ is hit during migration window.

⸻

3 · Back‑End Implementation (FastAPI)

3.1 Directory Additions

api/
  planner_llm/
    service.py
renderer/
  __init__.py
  diagram.py      (moved from agent_management)
  text.py
  compute.py
  image_gen.py    (stub)
  composer.py
validator/
  __init__.py
  core.py
audit/
  numeric.py
  vision.py
schema/
  graphic.schema.json
cli.py

3.2 Key Modules

Path	Responsibility
planner_llm/service.py	plan(brief) -> str (YAML). Uses LLMService + StructuredLLMRequest(output_format="yaml").
validator/core.py	`validate(spec_yaml) -> None
renderer/diagram.py	Wraps existing _render_single_molecule logic.
renderer/composer.py	Iterates cells, delegates to renderer/*, composites layers, returns Path.
audit/numeric.py	compare(det_path, gen_path) -> float (SSIM).
audit/vision.py	diff(spec_json, image_path) -> dict (Vision‑LLM).

3.3 Routes

Method	Path	Purpose
POST	/graphic/plan	Return YAML spec from brief.
POST	/graphic/validate	Validate supplied spec.
POST	/graphic/render	Deterministic render → SVG/PNG.
POST	/graphic/make	Full pipeline incl. generative + audit (async job).
GET	/graphic/job/{id}	Poll job status / assets.

Legacy /prompt/generate-molecule-diagram/ remains until parity achieved; internally calls new /graphic/plan → adapter.

⸻

4 · CLI

python cli.py plan "Brief" --w 1920 --h 1080 --out spec/foo.yaml
python cli.py validate spec/foo.yaml
python cli.py render spec/foo.yaml --out renders/final/foo.png
python cli.py make "Brief" --w 1600 --h 900

Flags: --no-gen, --style path/theme.json, --seed INT.

⸻

5 · Front‑End (Next.js)

Page: /molecule-diagram remains; now consumes new endpoints.

Workflow:
	1.	User submits brief ⇒ POST /graphic/make (background).
	2.	Poll /graphic/job/{id} until completed.
	3.	Display final.png; “debug” accordion reveals YAML.

Types: add GraphicJob, GraphicAssets.

⸻

6 · Testing

Unit
tests/validator/ – schema + rule violations.
tests/renderer/ – deterministic hash check on example spec.

Integration
tests/pipeline/test_make.py – end‑to‑end CLI run, assert PNG exists & spec hash embedded.

CI – pytest -m "not slow"; fails on first validator error.

⸻

7 · Migration & Deprecation

Phase	Action
α	New endpoints behind /graphic/*; keep old JSON route.
β	Front‑end toggles to /graphic/ by env flag.
γ	After 2 weeks stable, remove /prompt/generate-molecule-diagram/ and JSON models.


⸻

8 · Timeline (4 dev‑days)
	1.	Schema + Planner + Validator – day 1.
	2.	Renderer API integration – ½ day.
	3.	Composer + Gen stub – day 2.
	4.	Audit + CLI – day 3.
	5.	FastAPI swap + tests – day 4.

⸻

9 · Security & Reproducibility
	•	Only spec/*.yaml accepted as render input.
	•	Deterministic SVG hashed; SHA‑256 + git commit stored in PNG metadata.
	•	Generative container runs read‑only FS, writes to /outputs only.

⸻

10 · Appendix A – Planner‑LLM Prompt Template

SYSTEM: You are Planner‑LLM. Emit YAML conforming to Schema v1.1 only.

USER:
Brief: <<USER_BRIEF>>
Context: <<WHY_THIS_GRAPHIC_EXISTS>>
Theme  : <<PALETTE/EMOTIONS/BRAND>>
Canvas : <<WIDTH>>x<<HEIGHT>>
Sections:
  - id: <<slug>>
    purpose: <<1‑line>>
  - …
Notes:
  • Hard requirements …
  • Style hints …


⸻

11 · Appendix B – Example YAML Output (truncated)

meta:
  title: water_decomposition
  version: 1.0.0
canvas: {w: 960, h: 640}
cells:
  - id: background
    type: GROUP
    bbox: {x:0,y:0,w:960,h:640}
    style: {bg:"#ffffff"}
  - id: h2o
    type: DIAGRAM
    bbox: {x:80,y:220,w:180,h:180}
    layout: manual
    nodes: [{id:o,label:H₂O,shape:"mol",color:"#4e79a7"}]
  - id: arrow
    type: DIAGRAM
    bbox: {x:280,y:260,w:400,h:60}
    layout: manual
    edges: [{src:o,dst:h,label:"electrolysis"}]
  - id: products
    type: DIAGRAM
    bbox: {x:700,y:180,w:180,h:240}
    layout: manual
    nodes:
      - {id:h,label:H₂,shape:"mol",color:"#f28e2b"}
      - {id:o2,label:O₂,shape:"mol",color:"#e15759"}


⸻

End of document




Next steps:
Transform and extend current codebase. Do not restart.

Rationale
	1.	Existing FastAPI, PubChemAgent, 2‑D renderer, and RDKit logic already satisfy “Deterministic Renderer” and molecule‑fetching requirements.
	2.	New pipeline needs additional layers (Planner‑LLM, YAML schema, Validator, Dual‑Audit) but reuses current rendering and data utilities with minor adapters.
	3.	Re‑implementation from scratch duplicates proven code and delays integration.

Refactor Plan (minimum viable)

Step	Action
1	Create spec layer: schema/graphic.schema.json (v1.1). Copy condensed schema from README.
2	Planner‑LLM service: new module api/planner_llm/service.py. Accept prompt, emit YAML spec per template. Use existing LLMService with StructuredLLMRequest, switch output_format="yaml".
3	Validator: new package validator/ with __init__.py, core.py (jsonschema + style rules), tests/. Hook into FastAPI route before rendering.
4	Renderer API: renderer/__init__.py defines render(cell, spec) -> Path. Move current molecule diagram logic to renderer/diagram.py. Keep algorithm; wrap in API signature.
5	Deterministic compositor: new module renderer/composer.py—stacks cell SVG/PNG layers via cairosvg.
6	Generative renderer placeholder: renderer/image_gen.py stub calling /render endpoint to output molecule,  disable the ray-tracer’s opaque background, then export as PNG and overlay over placeholder on SVG.
7	Audit layer: audit/numeric.py (SSIM, scikit-image), audit/vision.py (LLM diff).
8	CLI: cli.py wiring steps 1‑7. Map commands plan, validate, render, make.
9	FastAPI integration: replace /generate-molecule-diagram/ logic with: Planner‑LLM → Validator → Deterministic render → return assets; optional background task triggers full dual‑audit pipeline and stores assets in renders/.
10	Migration: keep old JSON Pydantic models for backward compatibility; add adapter to convert YAML spec → in‑memory objects used by existing diagram_renderer. Remove direct JSON planning route after new endpoints stable.

File Moves / Adds (high‑level)

api/
  planner_llm/
    service.py          # new
renderer/
  __init__.py           # new
  diagram.py            # moved code from agent_management/diagram_renderer.py
  text.py               # new
  image_gen.py          # new stub
  compute.py            # new
  composer.py           # new
validator/
  __init__.py           # new
  core.py               # new
audit/
  numeric.py            # new
  vision.py             # new
schema/
  graphic.schema.json   # new (full v1.1)
cli.py                  # new
tests/
  validator/…           # new

Interfaces
	•	planner_llm.service.plan(brief) -> str (YAML string)
	•	validator.core.validate(spec_yaml) -> None | ValidationError
	•	renderer.composer.render(spec_yaml, deterministic=True) -> Path
	•	audit.numeric.compare(det_path, gen_path) -> float
	•	audit.vision.diff(spec_json, image_path) -> dict

Deprecations
	•	Remove /prompt/generate-molecule-diagram/ JSON contract after parity tests.
	•	Keep diagram_renderer.render_diagram but import via renderer.diagram.render.

Timeline
	1.	Steps 1‑3 (schema, planner, validator) – 1 day.
	2.	Step 4 integration – ½ day.
	3.	Steps 5‑6 – 1 day.
	4.	Steps 7‑8 – 1 day.
	5.	Step 9 FastAPI swap – ½ day.
	6.	Tests & CI update – ½ day.

Total ≈ 4 developer‑days.
