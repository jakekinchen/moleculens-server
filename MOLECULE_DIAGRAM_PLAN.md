# Molecule Diagram Feature Implementation Plan

## Overview
Add a new page at https://www.moleculens.com/molecule-diagram that lets users enter a prompt describing a 2-D molecular diagram (single or multiple molecules arranged in a layout). The backend will use an LLM to interpret the prompt, determine molecule placement, and return a rendered diagram.

Example prompt:
> "show a water decomposition reaction for a group of middle schoolers"

## A. Back-end (FastAPI)

### 1. Router
- File: `api/routers/prompt/routes.py`
- New endpoint: `POST /prompt/generate-molecule-diagram/`

### 2. Request/Response Models
```python
class MoleculePlacement(BaseModel):
    molecule: str               # e.g. "H2O"  (can be common name)
    x: float                    # centre-point or top-left, we'll document
    y: float
    width: Optional[float] = None   # optional bounding box
    height: Optional[float] = None
    label: Optional[str] = None     # override label
    label_position: Literal["above","below","left","right"] = "below"

class Arrow(BaseModel):               # optional
    start: Tuple[float, float]
    end: Tuple[float, float]
    style: Literal["straight","curved"] = "straight"
    text: Optional[str] = None        # e.g. "+ energy"

class DiagramPromptRequest(BaseModel):
    prompt: str
    canvas_width: int = 960
    canvas_height: int = 640
    model: Optional[str] = None
    preferred_model_category: Optional[str] = None

class DiagramPlan(BaseModel):
    plan: str
    molecule_list: List[MoleculePlacement]
    arrows: Optional[List[Arrow]] = None

class DiagramResponse(BaseModel):
    diagram_image: str           # base-64 PNG or SVG string
    diagram_plan: DiagramPlan
    status: Literal["completed","failed","processing"] = "completed"
    job_id: Optional[str] = None
    error: Optional[str] = None
```

### 3. Implementation Flow
a. Parse request
b. Compose an LLM **Structured** prompt (via `LLMService` + `StructuredLLMRequest`) that instructs the model to return strictly-valid JSON in the above schema
c. Validate JSON → `DiagramPlan`
d. Build a list of `{query, box}` where box = `{x,y,width,height}` (fallback width/height if None)
e. Fetch 2-D atom/bond data via `PubChemAgent.get_molecules_2d_layout()`
f. Render SVG:
   - For each molecule, use RDKit's `Draw` module (`rdkit.Chem.Draw.MolToImage` or `MolDraw2DSVG`) scaled into the supplied box
   - Compose onto a single SVG canvas (`svgwrite`)
   - Draw arrows & labels
g. Convert SVG to base-64 PNG (optional) or just return the SVG string
h. Return `DiagramResponse`
i. If rendering proves slow, off-load to a background task and return `job_id` + `processing`, re-using the existing job/poll helpers

### 4. Unit Tests
File: `api/tests/test_diagram.py`
- Mock LLM response for deterministic output
- Assert JSON schema validation and SVG generation succeed

## B. Front-end (Next.js)

### 1. New Page
- `src/pages/molecule-diagram.tsx`
- Clone header/footer/layout from `pages/index.tsx`

### 2. Panels

#### a. DiagramInputPanel.tsx
- Textarea for prompt, submit button, loading indicator
- On submit: `generateMoleculeDiagram(request)`
- On success: calls `onDiagramUpdate`

#### b. DiagramViewer.tsx
- Props: `isLoading`, `diagramImage` (SVG string or base-64 png), `diagramPlan` (optional)
- Renders the SVG/PNG; optionally collapsible "Plan (debug)" panel showing the JSON

### 3. API Helper
Extend `src/services/api.ts`:
```typescript
export const generateMoleculeDiagram = async (
   request: { prompt: string; model?: string; preferred_model_category?: string; }
): Promise<DiagramResponse> => { … }          // uses /prompt/generate-molecule-diagram/
```

### 4. Types
- Add `DiagramResponse`, `MoleculePlacement`, `Arrow` TS interfaces in `src/types`

### 5. Routing & Navigation
- Page auto-maps to `/molecule-diagram`
- (Optional) add link in header nav later

## C. LLM Prompt Template (backend)

```
SYSTEM:
You are a chemistry teaching assistant.
Return ONLY valid JSON conforming to the schema below.
Do NOT wrap in markdown.

USER_PROMPT: {user prompt}

JSON_SCHEMA:
{
 "plan": "<plain language reasoning>",
 "molecule_list": [
   { "molecule": "H2O", "x": 320, "y": 360, "width": 200, "height": 200, "label": "2H₂O", "label_position": "below" },
   …
 ],
 "arrows": [
   { "start": [320,360], "end": [640,360], "style": "straight", "text": "electricity" }
 ]
}
```

`LLMService` validates and retries on invalid JSON as done elsewhere.

## D. Milestones / Sequencing

1. Add Pydantic models & endpoint skeleton 🡒 unit tests (LLM mocked)
2. LLM prompt + JSON validation
3. Diagram renderer (SVG placeholder first)
4. Front-end page + API helper + components
5. Wire up, local test
6. Optional refinements: job polling, nicer SVG styling, header nav
