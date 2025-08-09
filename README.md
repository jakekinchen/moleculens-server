# MoleculeLens Server API

A FastAPI-based server for molecular visualization, diagram generation, and protein structure analysis. This server provides comprehensive endpoints for working with molecular data, generating scientific diagrams, and rendering 3D molecular structures.

## Overview

The MoleculeLens Server provides four main route groups:

- **`/graphic`** - Generate scientific diagrams from natural language descriptions
- **`/prompt`** - Process molecular prompts and generate visualizations
- **`/rcsb`** - Fetch protein structures and metadata from RCSB PDB
- **`/render`** - Render 3D molecular structures using PyMOL
- **`/v1/figure`** - Deterministic, content-addressed figure generation (FigureSpec v1)

## Getting Started

### Prerequisites

- Python 3.8+
- PyMOL (for 3D rendering)
- OpenAI API key (for LLM-powered features)

### Installation

```bash
pip install -r api/requirements.txt
```

### Running the Server

```bash
cd api
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

The API will be available at `http://localhost:8000` with interactive documentation at `http://localhost:8000/docs`.

## 🎯 Client Integration: 2D & 3D Data

### Quick Summary for Frontend Developers

| Data Type | Endpoint | Format | Best For |
|-----------|----------|--------|----------|
| **2D Transparent PNG** | `/render` | `format: "image"` | UI overlays, diagrams, thumbnails |
| **3D PDB Data** | `/render` | `format: "model"` | Three.js reconstruction |
| **3D PDB (Known)** | `/rcsb/fetch-structure/` | `format: "pdb"` | Known protein structures |

### Getting 2D Transparent PNGs

Perfect for UI elements, molecule thumbnails, and diagram overlays:

```javascript
const get2DTransparentPNG = async (moleculeName) => {
  const response = await fetch('/render', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      description: `Show ${moleculeName} molecule with transparent background`,
      format: "image",
      transparent_background: true,
      ray_trace: true,
      resolution: [512, 512],
      dpi: 150
    })
  });

  // Handle both direct PNG and large file URL responses
  if (response.headers.get('content-type')?.includes('image/png')) {
    return response.blob(); // Direct PNG blob
  } else {
    const data = await response.json();
    const imageResponse = await fetch(data.url);
    return imageResponse.blob(); // PNG from URL
  }
};

// Usage
const pngBlob = await get2DTransparentPNG('caffeine');
const imageUrl = URL.createObjectURL(pngBlob);
document.getElementById('molecule-img').src = imageUrl;
```

### Getting 3D PDB Data for Three.js

PDB format is recommended over GLB because it's:
- ✅ Lightweight text format
- ✅ Natively supported by Three.js PDBLoader
- ✅ Contains all atomic coordinates and bonds
- ✅ Easier to manipulate and style

```javascript
const get3DPDBData = async (moleculeName) => {
  const response = await fetch('/render', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      description: `Load ${moleculeName} molecule structure`,
      format: "model" // Returns PDB format
    })
  });

  if (response.headers.get('content-type')?.includes('chemical/x-pdb')) {
    return response.text(); // Direct PDB text
  } else {
    const data = await response.json();
    const pdbResponse = await fetch(data.url);
    return pdbResponse.text(); // PDB from URL
  }
};

// Usage with Three.js PDB loader
import { PDBLoader } from 'three/examples/jsm/loaders/PDBLoader.js';

const load3DMolecule = async (moleculeName) => {
  const pdbData = await get3DPDBData(moleculeName);

  const loader = new PDBLoader();
  const pdb = loader.parse(pdbData);

  // Add to scene
  scene.add(pdb);
  return pdb;
};
```

### Enhanced Endpoints (Recommended)

Get both 2D and 3D data in a single request:

```javascript
const getMoleculeData = async (moleculeName) => {
  const response = await fetch('/render/molecule', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      molecule_name: moleculeName,
      render_type: "both", // "2d_transparent", "3d_pdb", or "both"
      size: "medium", // "small", "medium", "large"
      quality: "high" // "fast", "high", "publication"
    })
  });

  const data = await response.json();
  return {
    name: data.molecule_name,
    png: data.png_base64 ? `data:image/png;base64,${data.png_base64}` : data.png_url,
    pdb: data.pdb_data,
    metadata: data.metadata
  };
};

// Usage
const { png, pdb } = await getMoleculeData('caffeine');
document.getElementById('molecule-img').src = png;
// Load PDB into Three.js...
```

### Complete Three.js Integration Example

```javascript
import * as THREE from 'three';
import { PDBLoader } from 'three/examples/jsm/loaders/PDBLoader.js';

class MoleculeViewer {
  constructor(containerId) {
    this.container = document.getElementById(containerId);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    this.renderer = new THREE.WebGLRenderer({ alpha: true });

    this.renderer.setSize(window.innerWidth, window.innerHeight);
    this.container.appendChild(this.renderer.domElement);

    this.loader = new PDBLoader();
  }

  async loadMolecule(moleculeName) {
    try {
      // Get PDB data from API
      const pdbData = await get3DPDBData(moleculeName);

      // Parse with Three.js
      const pdb = this.loader.parse(pdbData);

      // Style the molecule
      pdb.traverse((child) => {
        if (child.isMesh) {
          child.material.transparent = true;
          child.material.opacity = 0.8;
        }
      });

      // Add to scene
      this.scene.add(pdb);

      // Position camera
      const box = new THREE.Box3().setFromObject(pdb);
      const center = box.getCenter(new THREE.Vector3());
      const size = box.getSize(new THREE.Vector3());

      this.camera.position.copy(center);
      this.camera.position.z += size.length();
      this.camera.lookAt(center);

      return pdb;
    } catch (error) {
      console.error('Failed to load molecule:', error);
    }
  }

  render() {
    this.renderer.render(this.scene, this.camera);
  }
}

// Usage
const viewer = new MoleculeViewer('molecule-container');
await viewer.loadMolecule('caffeine');
viewer.render();
```

## API Routes

### 5. FigureSpec v1 (`/v1/figure`)

Deterministic, content-addressed figure generation API. Produces stable `spec_id` for the same spec using canonical JSON and SHA-256, and serves immutable assets at content-addressed URLs.

#### Endpoints
- `POST /v1/figure`
  - Body: FigureSpec v1 JSON
  - Behavior:
    - Normalize deterministically (see Canonicalization) and compute `spec_id = sha256(canonical_json)`
    - If known, return current `{ spec_id, status, urls? }`
    - Else, enqueue work and return quickly with `{ spec_id, status: "queued" }`
  - 200 OK response:
    ```json
    { "spec_id": "<64-hex>", "status": "queued|processing|completed|failed", "urls": { "svg2d": "...", "png2d": "...", "png3d": "...", "glb": "...", "meta": "..." } }
    ```
- `GET /v1/figure/{spec_id}`
  - Returns current `{ spec_id, status, urls? }` or 404 if unknown

Status values: `queued | processing | completed | failed`

#### FigureSpec v1 schema (shape)
- Required top-level keys: `version`, `input`, `render`, `style_preset`, `annotations`, `"3d"` (the key must literally be `"3d"`)
- `version`: `1`
- `input`: `{ kind: "smiles"|"pdb"|"name", value: string, protonation_pH: number, conformer_method: "etkdg"|"none" }`
- `render`: `{ modes: ["2d"|"3d"+], outputs: ["svg"|"png"+], width: int>0, height: int>0, transparent: boolean, dpi: int>0 }`
- `style_preset`: string (e.g., `"nature-2025"`)
- `annotations`: `{ functional_groups: boolean, charge_labels: "none"|"minimal"|"all", atom_numbering: boolean, scale_bar: boolean, legend: "none"|"auto" }`
- `3d`: `{ representation: "cartoon+licorice"|"surface"|"licorice", bg: "transparent"|"black"|"white", camera: { target: "auto", distance: "auto"|number>0, azimuth: number, elevation: number }, lighting: "three_point_soft", quality: "raytrace_high" }`

#### Deterministic spec_id (must match `@moleculens/chem`)
Canonicalization rules:
- Round all finite numbers to 6 decimal places
- Recursively sort all object keys lexicographically
- Preserve array order (but canonicalize elements)
- Serialize JSON with no whitespace: separators (",", ":")

Reference Python (used verbatim):
```python
import json, hashlib

def _round(x):
    if isinstance(x, float): return round(x, 6)
    if isinstance(x, list): return [_round(i) for i in x]
    if isinstance(x, dict): return {k: _round(x[k]) for k in sorted(x)}
    return x

def canonical_json(obj: dict) -> str:
    return json.dumps(_round(obj), separators=(",", ":"), ensure_ascii=False)

def spec_id(spec: dict) -> str:
    return hashlib.sha256(canonical_json(spec).encode('utf-8')).hexdigest()
```

Important: The JSON key must be `"3d"` (not `_3d`). If any model layer uses an alias, it is re-emitted as `"3d"` before hashing.

#### Assets and URLs (content-addressed)
- Stored under:
  - `figures/{spec_id}/2d.svg`
  - `figures/{spec_id}/2d.png`
  - `figures/{spec_id}/3d.png`
  - `figures/{spec_id}/scene.glb`
  - `figures/{spec_id}/meta.json`
- Response `urls` keys: `svg2d`, `png2d`, `png3d`, `glb`, `meta` (absolute URLs)

#### Operational notes
- Idempotent: Posting the same spec returns the same `spec_id` and existing status/urls
- Responds in <2s; heavy work is done asynchronously
- Write-once assets: once produced, not mutated for the same `spec_id`
- CORS: requests originate from our Next.js server via proxy

#### Example
POST request
```json
{
  "version": 1,
  "input": { "kind": "smiles", "value": "CC(=O)Oc1ccccc1C(=O)O", "protonation_pH": 7.4, "conformer_method": "none" },
  "render": { "modes": ["2d","3d"], "outputs": ["svg","png"], "width": 1024, "height": 768, "transparent": true, "dpi": 300 },
  "style_preset": "nature-2025",
  "annotations": { "functional_groups": true, "charge_labels": "minimal", "atom_numbering": false, "scale_bar": true, "legend": "auto" },
  "3d": { "representation": "cartoon+licorice", "bg": "transparent", "camera": { "target": "auto", "distance": "auto", "azimuth": 30, "elevation": 15 }, "lighting": "three_point_soft", "quality": "raytrace_high" }
}
```

Possible responses
```json
{ "spec_id": "<64-hex>", "status": "queued" }
```
```json
{ "spec_id": "<64-hex>", "status": "completed", "urls": { "svg2d": "https://cdn/.../2d.svg", "png2d": "https://cdn/.../2d.png", "png3d": "https://cdn/.../3d.png", "glb": "https://cdn/.../scene.glb", "meta": "https://cdn/.../meta.json" } }
```

### 1. Graphic Routes (`/graphic`)

Generate scientific diagrams from natural language descriptions using YAML specifications.

#### `POST /graphic/plan`

Generate a YAML specification from a natural language description.

**Request Schema:**
```json
{
  "brief": "string (required) - Natural language description of the desired diagram",
  "context": "string (default: 'Educational molecular diagram') - Context for the graphic",
  "theme": "string (default: 'Clean scientific visualization') - Visual theme",
  "width": "integer (default: 960) - Canvas width in pixels",
  "height": "integer (default: 640) - Canvas height in pixels",
  "sections": "string (optional) - Section breakdown",
  "notes": "string (optional) - Additional requirements",
  "model_name": "string (default: 'o3-mini') - LLM model to use"
}
```

**Response Schema:**
```json
{
  "yaml_spec": "string - Generated YAML specification",
  "status": "string - Generation status (completed/failed)",
  "error": "string (optional) - Error message if failed"
}
```

#### `POST /graphic/validate`

Validate a YAML specification against the schema.

**Request Schema:**
```json
{
  "yaml_spec": "string (required) - YAML specification to validate"
}
```

**Response Schema:**
```json
{
  "valid": "boolean - Whether the spec is valid",
  "errors": "array - List of validation errors if any",
  "status": "string - Validation status"
}
```

#### `POST /graphic/render`

Render a graphic from a YAML specification.

**Request Schema:**
```json
{
  "yaml_spec": "string (required) - YAML specification to render",
  "deterministic": "boolean (default: true) - Use deterministic rendering",
  "output_format": "string (default: 'svg') - Output format: svg, png, or both"
}
```

**Response Schema:**
```json
{
  "svg_content": "string - Rendered SVG content",
  "png_base64": "string (optional) - Base64 encoded PNG if requested",
  "status": "string - Rendering status",
  "error": "string (optional) - Error message if failed"
}
```

#### `POST /graphic/make`

Full pipeline: plan, validate, and render a graphic in one request.

**Request Schema:**
```json
{
  "brief": "string (required) - Natural language description",
  "width": "integer (default: 960) - Canvas width",
  "height": "integer (default: 640) - Canvas height",
  "model_name": "string (default: 'o3-mini') - LLM model to use",
  "output_format": "string (default: 'both') - Output format: svg, png, or both"
}
```

**Response Schema:**
```json
{
  "yaml_spec": "string - Generated YAML specification",
  "svg_content": "string - Rendered SVG content",
  "png_base64": "string (optional) - Base64 encoded PNG if requested",
  "job_id": "string (optional) - Job ID for async processing",
  "status": "string - Processing status",
  "error": "string (optional) - Error message if failed"
}
```

#### `GET /graphic/job/{job_id}`

Get status of an async job (placeholder for future async processing).

**Response Schema:**
```json
{
  "job_id": "string - Job identifier",
  "status": "string - Job status",
  "message": "string - Status message"
}
```

### 2. Prompt Routes (`/prompt`)

Process molecular prompts and generate various types of visualizations.

#### `POST /prompt/generate-from-pubchem/`

Generate molecule names from a natural language prompt.

**Request Schema:**
```json
{
  "prompt": "string (required) - Natural language prompt",
  "model": "string (optional) - Preferred model type"
}
```

**Response Schema:**
```json
{
  "molecule_names": "array[string] - List of generated molecule names"
}
```

#### `POST /prompt/validate-scientific/`

Validate if a prompt is scientific in nature.

**Request Schema:**
```json
{
  "prompt": "string (required) - Prompt to validate",
  "model": "string (optional) - Preferred model type"
}
```

**Response Schema:**
```json
{
  "is_molecular": "boolean - Whether the prompt is scientific/molecular"
}
```

#### `POST /prompt/generate-molecule-diagram/`

Generate a 2D molecule diagram from a text prompt.

**Request Schema:**
```json
{
  "prompt": "string (required) - Text description of the diagram",
  "canvas_width": "integer (default: 960) - Canvas width in pixels",
  "canvas_height": "integer (default: 640) - Canvas height in pixels",
  "model": "string (optional) - LLM model to use"
}
```

**Response Schema:**
```json
{
  "diagram_image": "string - Base64 PNG or SVG string",
  "diagram_plan": {
    "plan": "string - Description of the diagram plan",
    "molecule_list": "array - List of molecules and their positions",
    "arrows": "array - List of arrows connecting molecules",
    "canvas_width": "integer - Canvas width",
    "canvas_height": "integer - Canvas height"
  },
  "status": "string - Processing status (completed/failed/processing)",
  "job_id": "string (optional) - Job ID for async processing",
  "error": "string (optional) - Error message if failed"
}
```

#### `POST /prompt/fetch-molecule-data/`

Fetch molecule information and PDB data.

**Request Schema:**
```json
{
  "query": "string (required) - Molecule name or identifier"
}
```

**Response Schema:**
```json
{
  "molecule_data": "object - Molecule information and PDB data"
}
```

#### `POST /prompt/fetch-molecule-2d/`

Return 2D coordinate information for a molecule.

**Request Schema:**
```json
{
  "query": "string (required) - Molecule name or identifier"
}
```

**Response Schema:**
```json
{
  "query": "string - Original query",
  "status": "string - Processing status"
}
```

#### `POST /prompt/fetch-molecule-layout/`

Return 2D info for multiple molecules with layout boxes.

**Request Schema:**
```json
{
  "molecules": [
    {
      "query": "string - Molecule name or identifier",
      "box": {
        "x": "float - X coordinate",
        "y": "float - Y coordinate",
        "width": "float - Box width",
        "height": "float - Box height"
      }
    }
  ]
}
```

**Response Schema:**
```json
{
  "molecules": "array - List of molecules with layout information"
}
```

#### `POST /prompt/sdf-to-pdb/`

Convert SDF text to PDB format using RDKit.

**Request Schema:**
```json
{
  "sdf": "string (required) - SDF text to convert"
}
```

**Response Schema:**
```json
{
  "pdb_data": "string - Converted PDB data"
}
```

### 3. RCSB Routes (`/rcsb`)

Fetch protein structures and metadata from the RCSB Protein Data Bank.

#### `POST /rcsb/fetch-structure/`

Fetch a protein structure by identifier.

**Request Schema:**
```json
{
  "identifier": "string (required) - PDB identifier (e.g., '1ubq')",
  "format": "string (default: 'pdb') - Format: pdb or cif"
}
```

**Response Schema:**
```json
{
  "data": "string - Structure data in requested format"
}
```

#### `POST /rcsb/fetch-model/`

Fetch an AlphaFold model by UniProt ID.

**Request Schema:**
```json
{
  "uniprot_id": "string (required) - UniProt identifier",
  "format": "string (default: 'pdb') - Format: pdb or cif"
}
```

**Response Schema:**
```json
{
  "data": "string - Model data in requested format"
}
```

#### `GET /rcsb/entry/{identifier}`

Get metadata for a PDB entry.

**Response Schema:**
```json
{
  "metadata": "object - Entry metadata"
}
```

#### `GET /rcsb/annotations/{identifier}`

Get sequence annotations for a PDB entry.

**Response Schema:**
```json
{
  "annotations": "object - Sequence annotation data"
}
```

#### `POST /rcsb/computed-model/`

Fetch computed model information via GraphQL.

**Request Schema:**
```json
{
  "identifier": "string (required) - Structure identifier",
  "model_id": "string (required) - Model identifier"
}
```

**Response Schema:**
```json
{
  "metadata": "object - Computed model metadata"
}
```

#### `POST /rcsb/fetch-esm-model/`

Fetch an ESM (Evolutionary Scale Modeling) model.

**Request Schema:**
```json
{
  "uniprot_id": "string (required) - UniProt identifier",
  "format": "string (default: 'pdb') - Format: pdb or cif"
}
```

**Response Schema:**
```json
{
  "data": "string - ESM model data"
}
```

#### `POST /rcsb/upload-structure/`

Upload a structure for processing.

**Request Schema:**
```json
{
  "data": "string (required) - Structure data to upload",
  "filename": "string (default: 'upload.pdb') - Filename for the upload"
}
```

**Response Schema:**
```json
{
  "upload_id": "string - Unique identifier for the uploaded structure"
}
```

#### `POST /rcsb/align/`

Perform pairwise alignment between two structures.

**Request Schema:**
```json
{
  "identifier1": "string (required) - First structure identifier",
  "identifier2": "string (required) - Second structure identifier"
}
```

**Response Schema:**
```json
{
  "metadata": "object - Alignment results and metadata"
}
```

#### `GET /rcsb/group/{group_id}`

Get entries belonging to a specific group.

**Response Schema:**
```json
{
  "metadata": "object - Group entries and metadata"
}
```

#### `GET /rcsb/feature-annotations/{identifier}`

Get feature annotations for a structure.

**Response Schema:**
```json
{
  "annotations": "object - Feature annotation data"
}
```

### 4. Render Routes (`/render`)

Render 3D molecular structures using PyMOL with advanced rendering options.

#### `POST /render`

Render a molecular structure from a natural language description.

**Request Schema:**
```json
{
  "description": "string (required) - Natural language description of what to render",
  "format": "string (default: 'image') - Output format: image, model, or animation",
  "transparent_background": "boolean (default: false) - Use transparent background",
  "ray_trace": "boolean (default: true) - Enable ray tracing for high quality",
  "resolution": "array[int, int] (default: [1920, 1080]) - Output resolution",
  "dpi": "integer (default: 300) - DPI for image output",
  "ray_trace_mode": "string (default: 'default') - Ray trace mode: default, cartoon_outline, bw, poster",
  "antialias": "boolean (default: true) - Enable antialiasing",
  "ray_shadow": "boolean (default: true) - Enable ray shadows",
  "depth_cue": "boolean (default: true) - Enable depth cueing",
  "background_color": "string (default: 'white') - Background color"
}
```

**Response:**
- For small files: Direct file response with appropriate media type
- For large files (>25MB): JSON response with URL to static file

**Large File Response Schema:**
```json
{
  "url": "string - URL to access the rendered file",
  "metadata": "object - Rendering metadata including camera position, center, and bounding box"
}
```

**Headers:**
- `X-Metadata`: JSON string containing rendering metadata (for direct file responses)

## Data Models

### DiagramPlan
```json
{
  "plan": "string - Description of the diagram plan",
  "molecule_list": [
    {
      "molecule": "string - Molecule name or formula",
      "x": "float - X coordinate",
      "y": "float - Y coordinate",
      "width": "float (optional) - Width of molecule bounds",
      "height": "float (optional) - Height of molecule bounds",
      "label": "string (optional) - Display label",
      "label_position": "string (optional) - Label position: above, below, left, right"
    }
  ],
  "arrows": [
    {
      "start": "array[float, float] - Starting coordinates",
      "end": "array[float, float] - Ending coordinates",
      "text": "string (optional) - Arrow label text"
    }
  ],
  "canvas_width": "integer - Canvas width in pixels",
  "canvas_height": "integer - Canvas height in pixels"
}
```

### YAML Graphic Schema
The graphic routes use a YAML specification format with the following structure:

```yaml
meta:
  title: "string - Diagram title"
  version: "string - Version (e.g., '1.0.0')"

canvas:
  w: integer  # Width in pixels
  h: integer  # Height in pixels
  dpi: integer (optional)  # DPI setting

cells:
  - id: "string - Unique cell identifier"
    type: "string - Cell type: TEXT, DIAGRAM, IMAGE_GEN, COMPUTE, GROUP"
    bbox:
      x: number  # X coordinate
      y: number  # Y coordinate
      w: number  # Width
      h: number  # Height
    # Additional properties based on cell type
    nodes: []    # For DIAGRAM cells
    edges: []    # For DIAGRAM cells
```

## Error Handling

All endpoints return appropriate HTTP status codes:

- `200` - Success
- `400` - Bad Request (invalid input)
- `404` - Not Found
- `422` - Validation Error
- `429` - Rate Limited
- `500` - Internal Server Error

Error responses include a `detail` field with a descriptive error message.

## Environment Variables

- `OPENAI_API_KEY` - Required for LLM-powered features
- `ENVIRONMENT` - Set to "development" or "production" for CORS configuration
- `PYMOL_QUIET` - Set to "1" for headless PyMOL operation
- `PYMOL_HEADLESS` - Set to "1" for headless PyMOL operation

## Rate Limiting

The server implements rate limiting for external API calls, particularly when fetching molecular data from PubChem and other external services.

## Caching

The render endpoint implements caching to improve performance for repeated requests. Cached files are stored in `/tmp/moleculens_cache/`.

## Security

- All PyMOL commands are validated before execution to prevent code injection
- CORS is configured based on environment settings
- File uploads are validated and sanitized

## Examples

### Generate a Simple Diagram

```bash
curl -X POST "http://localhost:8000/graphic/make" \
  -H "Content-Type: application/json" \
  -d '{
    "brief": "Show water H2O splitting into hydrogen H2 and oxygen O2",
    "width": 800,
    "height": 400
  }'
```

### Render a Protein Structure

```bash
curl -X POST "http://localhost:8000/render" \
  -H "Content-Type: application/json" \
  -d '{
    "description": "Show the structure of ubiquitin with cartoon representation",
    "format": "image",
    "ray_trace": true,
    "resolution": [1920, 1080]
  }'
```

### Fetch Protein Structure

```bash
curl -X POST "http://localhost:8000/rcsb/fetch-structure/" \
  -H "Content-Type: application/json" \
  -d '{
    "identifier": "1ubq",
    "format": "pdb"
  }'
```

### Validate Scientific Prompt

```bash
curl -X POST "http://localhost:8000/prompt/validate-scientific/" \
  -H "Content-Type: application/json" \
  -d '{
    "prompt": "Show me the structure of caffeine"
  }'
```
