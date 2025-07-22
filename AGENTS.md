# Moleculens Server: Past ▸ Present ▸ Future

## Past (previous_solutions)
* 2025-07-18  Initial RCSBAgent + `/fetch-structure` endpoint ✅
* 2025-07-25  Expanded RCSB endpoints (AlphaFold, annotations, uploads) ✅
* 2025-07-30  Alignment / group / feature routes added + tests ✅
* 2025-08-01  PyMOL template helpers (`overview_scene`, `binding_site_scene`) ✅
* 2025-08-02  Celery + Redis skeleton & echo task ✅
* 2025-08-03  Keyword prompt→template translator, lint baseline ✅
* 2025-07-21  Translator wired into render flow ✅
* 2025-07-21  Added `mutation_focus_scene` template ✅
* 2025-07-22  Celery render_scene task → produce glTF/USDZ ✅
* 2025-08-04  Tighten flake8/mypy rules; re-enable strict hook ✅

---

## Present (current_problem)
Exactly one open item keeps all agents aligned:

| Area | Task | Owner | Status |
|------|------|-------|--------|
| Auth | **Token mint & verify** | Backend | 🛠 |

*CI fails if this table becomes empty.*

---

## Future (backlog + next_steps)
1. Assess whether any of the following should be integrated or are already integrated "Backend to-do list—keep it lean:
	1.	POST /scenes
	•	Parse payload, run PyMOL render.
	•	Free plan → skip Supabase upload; write bytes to TTLCache.
	•	Pro plan → supabase.storage.from("scenes").upload(...); return signed URL.
	2.	GET /viewer/asset
	•	Query param token; verify.
	•	If free → cache.get(fp) or regenerate(); stream StreamingResponse.
	•	If pro  → 302 redirect to Supabase signed URL.
	3.	Cache layer
	•	Simple TTLCache(maxsize=500, ttl=6*3600) in module scope.
	•	Wrap regeneration in asyncio.Lock keyed by fingerprint.
	4.	Rate-limit & auth
	•	FastAPI dependency: limiter = Limiter(key_func=ip).
	•	30 req/min IP cap on asset route.
	•	Bearer JWT (user auth) only on /scenes; asset route solely token-based.
	5.	Headers
	•	Cache-Control: public,max-age=21600 on asset response.
	•	Content-Type dynamic: model/gltf+json or application/json (Mol*).
	•	CSP/HSTS added in middleware/security.py.
	6.	Supabase helper (if enabled)
	•	supabase_client = create_client(url, service_key) in a singleton.
	•	get_signed_url(path, 86400) for pro tier.
	7.	Metrics & logs
	•	Prometheus counters: scene_render_seconds, asset_requests_total, token_invalid_total.
	•	Log structured JSON (fp, plan, ms, status) to stdout.

# AGENTS.md Instructions
Once a task is completed, it needs to be moved from the "Present" section to the "Past" section below. If there are no remaining tasks in the "Present" section, then move one or more tasks from the "Future" section to the "Present" section. If there are no tasks in the "Future" section or in the "Present" section, then the document is considered complete and you should return it as is. If you notice that a step or tasks has already been completed in the "Present" or "Future" section, then you should change it to the "Past" section and update the date to the current date.

## Appendix A – At-a-Glance Status (auto-generated monthly)

| Area | Task | Status |
|------|------|--------|
| RCSB Data | fetch_structure endpoint | ✅ |
| RCSB Data | fetch_alphafold_model endpoint | ✅ |
| RCSB Data | fetch_entry_metadata endpoint | ✅ |
| RCSB Data | fetch_sequence_annotations endpoint | ✅ |
| RCSB Data | fetch_graphql_model endpoint | ✅ |
| RCSB Data | fetch_esmf_model endpoint | ✅ |
| RCSB Data | upload_structure endpoint | ✅ |
| RCSB Data | fetch_pairwise_alignment endpoint | ✅ |
| RCSB Data | fetch_group_entries endpoint | ✅ |
| RCSB Data | fetch_feature_annotations endpoint | ✅ |
| PyMOL Scenes | Scene template library (overview, binding_site, mutation) | ✅ |
| Testing | Coverage for all RCSB endpoints | ✅ |

Legend: ✅ completed 🛠 in-progress ❑ todo

Current Codebase

The existing backend is a FastAPI service designed for molecular visualization.  Its key components include:
	•	Headless PyMOL integration.  The Docker image installs the open-source PyMOL (via pymol-open-source), assimp and ffmpeg, and the server launches PyMOL headlessly on startup .  A render endpoint uses a description-to-PyMOL-commands translator to convert free-form text into deterministic PyMOL commands.  Allowed commands include fetching structures, changing representations, colouring, orienting, saving PNGs and PDBs, etc. .  The endpoint holds a lock around PyMOL commands and renders either images, PDB models or simple animations .
	•	LLM-driven command generation.  A small LLM wrapper converts natural language descriptions into PyMOL commands.  If the LLM call fails, a fallback script loads ubiquitin and renders a coloured cartoon .  The security module whitelists only non-destructive PyMOL commands .
	•	PubChem integration and RDKit utilities.  The repository contains a comprehensive PubChemAgent which uses PubChem's REST API and RDKit to retrieve chemical structures, convert SDF files to PDB blocks and compute chemical properties.  This agent can convert molecules into 3D conformers and embed them into PDB format using RDKit's MolToPDBBlock .  An in-progress plan (MOLECULE_DIAGRAM_PLAN.md) describes adding an endpoint for 2-D molecular diagrams rendered by RDKit, with Pydantic models and front-end components .
	•	Modular agents and orchestration.  Additional agents handle prompt orchestration, geometry generation, animation and diagram rendering.  These orchestrate RDKit, PyMOL and LLM services to produce multi-modal outputs.  However, there is no direct integration with RCSB's services or web-based interactive viewers.

Additional Research: Emerging Tools and APIs

RCSB / Mol* ecosystem
	•	Mol replaces NGL.*  RCSB is deprecating the NGL viewer; users must switch to Mol*, which provides enhanced features and improved performance for 3D visualization .
	•	Improved rendering and geometry export.  Mol* now supports artifact-free transparency, improved depth perception and cleaner outlines.  It can export molecular scenes as glTF, STL and OBJ files for use in external rendering programs and 3-D printing .  This capability paves the way for AR/VR experiences and 3-D object generation.
	•	User-upload and structure motif search.  The RCSB Mol* viewer can upload user-provided structures (PDB, mmCIF or BinaryCIF) via a dedicated API.  Uploaded files receive a shareable URL and can be used for structure similarity and motif searches; the resulting alignment can be exported as a ZIP archive .
	•	Sequence Coordinates Service.  A new Sequence Coordinates Service, replacing the 1D Coordinates service by May 31 2025, integrates growing numbers of structures and makes minor API changes .  Developers are advised to migrate to this service.
	•	Enhanced pairwise alignment and feature mapping.  The pairwise structure alignment tool now synchronizes highlighting between aligned sequences and superposed structures and allows toggling of polymeric chains; the alignment API is accessible from the UI .  Users can also search by UniProt IDs to fetch AlphaFold or ESM predicted models for comparison .
	•	Grouping and sequence annotation tools.  RCSB's Advanced Search groups similar structures and exposes API endpoints for group exploration .  The Sequence Annotations Viewer maps domains, motifs and modifications onto structures and is integrated with Mol* for bi-directional exploration.

PyMOL community enhancements
	•	APBS and interaction plugins.  PyMOL has plugins for electrostatics (APBS), π–π and cation-π interaction detection and water network mapping, which can automatically annotate structures.
	•	Headless rendering.  PyMOL's Python API supports ray-traced, publication-quality rendering of images and movies.  It can export surfaces and electron density maps via plugins.
	•	Surface extraction and glTF conversion.  PyMOL can export OBJ meshes; using assimp and tools such as gltfpack, these meshes can be converted to glTF and USDZ, enabling AR/VR deployment.

RDKit and Cheminformatics
	•	RDKit generates 2-D depictions (MolDraw2D), calculates physicochemical properties and creates 3-D conformers when experimental coordinates are unavailable.  It also converts molecules to PDB or mol2 formats, which PyMOL can read.

XR/VR opportunities
	•	glTF/GLB and USDZ exports allow molecules to be viewed natively on Meta Quest or Apple Vision Pro headsets.  Multi-resolution meshes and metadata can enable interactive highlighting and selection in XR.

Path to Expand the Codebase

The goal is to move from an LLM-driven PyMOL renderer to a full-service platform that delivers high-quality static graphics, interactive web scenes and XR-ready assets on demand.  The following steps outline this expansion.

1. Robust Data Retrieval Layer
	1.	Integrate RCSB APIs.  Add modules to fetch PDB/mmCIF files, computed structure models and annotations via RCSB's REST/GraphQL endpoints.  Use the new Sequence Coordinates Service to retrieve residue-level annotations.  Provide functions to fetch AlphaFold and ESM models using UniProt or MGnify IDs as described in the updated alignment tool .
	2.	Enhance PubChem agent.  Extend the existing PubChemAgent to also query RCSB's Ligand module (until it is retired) and the Chemical Component Dictionary.  Use RDKit to sanitize and prepare ligand structures for 3-D depiction.
	3.	User upload support.  Implement an endpoint that accepts user-uploaded PDB/mmCIF files and uses RCSB's user-upload API to obtain shareable URLs for motif and similarity searches .

2. Static Rendering Pipeline (PyMOL)
	1.	Scene templates and modules.  Replace the current generic LLM translator with a library of parameterized PyMOL scene functions (e.g., overview, binding site, mutation focus) that accept a structure ID, selection strings and styling options.  These functions encapsulate best-practice settings (cartoon + transparent surface, hydrophobic coloring, etc.).
	2.	Feature annotation.  Incorporate PyMOL plugins for electrostatics (APBS) and interaction detection.  Automatically compute distances, hydrogen bonds and π interactions and annotate them in the scene.
	3.	Model export options.  When the user requests a 3-D model, call cmd.save to export an OBJ mesh and convert it to glTF (using assimp → glTF or gltfpack) and USDZ.  Provide options for mesh decimation for XR use.
	4.	Ray-traced movies.  Use PyMOL's scene and ray commands to generate smooth fly-through animations, optionally guided by user prompts.  Leverage ffmpeg (already installed) to assemble PNG frames into MP4.

3. Web-Based Interactive Visualization
	1.	Embed Mol*.  Integrate Mol* into the front end (Next.js app) as the default viewer.  Provide an API that returns a pre-configured Mol* state in JSON or encoded URL parameters (structure ID, selection, coloring, surfaces).  Leverage Mol*'s ability to export scenes as glTF and to align multiple structures interactively .
	2.	Sequence annotations and ligand viewer.  Use RCSB's Sequence Annotations Viewer and Ligand Explorer within the app to allow users to click on features and see them highlighted in 3-D.  When the backend computes features (mutations, domains, post-translational modifications), encode them in the Mol* state.
	3.	Interactive alignment and group exploration.  Add endpoints that call the RCSB alignment API and return alignment results.  Provide UI components that display aligned sequences and allow toggling of chains, replicating the improved alignment tool .
	4.	High-resolution export.  Surface a client-side "Save Image" feature that uses Mol*'s built-in image export to produce PNG/TIFF images at user-specified DPI.

4. Automated Pipeline and Orchestration
	1.	Job orchestration.  Use an asynchronous task queue (e.g., Celery) to handle long-running rendering and conversion jobs.  Jobs accept JSON payloads specifying structure IDs, ligand selections, view templates and output formats.
	2.	Caching and storage.  Extend the existing caching mechanism to store not only PNGs but also glTF/GLB, USDZ and MP4 outputs.  Use object storage (e.g., S3) and return signed URLs for retrieval.
	3.	Metadata enrichment.  Include in responses metadata such as camera view matrices, bounding boxes and sequence feature mappings.  This helps clients replicate views and implement XR interactions.

5. XR/VR Delivery (Most Lucrative Path)

The ability to deliver immersive molecular experiences differentiates the platform and opens new revenue streams (educational licensing, pharma presentations and VR analytics).  The following actions provide the highest return:
	1.	Mesh export and LOD generation.  Automate extraction of molecular surfaces from PyMOL (cmd.save('obj')) and convert them to glTF/GLB with Draco compression for Meta Quest and USDZ for Apple Vision Pro.  Generate multiple levels of detail and include vertex-to-atom index mappings so clients can highlight residues in XR.
	2.	XR viewer integration.  Develop front-end modules using WebXR (Babylon.js or Three.js) to load glTF assets and implement basic interactions (rotate, zoom, residue pick).  On Apple devices, use RealityKit to load USDZ files.  Provide fallback to Mol* for browsers.
	3.	Remote streaming for huge structures.  For assemblies with >200 k atoms, implement cloud rendering and stream H.265 video to the headset (CloudXR).  Users still send selection events back to the server.
	4.	AR annotations and analytics.  Overlay labels, interaction distances and domain boundaries in XR.  Enable measurement tools similar to Mol* so users can interrogate structures in 3-D space.

6. Optional Enhancements
	•	ChimeraX (or regular Chimera) integration.  Offer high-end GPU rendering and VR features by running ChimeraX headless for certain scenes (e.g., cryo-EM segmentation).  Use --cmd to load a structure and export PNGs or glTF when PyMOL's visuals are insufficient.
	•	Cheminformatics metrics.  Annotate ligand 2-D diagrams with properties (molecular weight, logP) computed by RDKit.  Provide cross-links between 2-D diagrams and 3-D views.
	•	Automated reporting.  Generate PDF/PowerPoint reports summarizing binding site features, interactions and comparative analyses using Python libraries.

Conclusion

The Moleculens server already leverages PyMOL and RDKit to render molecular images, but it lacks integration with modern RCSB tools and does not support interactive or XR-ready outputs.  Recent improvements to RCSB's Mol* viewer—improved rendering, glTF export and user-upload functionality—provide an opportunity to build a powerful on-demand graphics service.  By systematically adding data retrieval layers, modular rendering templates, Mol* integration and XR/VR export capabilities, the platform can evolve into a full-fledged molecular visualization pipeline.  The most lucrative path involves investing in XR/VR outputs via glTF and USDZ export, enabling immersive experiences for education, research and corporate clients.

### July 2025 Progress: Initial RCSB Integration

* Implemented a simple `RCSBAgent` that downloads PDB or mmCIF files from the RCSB repository.
* Added a `/rcsb/fetch-structure/` endpoint returning the raw structure data.
* Registered this agent in `AgentFactory` and extended `AgentType` with a `RCSB` option.
* Created tests that stub heavy dependencies (RDKit, OpenAI) and validate both geometry and RCSB routes.
* Pinned `httpx` to a compatible version so Starlette's `TestClient` works correctly during testing.

## RCSB Agent expansion (July 18 2025)

Following the initial integration of the RCSB agent, the back-end now exposes additional capabilities to support computed structure models and metadata retrieval:

* Added `fetch_alphafold_model` and `fetch_entry_metadata` methods to `RCSBAgent`.  These allow downloading predicted structures from the AlphaFold database (using UniProt accessions) and retrieving JSON metadata for PDB entries via the RCSB Data API.  A common `_check_format` helper normalizes the PDB/mmCIF format argument.
* Extended `api/routers/rcsb/routes.py` with two new endpoints:
  * `POST /rcsb/fetch-model/` accepts a UniProt ID and format (pdb or cif) and returns a predicted model downloaded from the AlphaFold DB.
  * `GET /rcsb/entry/{identifier}` returns metadata for a PDB entry by delegating to the Data API.
* Added a default `AgentModelConfig` for the `rcsb` agent type to `agent_model_config.py` to ensure consistent registration in the factory.
* Added `tests/test_rcsb_extensions.py` which stubs external dependencies and verifies the new endpoints.  These tests follow the pattern used for geometry and RCSB structure tests by injecting fake modules and patching agent methods to avoid network calls.
* The original `fetch_structure` endpoint and tests continue to function unchanged.

This update progresses the data retrieval layer outlined in the expansion plan and lays groundwork for subsequent features like sequence annotations, alignments and user uploads.

## RCSB data layer enhancements (July 25 2025)

Building on the earlier work, the server now exposes additional endpoints that tap into the Sequence Coordinates Service and GraphQL API while paving the way for user contributed data:

* Added `fetch_sequence_annotations`, `fetch_graphql_model`, `fetch_esmf_model` and `upload_structure` methods to `RCSBAgent`.
* Updated `api/routers/rcsb/routes.py` with new routes:
  * `GET /rcsb/annotations/{identifier}` returns residue-level annotations.
  * `POST /rcsb/computed-model/` queries the GraphQL API for a computed model.
  * `POST /rcsb/fetch-esm-model/` retrieves an ESMFold prediction.
  * `POST /rcsb/upload-structure/` uploads user data and returns a shareable ID.
* Added integration tests that stub network access to verify the new endpoints.

Next steps are to integrate ligand data from the CCD via the `PubChemAgent`, define reusable PyMOL scene templates and begin Mol* embedding for interactive viewing.

## July 21 2025 Follow-Up

The data layer expansion outlined above has been fully implemented:

* Added `fetch_sequence_coordinates` to `RCSBAgent` and introduced the `/rcsb/sequence-coordinates/{identifier}` endpoint with tests.
* Implemented `/rcsb/align/`, `/rcsb/group/{group_id}` and `/rcsb/feature-annotations/{identifier}` routes backed by new agent methods.
* Added `/rcsb/upload-structure/` for user uploads and verified it with integration tests.
* Created a scene template library and a lightweight prompt translator with comprehensive unit tests.

These features pave the way for upcoming Mol* integration and XR/VR export capabilities [oai_citation:8‡GitHub](https://github.com/jakekinchen/moleculens-server/blob/c032a23995e32c9aafb72cd96d2c894e789b0bb8/AGENTS.md#L32-L65).

### pymol_prompt_parser
* Stateless; uses OpenAI function calling (`build_scene_request`) to convert free-form English into a `SceneSpec` JSON payload.

### pymol_command_builder
* Accepts `SceneSpec`, dispatches to template helpers (`overview_scene`, `binding_site_scene`, `mutation_scene`).
* Raises `SceneValidationError` (Pydantic) on invalid specs.
