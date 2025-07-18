Expanding Moleculens Server: Current State and Path to On‑Demand Molecular Graphics

Current Codebase

The existing backend is a FastAPI service designed for molecular visualization.  Its key components include:
	•	Headless PyMOL integration.  The Docker image installs the open‑source PyMOL (via pymol-open-source), assimp and ffmpeg, and the server launches PyMOL headlessly on startup ￼.  A render endpoint uses a description‑to‑PyMOL‑commands translator to convert free‑form text into deterministic PyMOL commands.  Allowed commands include fetching structures, changing representations, colouring, orienting, saving PNGs and PDBs, etc. ￼.  The endpoint holds a lock around PyMOL commands and renders either images, PDB models or simple animations ￼.
	•	LLM‑driven command generation.  A small LLM wrapper converts natural language descriptions into PyMOL commands.  If the LLM call fails, a fallback script loads ubiquitin and renders a coloured cartoon ￼.  The security module whitelists only non‑destructive PyMOL commands ￼.
	•	PubChem integration and RDKit utilities.  The repository contains a comprehensive PubChemAgent which uses PubChem’s REST API and RDKit to retrieve chemical structures, convert SDF files to PDB blocks and compute chemical properties.  This agent can convert molecules into 3D conformers and embed them into PDB format using RDKit’s MolToPDBBlock ￼.  An in‑progress plan (MOLECULE_DIAGRAM_PLAN.md) describes adding an endpoint for 2‑D molecular diagrams rendered by RDKit, with Pydantic models and front‑end components ￼ ￼.
	•	Modular agents and orchestration.  Additional agents handle prompt orchestration, geometry generation, animation and diagram rendering.  These orchestrate RDKit, PyMOL and LLM services to produce multi‑modal outputs.  However, there is no direct integration with RCSB’s services or web‑based interactive viewers.

Additional Research: Emerging Tools and APIs

RCSB / Mol* ecosystem
	•	Mol replaces NGL.*  RCSB is deprecating the NGL viewer; users must switch to Mol*, which provides enhanced features and improved performance for 3D visualization ￼.
	•	Improved rendering and geometry export.  Mol* now supports artifact‑free transparency, improved depth perception and cleaner outlines.  It can export molecular scenes as glTF, STL and OBJ files for use in external rendering programs and 3‑D printing ￼.  This capability paves the way for AR/VR experiences and 3‑D object generation.
	•	User‑upload and structure motif search.  The RCSB Mol* viewer can upload user‑provided structures (PDB, mmCIF or BinaryCIF) via a dedicated API.  Uploaded files receive a shareable URL and can be used for structure similarity and motif searches; the resulting alignment can be exported as a ZIP archive ￼.
	•	Sequence Coordinates Service.  A new Sequence Coordinates Service, replacing the 1D Coordinates service by May 31 2025, integrates growing numbers of structures and makes minor API changes ￼.  Developers are advised to migrate to this service.
	•	Enhanced pairwise alignment and feature mapping.  The pairwise structure alignment tool now synchronizes highlighting between aligned sequences and superposed structures and allows toggling of polymeric chains; the alignment API is accessible from the UI ￼.  Users can also search by UniProt IDs to fetch AlphaFold or ESM predicted models for comparison ￼.
	•	Grouping and sequence annotation tools.  RCSB’s Advanced Search groups similar structures and exposes API endpoints for group exploration ￼.  The Sequence Annotations Viewer maps domains, motifs and modifications onto structures and is integrated with Mol* for bi‑directional exploration.

PyMOL community enhancements
	•	APBS and interaction plugins.  PyMOL has plugins for electrostatics (APBS), π–π and cation‑π interaction detection and water network mapping, which can automatically annotate structures.
	•	Headless rendering.  PyMOL’s Python API supports ray‑traced, publication‑quality rendering of images and movies.  It can export surfaces and electron density maps via plugins.
	•	Surface extraction and glTF conversion.  PyMOL can export OBJ meshes; using assimp and tools such as gltfpack, these meshes can be converted to glTF and USDZ, enabling AR/VR deployment.

RDKit and Cheminformatics
	•	RDKit generates 2‑D depictions (MolDraw2D), calculates physicochemical properties and creates 3‑D conformers when experimental coordinates are unavailable.  It also converts molecules to PDB or mol2 formats, which PyMOL can read.

XR/VR opportunities
	•	glTF/GLB and USDZ exports allow molecules to be viewed natively on Meta Quest or Apple Vision Pro headsets.  Multi‑resolution meshes and metadata can enable interactive highlighting and selection in XR.

Path to Expand the Codebase

The goal is to move from an LLM‑driven PyMOL renderer to a full‑service platform that delivers high‑quality static graphics, interactive web scenes and XR‑ready assets on demand.  The following steps outline this expansion.

1. Robust Data Retrieval Layer
	1.	Integrate RCSB APIs.  Add modules to fetch PDB/mmCIF files, computed structure models and annotations via RCSB’s REST/GraphQL endpoints.  Use the new Sequence Coordinates Service to retrieve residue‑level annotations.  Provide functions to fetch AlphaFold and ESM models using UniProt or MGnify IDs as described in the updated alignment tool ￼.
	2.	Enhance PubChem agent.  Extend the existing PubChemAgent to also query RCSB’s Ligand module (until it is retired) and the Chemical Component Dictionary.  Use RDKit to sanitize and prepare ligand structures for 3‑D depiction.
	3.	User upload support.  Implement an endpoint that accepts user‑uploaded PDB/mmCIF files and uses RCSB’s user‑upload API to obtain shareable URLs for motif and similarity searches ￼.

2. Static Rendering Pipeline (PyMOL)
	1.	Scene templates and modules.  Replace the current generic LLM translator with a library of parameterized PyMOL scene functions (e.g., overview, binding site, mutation focus) that accept a structure ID, selection strings and styling options.  These functions encapsulate best‑practice settings (cartoon + transparent surface, hydrophobic coloring, etc.).
	2.	Feature annotation.  Incorporate PyMOL plugins for electrostatics (APBS) and interaction detection.  Automatically compute distances, hydrogen bonds and π interactions and annotate them in the scene.
	3.	Model export options.  When the user requests a 3‑D model, call cmd.save to export an OBJ mesh and convert it to glTF (using assimp → glTF or gltfpack) and USDZ.  Provide options for mesh decimation for XR use.
	4.	Ray‑traced movies.  Use PyMOL’s scene and ray commands to generate smooth fly‑through animations, optionally guided by user prompts.  Leverage ffmpeg (already installed) to assemble PNG frames into MP4.

3. Web‑Based Interactive Visualization
	1.	Embed Mol*.  Integrate Mol* into the front end (Next.js app) as the default viewer.  Provide an API that returns a pre‑configured Mol* state in JSON or encoded URL parameters (structure ID, selection, coloring, surfaces).  Leverage Mol*’s ability to export scenes as glTF and to align multiple structures interactively ￼.
	2.	Sequence annotations and ligand viewer.  Use RCSB’s Sequence Annotations Viewer and Ligand Explorer within the app to allow users to click on features and see them highlighted in 3‑D.  When the backend computes features (mutations, domains, post‑translational modifications), encode them in the Mol* state.
	3.	Interactive alignment and group exploration.  Add endpoints that call the RCSB alignment API and return alignment results.  Provide UI components that display aligned sequences and allow toggling of chains, replicating the improved alignment tool ￼.
	4.	High‑resolution export.  Surface a client‑side “Save Image” feature that uses Mol*’s built‑in image export to produce PNG/TIFF images at user‑specified DPI.

4. Automated Pipeline and Orchestration
	1.	Job orchestration.  Use an asynchronous task queue (e.g., Celery) to handle long‑running rendering and conversion jobs.  Jobs accept JSON payloads specifying structure IDs, ligand selections, view templates and output formats.
	2.	Caching and storage.  Extend the existing caching mechanism to store not only PNGs but also glTF/GLB, USDZ and MP4 outputs.  Use object storage (e.g., S3) and return signed URLs for retrieval.
	3.	Metadata enrichment.  Include in responses metadata such as camera view matrices, bounding boxes and sequence feature mappings.  This helps clients replicate views and implement XR interactions.

5. XR/VR Delivery (Most Lucrative Path)

The ability to deliver immersive molecular experiences differentiates the platform and opens new revenue streams (educational licensing, pharma presentations and VR analytics).  The following actions provide the highest return:
	1.	Mesh export and LOD generation.  Automate extraction of molecular surfaces from PyMOL (cmd.save('obj')) and convert them to glTF/GLB with Draco compression for Meta Quest and USDZ for Apple Vision Pro.  Generate multiple levels of detail and include vertex‑to‑atom index mappings so clients can highlight residues in XR.
	2.	XR viewer integration.  Develop front‑end modules using WebXR (Babylon.js or Three.js) to load glTF assets and implement basic interactions (rotate, zoom, residue pick).  On Apple devices, use RealityKit to load USDZ files.  Provide fallback to Mol* for browsers.
	3.	Remote streaming for huge structures.  For assemblies with >200 k atoms, implement cloud rendering and stream H.265 video to the headset (CloudXR).  Users still send selection events back to the server.
	4.	AR annotations and analytics.  Overlay labels, interaction distances and domain boundaries in XR.  Enable measurement tools similar to Mol* so users can interrogate structures in 3‑D space.

6. Optional Enhancements
	•	ChimeraX (or regular Chimera) integration.  Offer high‑end GPU rendering and VR features by running ChimeraX headless for certain scenes (e.g., cryo‑EM segmentation).  Use --cmd to load a structure and export PNGs or glTF when PyMOL’s visuals are insufficient.
	•	Cheminformatics metrics.  Annotate ligand 2‑D diagrams with properties (molecular weight, logP) computed by RDKit.  Provide cross‑links between 2‑D diagrams and 3‑D views.
	•	Automated reporting.  Generate PDF/PowerPoint reports summarizing binding site features, interactions and comparative analyses using Python libraries.

Conclusion

The Moleculens server already leverages PyMOL and RDKit to render molecular images, but it lacks integration with modern RCSB tools and does not support interactive or XR‑ready outputs.  Recent improvements to RCSB’s Mol* viewer—improved rendering, glTF export and user‑upload functionality—provide an opportunity to build a powerful on‑demand graphics service.  By systematically adding data retrieval layers, modular rendering templates, Mol* integration and XR/VR export capabilities, the platform can evolve into a full‑fledged molecular visualization pipeline.  The most lucrative path involves investing in XR/VR outputs via glTF and USDZ export, enabling immersive experiences for education, research and corporate clients.

### July 2025 Progress: Initial RCSB Integration

* Implemented a simple `RCSBAgent` that downloads PDB or mmCIF files from the RCSB repository.
* Added a `/rcsb/fetch-structure/` endpoint returning the raw structure data.
* Registered this agent in `AgentFactory` and extended `AgentType` with a `RCSB` option.
* Created tests that stub heavy dependencies (RDKit, OpenAI) and validate both geometry and RCSB routes.
* Pinned `httpx` to a compatible version so Starlette's `TestClient` works correctly during testing.
