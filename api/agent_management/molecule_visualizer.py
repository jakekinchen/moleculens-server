#!/usr/bin/env python3
"""
Molecule Visualizer Module

Refactored to:
- Convert SDF data (as a string) into PDB data in-memory
- Generate minimal JavaScript code that can be embedded in another Three.js scene
- Generate a fully standalone HTML that initializes a Three.js scene and loads the molecule
"""

import re
from rdkit import Chem
from rdkit.Chem import AllChem

class MoleculeVisualizer:
    """
    A class to convert and visualize molecular structures by:
      1. Converting SDF text to PDB data in-memory
      2. Generating minimal JavaScript for embedding in an existing Three.js scene
      3. Generating a full standalone HTML viewer
    """

    @staticmethod
    def _escape_js_string(text: str) -> str:
        """
        Safely escape backslashes, backticks, and dollar signs for embedding in JS strings.
        """
        return text.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')

    @staticmethod
    def _get_molecule_label_style() -> str:
        """
        Returns the CSS styling for the molecule label.
        Used by both standalone HTML and embedded JS versions.
        """
        return """
            #container {
                position: relative;
                width: 100%;
                height: 100vh;
            }
            #molecule-label {
                position: absolute;
                bottom: 20px;
                left: 50%;
                transform: translateX(-50%);
                color: white;
                font-family: Arial, sans-serif;
                font-size: 18px;
                background-color: rgba(0, 0, 0, 0.7);
                padding: 10px 20px;
                border-radius: 5px;
                z-index: 1000;
                pointer-events: none;
                text-align: center;
                max-width: 80%;
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
            }
        """

    @staticmethod
    def _get_molecule_label_html(name: str) -> str:
        """
        Returns the HTML for the molecule label.
        Used by both standalone HTML and embedded JS versions.
        """
        return f'<div id="molecule-label">{name}</div>'

    @staticmethod
    def _sdf_to_pdb_block(sdf_data: str) -> str:
        """
        Convert SDF data (string) to a single PDB block using RDKit in-memory.

        Returns an empty string if conversion fails.
        """
        mol = Chem.MolFromMolBlock(sdf_data, sanitize=True, removeHs=False)
        if mol is None:
            return ""

        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        AllChem.MMFFOptimizeMolecule(mol)

        pdb_data = Chem.MolToPDBBlock(mol)
        return pdb_data if pdb_data else ""

    @classmethod
    def generate_js_code_from_sdf(cls, sdf_data: str, name: str = "Molecule") -> str:
        """
        Return minimal JavaScript snippet that:
        - Converts SDF to PDB in Python, embed it here
        - Loads with PDBLoader (assumes THREE.PDBLoader is available)
        - Creates atoms & bonds, adds them to the scene
        - No custom cameras or lights; relies on the parent environment
        - Perfect for embedding in React or other Three.js setups
        """

        pdb_data = cls._sdf_to_pdb_block(sdf_data)
        if not pdb_data:
            return "// GeometryAgent LLM-generated code\n// Error: Could not convert SDF to PDB"

        escaped_pdb = cls._escape_js_string(pdb_data)
        label_style = cls._escape_js_string(cls._get_molecule_label_style())
        label_html = cls._escape_js_string(cls._get_molecule_label_html(name))

        # Use double braces for JS blocks
        return f"""// GeometryAgent LLM-generated code
function createMoleculeVisualization(THREE, scene) {{
    // Create a group for the molecule
    const root = new THREE.Group();
    scene.add(root);

    // Add molecule label styles if not already present
    if (!document.getElementById('molecule-label-style')) {{
        const style = document.createElement('style');
        style.id = 'molecule-label-style';
        style.textContent = `{label_style}`;
        document.head.appendChild(style);
    }}

    // Add molecule label if not already present
    if (!document.getElementById('molecule-label')) {{
        const container = document.querySelector('#container');
        if (container) {{
            const labelContainer = document.createElement('div');
            labelContainer.innerHTML = `{label_html}`;
            container.appendChild(labelContainer.firstChild);
        }}
    }}

    // Convert SDF -> PDB in Python, embed it here
    const pdbData = `{escaped_pdb}`;
    
    // Create and configure the PDB loader
    let loader;
    if (typeof THREE !== 'undefined' && THREE.PDBLoader) {{
        loader = new THREE.PDBLoader();
    }} else if (typeof window !== 'undefined' && window.PDBLoader) {{
        // If we manually attached PDBLoader to the window
        loader = new window.PDBLoader();
    }} else if (typeof PDBLoader !== 'undefined') {{
        loader = new PDBLoader();
    }} else {{
        console.error('PDBLoader not found. Make sure it is loaded first.');
        return root;
    }}
    
    const pdbBlob = new Blob([pdbData], {{ type: 'text/plain' }});
    const pdbUrl = URL.createObjectURL(pdbBlob);

    // Load and process the PDB data
    loader.load(pdbUrl, function (pdb) {{
        const geometryAtoms = pdb.geometryAtoms;
        const geometryBonds = pdb.geometryBonds;
        const json = pdb.json;

        const sphereGeometry = new THREE.IcosahedronGeometry(1, 3);
        const boxGeometry = new THREE.BoxGeometry(1, 1, 1);
        const offset = new THREE.Vector3();

        geometryAtoms.computeBoundingBox();
        geometryAtoms.boundingBox.getCenter(offset).negate();

        geometryAtoms.translate(offset.x, offset.y, offset.z);
        geometryBonds.translate(offset.x, offset.y, offset.z);

        let positions = geometryAtoms.getAttribute('position');
        const colors = geometryAtoms.getAttribute('color');

        const position = new THREE.Vector3();
        const color = new THREE.Color();

        // Add atoms
        for (let i = 0; i < positions.count; i++) {{
            position.x = positions.getX(i);
            position.y = positions.getY(i);
            position.z = positions.getZ(i);

            color.r = colors.getX(i);
            color.g = colors.getY(i);
            color.b = colors.getZ(i);

            const material = new THREE.MeshPhongMaterial({{ color: color }});
            const object = new THREE.Mesh(sphereGeometry, material);
            object.position.copy(position);
            object.position.multiplyScalar(20);
            object.scale.multiplyScalar(7.5);
            root.add(object);
        }}

        // Add bonds
        positions = geometryBonds.getAttribute('position');
        const start = new THREE.Vector3();
        const end = new THREE.Vector3();

        for (let i = 0; i < positions.count; i += 2) {{
            start.x = positions.getX(i);
            start.y = positions.getY(i);
            start.z = positions.getZ(i);

            end.x = positions.getX(i + 1);
            end.y = positions.getY(i + 1);
            end.z = positions.getZ(i + 1);

            start.multiplyScalar(20);
            end.multiplyScalar(20);

            const object = new THREE.Mesh(
                boxGeometry,
                new THREE.MeshPhongMaterial({{ color: 0xffffff }})
            );
            object.position.copy(start);
            object.position.lerp(end, 0.5);
            object.scale.set(5, 5, start.distanceTo(end));
            object.lookAt(end);
            root.add(object);
        }}

        // Clean up
        URL.revokeObjectURL(pdbUrl);
    }});

    // Return the root group for external control
    return root;
}}

// Execute the function to create the visualization
createMoleculeVisualization(THREE, scene);
"""

    @classmethod
    def generate_html_viewer_from_sdf(cls, sdf_data: str, name: str = "Molecule") -> str:
        """
        Generate a standalone HTML file that:
        - Embeds Three.js (ESM)
        - Sets up the camera, lights, and controls
        - Converts SDF -> PDB in-memory, loads it, and displays the molecule
        - Provides a full animation loop with auto-rotation
        - Fills the browser window
        - Displays the molecule name
        """

        pdb_data = cls._sdf_to_pdb_block(sdf_data)
        if not pdb_data:
            return f"<!DOCTYPE html><html><body><h1>Error converting SDF to PDB data for {name}</h1></body></html>"

        escaped_pdb = cls._escape_js_string(pdb_data)

        # Use double braces for JS blocks
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <title>Molecule Viewer - {name}</title>
    <meta charset="utf-8">
    <meta name="viewport"
          content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
    <style>
        body {{ margin: 0; padding: 0; overflow: hidden; }}
        {cls._get_molecule_label_style()}
    </style>
</head>
<body>
<div id="container"></div>
{cls._get_molecule_label_html(name)}
<script async src="https://unpkg.com/es-module-shims@1.8.0/dist/es-module-shims.js"></script>
<script type="importmap">
{{
  "imports": {{
    "three": "https://unpkg.com/three@0.156.1/build/three.module.js",
    "three/addons/": "https://unpkg.com/three@0.156.1/examples/jsm/"
  }}
}}
</script>
<script type="module">

import * as THREE from 'three';
import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';
import {{ PDBLoader }} from 'three/addons/loaders/PDBLoader.js';

let camera, scene, renderer, controls;
let root;
const rotationSpeed = 0.001; // Speed of auto-rotation (radians per frame)

init();
animate();

function init() {{
    const container = document.getElementById('container');

    // Scene
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x050505);

    // Camera
    camera = new THREE.PerspectiveCamera(70, window.innerWidth / window.innerHeight, 1, 5000);
    camera.position.z = 1000;

    // Lights
    const light1 = new THREE.DirectionalLight(0xffffff, 2.5);
    light1.position.set(1, 1, 1);
    scene.add(light1);

    const light2 = new THREE.DirectionalLight(0xffffff, 1.5);
    light2.position.set(-1, -1, 1);
    scene.add(light2);

    // Root group
    root = new THREE.Group();
    scene.add(root);

    // Renderer
    renderer = new THREE.WebGLRenderer({{ antialias: true }});
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth, window.innerHeight);
    container.appendChild(renderer.domElement);

    // Controls
    controls = new OrbitControls(camera, renderer.domElement);
    controls.minDistance = 200;
    controls.maxDistance = 3000;
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;

    // Handle resize
    window.addEventListener('resize', onWindowResize);

    // Load molecule
    const pdbData = `{escaped_pdb}`;
    const pdbBlob = new Blob([pdbData], {{ type: 'text/plain' }});
    const pdbUrl = URL.createObjectURL(pdbBlob);

    // Make PDBLoader available globally for compatibility with both approaches
    window.PDBLoader = PDBLoader;
    const loader = new PDBLoader();
    loader.load(pdbUrl, (pdb) => {{
        const geometryAtoms = pdb.geometryAtoms;
        const geometryBonds = pdb.geometryBonds;
        const boxGeometry = new THREE.BoxGeometry(1, 1, 1);
        const sphereGeometry = new THREE.IcosahedronGeometry(1, 3);
        const offset = new THREE.Vector3();

        geometryAtoms.computeBoundingBox();
        geometryAtoms.boundingBox.getCenter(offset).negate();

        geometryAtoms.translate(offset.x, offset.y, offset.z);
        geometryBonds.translate(offset.x, offset.y, offset.z);

        let positions = geometryAtoms.getAttribute('position');
        const colors = geometryAtoms.getAttribute('color');
        const position = new THREE.Vector3();
        const color = new THREE.Color();

        // Add atoms
        for (let i = 0; i < positions.count; i++) {{
            position.set(
                positions.getX(i),
                positions.getY(i),
                positions.getZ(i)
            );
            color.setRGB(
                colors.getX(i),
                colors.getY(i),
                colors.getZ(i)
            );

            const material = new THREE.MeshPhongMaterial({{ color: color }});
            const atom = new THREE.Mesh(sphereGeometry, material);
            atom.position.copy(position).multiplyScalar(75);
            atom.scale.setScalar(25);
            root.add(atom);
        }}

        // Add bonds
        positions = geometryBonds.getAttribute('position');
        const start = new THREE.Vector3();
        const end = new THREE.Vector3();

        for (let i = 0; i < positions.count; i += 2) {{
            start.set(
                positions.getX(i),
                positions.getY(i),
                positions.getZ(i)
            ).multiplyScalar(75);

            end.set(
                positions.getX(i+1),
                positions.getY(i+1),
                positions.getZ(i+1)
            ).multiplyScalar(75);

            const bondMesh = new THREE.Mesh(
                boxGeometry,
                new THREE.MeshPhongMaterial({{ color: 0xffffff }})
            );
            bondMesh.position.copy(start).lerp(end, 0.5);
            bondMesh.scale.set(5, 5, start.distanceTo(end));
            bondMesh.lookAt(end);
            root.add(bondMesh);
        }}

        // Clean up
        URL.revokeObjectURL(pdbUrl);
    }});
}}

function onWindowResize() {{
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}}

function animate() {{
    requestAnimationFrame(animate);
    
    // Auto-rotate the molecule
    if (root) {{
        root.rotation.y += rotationSpeed;
    }}
    
    controls.update();
    renderer.render(scene, camera);
}}
</script>
</body>
</html>
"""

def main():
    """
    Simple CLI usage:
      python molecule_visualizer.py <input.sdf> [html|js]
    """
    import sys
    if len(sys.argv) < 2:
        print("Usage: python molecule_visualizer.py <input.sdf> [html|js]")
        sys.exit(1)

    sdf_file = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "html"

    with open(sdf_file, "r") as f:
        sdf_data = f.read()

    name = "Molecule"
    if mode == "js":
        code = MoleculeVisualizer.generate_js_code_from_sdf(sdf_data, name)
        print(code)
    else:
        html = MoleculeVisualizer.generate_html_viewer_from_sdf(sdf_data, name)
        print(html)

if __name__ == "__main__":
    main()