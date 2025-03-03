#!/usr/bin/env python3
"""
Molecule Visualizer Module

Refactored to:
- Convert SDF data (as a string) into PDB data in-memory
- Generate minimal JavaScript code that can be embedded in another Three.js scene
- Generate a fully standalone HTML that initializes a Three.js scene and loads the molecule
- Add toggleable atomic annotations that face the camera
"""

import re
from typing import Optional, Dict, Any

class MoleculeVisualizer:
    """
    A class to convert and visualize molecular structures by:
      1. Converting SDF text to PDB data in-memory
      2. Generating minimal JavaScript for embedding in an existing Three.js scene
      3. Generating a full standalone HTML viewer
      4. Showing toggleable atomic annotations that always face the camera
    """

    enableAnnotations = True

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
            .atom-label {
                text-shadow: -1px 1px 1px rgb(0,0,0);
                margin-left: 5px;
                font-size: 14px;
                color: white;
                pointer-events: none;
            }
        """

    @staticmethod
    def _get_molecule_label_html(name: str) -> str:
        """
        Returns the HTML for the molecule label.
        Used by both standalone HTML and embedded JS versions.
        """
        return f'<div id="molecule-label">{name}</div>'

    

    @classmethod
    def generate_js_code_from_pdb(cls, pdb_data: str, name: str = "Molecule") -> str:
        """
        Return minimal JavaScript snippet that:
        - Converts SDF to PDB in Python, embed it here
        - Loads with PDBLoader (assumes THREE.PDBLoader is available)
        - Creates atoms & bonds, adds them to the scene
        - Adds toggleable atom annotations that face the camera
        - No custom cameras or lights; relies on the parent environment
        - Perfect for embedding in React or other Three.js setups
        """


        escaped_pdb = cls._escape_js_string(pdb_data)
        label_style = cls._escape_js_string(cls._get_molecule_label_style())
        label_html = cls._escape_js_string(cls._get_molecule_label_html(name))

        # Use double braces for JS blocks
        return f"""
function createMoleculeVisualization(THREE, scene, options = {{}}) {{
    console.log('createMoleculeVisualization');
    // Configuration options with defaults
    const config = {{
        enableAnnotations: {str(cls.enableAnnotations).lower()},  // Toggle atomic annotations
        scaleFactor: 0.25,       // Scale factor to control molecule size (reduced from 0.6)
        camera: null,           // Camera instance (optional)
        controls: null,         // Controls instance (optional)
        ...options
    }};
    
    // Create a group for the molecule
    const root = new THREE.Group();
    scene.add(root);
    
    // Store labels in a separate group for easier toggling
    const labelsGroup = new THREE.Group();
    root.add(labelsGroup);
    
    // Set a public property to allow external toggling of annotations
    root.enableAnnotations = config.enableAnnotations;

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

    // Set up CSS2D renderer for atom labels if it doesn't exist yet
    if (config.enableAnnotations && !window.labelRenderer && typeof THREE.CSS2DRenderer !== 'undefined') {{
        window.labelRenderer = new THREE.CSS2DRenderer();
        const updateLabelRendererSize = () => {{
            const container = document.querySelector('#container');
            if (container) {{
                const rect = container.getBoundingClientRect();
                window.labelRenderer.setSize(rect.width, rect.height);
            }}
        }};
        updateLabelRendererSize();
        window.labelRenderer.domElement.style.position = 'absolute';
        window.labelRenderer.domElement.style.top = '0px';
        window.labelRenderer.domElement.style.pointerEvents = 'none';
        const container = document.querySelector('#container');
        if (container) {{
            container.appendChild(window.labelRenderer.domElement);
        }} else {{
            document.body.appendChild(window.labelRenderer.domElement);
        }}
        
        // Add resize listener for labelRenderer if not already present
        if (!window.labelRendererResizeListener) {{
            window.labelRendererResizeListener = true;
            window.addEventListener('resize', updateLabelRendererSize);
        }}
    }}

    // Convert SDF -> PDB in Python, embed it here
    const pdbData = `{escaped_pdb}`;
    
    // Create and configure the PDB loader
    let loader;
    if (typeof THREE.PDBLoader !== 'undefined') {{
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

        // Add atoms and their labels
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
            object.position.multiplyScalar(1.5 * config.scaleFactor);
            object.scale.multiplyScalar(0.75 * config.scaleFactor);
            root.add(object);
            
            // Create atom annotation using CSS2DObject if available
            if (config.enableAnnotations && typeof THREE.CSS2DObject !== 'undefined') {{
                const atom = json.atoms[i];
                const atomSymbol = atom ? (atom[4] || '') : '';
                
                if (atomSymbol) {{
                    const text = document.createElement('div');
                    text.className = 'atom-label';
                    text.textContent = atomSymbol;
                    text.style.color = `rgb(${{Math.round(color.r*255)}},${{Math.round(color.g*255)}},${{Math.round(color.b*255)}})`;
                    
                    // Create CSS2DObject and attach it directly to the scene (not labelsGroup)
                    // This ensures it's not affected by group transformations
                    const label = new THREE.CSS2DObject(text);
                    label.position.copy(object.position);
                    scene.add(label);
                    
                    // Add reference to the label in the labelsGroup array for toggling
                    if (!labelsGroup.labels) labelsGroup.labels = [];
                    labelsGroup.labels.push(label);
                }}
            }}
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

            start.multiplyScalar(1.5 * config.scaleFactor);
            end.multiplyScalar(1.5 * config.scaleFactor);

            const object = new THREE.Mesh(
                boxGeometry,
                new THREE.MeshPhongMaterial({{ color: 0xffffff }})
            );
            object.position.copy(start);
            object.position.lerp(end, 0.5);
            object.scale.set(0.25 * config.scaleFactor, 0.25 * config.scaleFactor, start.distanceTo(end));
            object.lookAt(end);
            root.add(object);
        }}

        // Clean up
        URL.revokeObjectURL(pdbUrl);
        
        // Set initial visibility based on config
        labelsGroup.visible = config.enableAnnotations;

        // Fit camera to the molecule after loading
        root.fitCameraToMolecule();
    }});
    
    // Add a method to toggle annotations visibility
    root.toggleAnnotations = function(enable) {{
        if (typeof enable === 'boolean') {{
            root.enableAnnotations = enable;
        }} else {{
            root.enableAnnotations = !root.enableAnnotations;
        }}
        
        // Toggle visibility of each label in the labels array
        if (labelsGroup.labels && Array.isArray(labelsGroup.labels)) {{
            labelsGroup.labels.forEach(label => {{
                label.visible = root.enableAnnotations;
            }});
        }}
        
        return root.enableAnnotations;
    }};

    // Add method to fit camera to molecule
    root.fitCameraToMolecule = function() {{
        const box = new THREE.Box3();
        root.children.forEach(child => {{
            if (!(child instanceof THREE.Light)) {{
                box.expandByObject(child);
            }}
        }});
        
        const size = box.getSize(new THREE.Vector3());
        const center = box.getCenter(new THREE.Vector3());
        
        // Calculate distance based on diagonal
        const diagonal = Math.sqrt(
            size.x * size.x + 
            size.y * size.y + 
            size.z * size.z
        );
        
        // Increase distance for larger molecules using log scale
        const scaleFactor = Math.max(1.2, Math.log10(diagonal) * 0.8);
        const distance = diagonal * scaleFactor;
        
        // Position camera using spherical coordinates
        const theta = Math.PI / 4; // 45 degrees
        const phi = Math.PI / 6;   // 30 degrees
        
        config.camera.position.set(
            center.x + distance * Math.sin(theta) * Math.cos(phi),
            center.y + distance * Math.sin(phi),
            center.z + distance * Math.cos(theta) * Math.cos(phi)
        );
        
        config.camera.lookAt(center);
        config.controls.target.copy(center);
        
        // Adjust near/far planes
        config.camera.near = distance * 0.01;
        config.camera.far = distance * 10;
        config.camera.updateProjectionMatrix();
        
        // Update controls min/max distance
        config.controls.minDistance = distance * 0.1;
        config.controls.maxDistance = distance * 5;
        config.controls.update();
    }};
    
    // Return the root group for external control
    return root;
}}

// Function to make sure CSS2DRenderer is included in render loop
function setupAnnotationRenderer(renderer, scene, camera) {{
    if (!renderer || !scene || !camera) {{
        console.error('setupAnnotationRenderer requires renderer, scene, and camera parameters');
        return;
    }}
    const originalRender = renderer.render.bind(renderer);
    renderer.render = function(scene, camera) {{
        originalRender(scene, camera);
        if (window.labelRenderer) {{
            window.labelRenderer.render(scene, camera);
        }}
    }};
}}

createMoleculeVisualization(THREE, scene);
"""

    @classmethod
    def generate_html_viewer_from_pdb(cls, pdb_data: str, name: str = "Molecule") -> str:
        """
        Generate a standalone HTML file that:
        - Embeds Three.js (ESM)
        - Sets up the camera, lights, and controls
        - Converts SDF -> PDB in-memory, loads it, and displays the molecule
        - Provides a full animation loop with auto-rotation
        - Fills the browser window
        - Displays the molecule name
        - Includes toggleable atomic annotations
        """


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
        
        #controls {{
            position: absolute;
            top: 20px;
            right: 20px;
            background-color: rgba(0, 0, 0, 0.7);
            padding: 10px;
            border-radius: 5px;
            z-index: 100;
        }}
        
        #controls button {{
            background-color: #333;
            color: white;
            border: 1px solid #666;
            padding: 5px 10px;
            cursor: pointer;
            border-radius: 3px;
            transition: background-color 0.3s;
        }}
        
        #controls button:hover {{
            background-color: #555;
        }}
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
import {{ CSS2DRenderer, CSS2DObject }} from 'three/addons/renderers/CSS2DRenderer.js';

let camera, scene, renderer, labelRenderer, controls;
let root, labelsGroup;
const rotationSpeed = 0.0015; // Speed of auto-rotation (radians per frame)

// Configuration settings
const config = {{
    enableAnnotations: {str(cls.enableAnnotations).lower()}  // Default setting - can be toggled
}};

init();
animate();

function init() {{
    const container = document.getElementById('container');

    // Scene
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x050505);

    // Camera
    camera = new THREE.PerspectiveCamera(70, window.innerWidth / window.innerHeight, 1, 5000);
    camera.position.z = 500;

    // Lights
    const light1 = new THREE.DirectionalLight(0xffffff, 2.5);
    light1.position.set(1, 1, 1);
    scene.add(light1);

    const light2 = new THREE.DirectionalLight(0xffffff, 1.5);
    light2.position.set(-1, -1, 1);
    scene.add(light2);

    // Root group and labels group
    root = new THREE.Group();
    scene.add(root);

    // Renderer
    renderer = new THREE.WebGLRenderer({{ antialias: true }});
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth, window.innerHeight);
    container.appendChild(renderer.domElement);
    
    // CSS2D Renderer for atom labels
    labelRenderer = new CSS2DRenderer();
    labelRenderer.setSize(window.innerWidth, window.innerHeight);
    labelRenderer.domElement.style.position = 'absolute';
    labelRenderer.domElement.style.top = '0px';
    labelRenderer.domElement.style.pointerEvents = 'none';
    container.appendChild(labelRenderer.domElement);

    // Controls
    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.05;
    controls.autoRotate = true;
    controls.autoRotateSpeed = 1.5;

    // Handle resize
    window.addEventListener('resize', onWindowResize);
    
    // Add labelRenderer resize tracking
    if (!window.labelRendererResizeListener) {{
        window.labelRendererResizeListener = true;
        window.addEventListener('resize', () => {{
            if (window.labelRenderer) {{
                const container = document.querySelector('#container');
                if (container) {{
                    window.labelRenderer.setSize(container.clientWidth, container.clientHeight);
                }}
            }}
        }});
    }}

    // Add camera fitting function
    function fitCameraToMolecule() {{
        const box = new THREE.Box3();
        root.children.forEach(child => {{
            if (!(child instanceof THREE.Light)) {{
                box.expandByObject(child);
            }}
        }});
        
        const size = box.getSize(new THREE.Vector3());
        const center = box.getCenter(new THREE.Vector3());
        
        // Calculate distance based on diagonal
        const diagonal = Math.sqrt(
            size.x * size.x + 
            size.y * size.y + 
            size.z * size.z
        );
        
        // Increase distance for larger molecules using log scale
        const scaleFactor = Math.max(1.2, Math.log10(diagonal) * 0.8);
        const distance = diagonal * scaleFactor;
        
        // Position camera using spherical coordinates
        const theta = Math.PI / 4; // 45 degrees
        const phi = Math.PI / 6;   // 30 degrees
        
        camera.position.set(
            center.x + distance * Math.sin(theta) * Math.cos(phi),
            center.y + distance * Math.sin(phi),
            center.z + distance * Math.cos(theta) * Math.cos(phi)
        );
        
        camera.lookAt(center);
        controls.target.copy(center);
        
        // Adjust near/far planes
        camera.near = distance * 0.01;
        camera.far = distance * 10;
        camera.updateProjectionMatrix();
        
        // Update controls min/max distance
        controls.minDistance = distance * 0.1;
        controls.maxDistance = distance * 5;
        controls.update();
    }}

    // Load molecule
    const pdbData = `{escaped_pdb}`;
    const pdbBlob = new Blob([pdbData], {{ type: 'text/plain' }});
    const pdbUrl = URL.createObjectURL(pdbBlob);

    // Make PDBLoader available globally for compatibility with both approaches
    window.PDBLoader = PDBLoader;
    const loader = new PDBLoader();
    const scaleFactor = .3;
    
    // Array to track all labels
    const labels = [];
    
    // Add toggle function to root
    root.toggleAnnotations = function(enable) {{
        if (typeof enable === 'boolean') {{
            root.enableAnnotations = enable;
        }} else {{
            root.enableAnnotations = !root.enableAnnotations;
        }}
        
        // Toggle visibility of each label
        labels.forEach(label => {{
            label.visible = root.enableAnnotations;
        }});
        
        return root.enableAnnotations;
    }};
    
    // Set initial annotation state
    root.enableAnnotations = config.enableAnnotations;

    loader.load(pdbUrl, (pdb) => {{
        const geometryAtoms = pdb.geometryAtoms;
        const geometryBonds = pdb.geometryBonds;
        const json = pdb.json;
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

        // Add atoms and their labels
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
            atom.position.copy(position).multiplyScalar(120 * scaleFactor);
            atom.scale.setScalar(40 * scaleFactor);
            root.add(atom);
            
            // Add atom labels using CSS2DObject
            if (config.enableAnnotations && json.atoms[i]) {{
                const atomSymbol = json.atoms[i][4];
                if (atomSymbol) {{
                    const text = document.createElement('div');
                    text.className = 'atom-label';
                    text.textContent = atomSymbol;
                    text.style.color = `rgb(${{Math.round(color.r*255)}},${{Math.round(color.g*255)}},${{Math.round(color.b*255)}})`;
                    
                    // Create CSS2DObject and attach it directly to the scene
                    const label = new CSS2DObject(text);
                    label.position.copy(atom.position);
                    scene.add(label);
                    
                    // Add reference to the label for toggling
                    labels.push(label);
                }}
            }}
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
            ).multiplyScalar(120 * scaleFactor);

            end.set(
                positions.getX(i+1),
                positions.getY(i+1),
                positions.getZ(i+1)
            ).multiplyScalar(120 * scaleFactor);

            const bondMesh = new THREE.Mesh(
                boxGeometry,
                new THREE.MeshPhongMaterial({{ color: 0xffffff }})
            );
            bondMesh.position.copy(start).lerp(end, 0.5);
            bondMesh.scale.set(8 * scaleFactor, 8 * scaleFactor, start.distanceTo(end));
            bondMesh.lookAt(end);
            root.add(bondMesh);
        }}

        // Clean up
        URL.revokeObjectURL(pdbUrl);

        // Fit camera to the molecule after loading
        fitCameraToMolecule();
    }});
}}

function onWindowResize() {{
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
    labelRenderer.setSize(window.innerWidth, window.innerHeight);
}}

function animate() {{
    requestAnimationFrame(animate);
    
    // Update controls
    controls.update();
    
    // Render both the 3D scene and labels
    renderer.render(scene, camera);
    labelRenderer.render(scene, camera);
}}
</script>
</body>
</html>
"""

    @classmethod
    def generate_interactive_html(cls, 
                                 pdb_data: str, 
                                 title: Optional[str] = None,
                                 script_data: Optional[Dict[str, Any]] = None,
                                 output_path: Optional[str] = None) -> str:
        """
        Generate an interactive HTML visualization by injecting PDB data and optional script data
        into the output.html template file.
        
        Args:
            pdb_data (str): PDB data to inject into the template
            title (str, optional): Title for the molecule visualization
            script_data (dict, optional): Script data for interactive animation
            output_path (str, optional): Path where the generated HTML file will be saved.
                                        If not provided, returns the HTML content as a string.
        
        Returns:
            str: HTML content as a string if output_path is None, otherwise the path to the generated file
        """
        import os
        import re
        import json
        from pathlib import Path
        
        # Get the path to the template file (output.html in the same directory as this script)
        template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output.html')
        
        # Read the template file
        with open(template_path, 'r') as f:
            template_content = f.read()
        
        # Create default script data if none provided
        if script_data is None and title is not None:
            script_data = {
                "title": title,
                "content": [
                    {
                        "timecode": "00:00",
                        "atoms": [],
                        "caption": f"This is {title}."
                    }
                ]
            }
        
        # Replace script data if provided
        if script_data:
            # Convert script data to JSON string with proper formatting
            script_json = json.dumps(script_data, indent=2)
            
            # Find the scriptData constant in the template
            script_pattern = r'const scriptData = \{[\s\S]*?\};'
            script_match = re.search(script_pattern, template_content)
            
            if script_match:
                # Replace the scriptData constant
                template_content = template_content[:script_match.start()] + \
                                f'const scriptData = {script_json};' + \
                                template_content[script_match.end():]
        
        # Replace PDB data if provided
        if pdb_data:
            # Find the pdbData constant in the template
            pdb_pattern = r'const pdbData = `[\s\S]*?`;'
            pdb_match = re.search(pdb_pattern, template_content)
            
            if pdb_match:
                # Replace the pdbData constant
                template_content = template_content[:pdb_match.start()] + \
                                f'const pdbData = `{pdb_data}`;' + \
                                template_content[pdb_match.end():]
        
        # Write the modified content to the output file if path provided
        if output_path:
            with open(output_path, 'w') as f:
                f.write(template_content)
            print(f"Interactive visualization generated successfully: {output_path}")
            return output_path
        
        # Otherwise return the HTML content as a string
        return template_content

def main():
    """
    Simple CLI usage:
      python molecule_visualizer.py <input.sdf> [html|js]
    """
    import sys
    if len(sys.argv) < 2:
        print("Usage: python molecule_visualizer.py <input.sdf> [html|js]")
        sys.exit(1)

    pdb_file = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "html"

    with open(pdb_file, "r") as f:
        pdb_data = f.read()

    name = "Molecule"
    if mode == "js":
        code = MoleculeVisualizer.generate_js_code_from_pdb(pdb_data, name)
        print(code)
    else:
        html = MoleculeVisualizer.generate_html_viewer_from_pdb(pdb_data, name)
        print(html)

if __name__ == "__main__":
    main()