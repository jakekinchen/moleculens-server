"""
Geometry Agent - Provides static geometry code snippets for the visualization.
"""

import json

from agent_management.llm_service import LLMRequest, LLMService
from agent_management.utils.code_extraction import extract_code_block


class GeometryAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service

    def get_geometry_snippet(self, user_prompt: str) -> str:
        """
        Generate Three.js geometry code using LLM.
        """
        prompt_for_llm = f"""You are a helpful assistant generating Three.js geometry code. The user wants geometry for a scene with the following prompt: '{user_prompt}'

Requirements:
1. Code Style and Syntax:
   - Use proper semicolons after ALL complete statements, including:
     * Variable declarations and assignments
     * Method calls
     * Multi-line method chains (after the last method)
     * Function expressions
   - Use consistent 4-space indentation
   - Add descriptive comments for each major section
   - Use proper variable declarations (const/let)

2. Three.js Specifics:
   - Return only JavaScript code that references 'scene' (a global variable)
   - Do NOT create or re-declare lights or cameras
   - Do NOT include text sprites or labels
   - Use window.* references for global objects (e.g., window.myMesh = ...)
   - For molecular structures, group atoms and bonds together first
   - Use THREE.Curve as an ES6 class with 'new' keyword

3. Code Organization:
   - Create materials first
   - Then create geometries
   - Group related objects together
   - Add cleanup/disposal functions when needed
   - Add error handling for undefined scene object

Here's an example of proper code style and semicolon usage, particularly for molecular structures:

```javascript
// Verify scene exists
if (typeof scene === 'undefined') {{
    throw new Error('Scene object is undefined.');
}}

// Create materials with proper semicolons
const carbonMaterial = new THREE.MeshPhongMaterial({{ color: 0x333333 }});
const hydrogenMaterial = new THREE.MeshPhongMaterial({{ color: 0xffffff }});

// Create geometries
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create a group for the molecule
const molecule = new THREE.Group();
window.molecule = molecule;  // Store global reference

// Create atoms with proper multi-line statements
const carbon = new THREE.Mesh(atomGeometry, carbonMaterial);
carbon.position.copy(
    new THREE.Vector3(0, 0, 0)
).normalize();  // Note the semicolon after multi-line method chain

const hydrogen = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen.position.set(
    Math.cos(109.5 * Math.PI / 180) * bondLength,
    Math.sin(109.5 * Math.PI / 180) * bondLength,
    0
);  // Note the semicolon after multi-line method call

// Bond creation function with proper returns
function createBond(start, end) {{
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();

    const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
    const bond = new THREE.Mesh(bondGeometry, bondMaterial);

    bond.position.copy(start).lerp(end, 0.5);  // Method chaining with semicolon

    bond.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );  // Multi-line method with semicolon

    return bond;
}}

// Add to scene with proper cleanup
molecule.add(carbon, hydrogen);
scene.add(molecule);

// Add cleanup function
if (typeof window !== 'undefined') {{
    window.disposeMolecule = function() {{
        if (window.molecule) {{
            window.molecule.traverse((object) => {{
                if (object.geometry) object.geometry.dispose();
                if (object.material) {{
                    if (Array.isArray(object.material)) {{
                        object.material.forEach(m => m.dispose());
                    }} else {{
                        object.material.dispose();
                    }}
                }}
            }});
            scene.remove(window.molecule);
            window.molecule = undefined;
        }}
    }};
}}
```

Return your code as a single JavaScript snippet, with no JSON wrapper.
"""

        # Create a proper LLMRequest with higher max_tokens
        request = LLMRequest(user_prompt=prompt_for_llm)

        # Get the response from the LLM service
        llm_response = self.llm_service.generate(request)

        # Extract the content from the LLMResponse and clean it
        geometry_code = extract_code_block(llm_response.content, "javascript")

        return f"""
// GeometryAgent LLM-generated code
{geometry_code}
"""
