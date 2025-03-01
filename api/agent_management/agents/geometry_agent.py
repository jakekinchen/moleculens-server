"""
Geometry Agent - Provides static geometry code snippets for the visualization.
"""

import json
from llm_service import LLMService

class GeometryAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service

    def get_geometry_snippet(self, user_prompt: str) -> str:
        """
        Generate Three.js geometry code using LLM.
        """
        prompt_for_llm = f"You are a helpful assistant generating Three.js geometry code. The user wants geometry for a scene with the following prompt: '{user_prompt}'\n\n" + r"""
Requirements:
- Return only JavaScript code that references 'scene' (a global variable).
- Do NOT create or re-declare lights or cameras. Only geometry or materials.
- Use window.* references if you need to store global references, e.g., window.myMesh = ...
- For molecular structures, group atoms and bonds together first and then add the group to the scene.

Here's an example of how to properly create an ethanol molecule when prompted with "make an ethanol molecule":

```javascript
// Create atom materials
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x333333 });
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999 });

// Create atom geometry
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create a group for the molecule first
const molecule = new THREE.Group();
window.molecule = molecule;

// Carbon atoms
const c1 = new THREE.Mesh(atomGeometry, carbonMaterial);
c1.position.set(0, 0, 0);
c1.scale.set(0.5, 0.5, 0.5);

const c2 = new THREE.Mesh(atomGeometry, carbonMaterial);
c2.position.set(1.5, 0, 0);
c2.scale.set(0.5, 0.5, 0.5);

// Oxygen atom
const o = new THREE.Mesh(atomGeometry, oxygenMaterial);
o.position.set(3, 0, 0);
o.scale.set(0.55, 0.55, 0.55);

// Hydrogen atoms
const h1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h1.position.set(-0.8, 0.8, 0);
h1.scale.set(0.3, 0.3, 0.3);

const h2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h2.position.set(-0.8, -0.8, 0);
h2.scale.set(0.3, 0.3, 0.3);

const h3 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h3.position.set(0, 0, -0.8);
h3.scale.set(0.3, 0.3, 0.3);

const h4 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h4.position.set(1.5, 0.8, 0);
h4.scale.set(0.3, 0.3, 0.3);

const h5 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h5.position.set(1.5, -0.8, 0);
h5.scale.set(0.3, 0.3, 0.3);

const h6 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h6.position.set(3.8, 0, 0);
h6.scale.set(0.3, 0.3, 0.3);

// Create bonds function
function createBond(start, end) {
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();
    
    const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
    const bond = new THREE.Mesh(bondGeometry, bondMaterial);
    
    // Position and rotate the bond
    bond.position.copy(start);
    bond.position.lerp(end, 0.5);
    
    // Orient the cylinder
    bond.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );
    
    return bond;
}

// Create all bonds
const bonds = [
    createBond(c1.position, c2.position),
    createBond(c2.position, o.position),
    createBond(c1.position, h1.position),
    createBond(c1.position, h2.position),
    createBond(c1.position, h3.position),
    createBond(c2.position, h4.position),
    createBond(c2.position, h5.position),
    createBond(o.position, h6.position)
];

// Add all atoms and bonds to the molecule group
molecule.add(c1, c2, o, h1, h2, h3, h4, h5, h6, ...bonds);

// Then add the molecule group to the scene
scene.add(molecule);
```

Return your code as a single JavaScript snippet, with no JSON wrapper.
"""

        llm_response = self.llm_service.generate(prompt_for_llm)
        geometry_code = llm_response.strip()
        
        return f"""
// GeometryAgent LLM-generated code
{geometry_code}
"""