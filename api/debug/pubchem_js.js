// GeometryAgent LLM-generated code
function createMoleculeVisualization(THREE, scene) {
    // Create a group for the molecule
    const root = new THREE.Group();
    scene.add(root);

    // Add molecule label styles if not already present
    if (!document.getElementById('molecule-label-style')) {
        const style = document.createElement('style');
        style.id = 'molecule-label-style';
        style.textContent = `
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
        `;
        document.head.appendChild(style);
    }

    // Add molecule label if not already present
    if (!document.getElementById('molecule-label')) {
        const container = document.querySelector('#container');
        if (container) {
            const labelContainer = document.createElement('div');
            labelContainer.innerHTML = `<div id="molecule-label">methane</div>`;
            container.appendChild(labelContainer.firstChild);
        }
    }

    // Convert SDF -> PDB in Python, embed it here
    const pdbData = `COMPND    297
HETATM    1  C1  UNL     1       2.514  -0.093   0.000  1.00  0.00           C  
HETATM    2  H1  UNL     1       3.363   0.595   0.000  1.00  0.00           H  
HETATM    3  H2  UNL     1       1.512  -0.496   0.000  1.00  0.00           H  
HETATM    4  H3  UNL     1       2.086   0.913   0.000  1.00  0.00           H  
HETATM    5  H4  UNL     1       3.209  -0.919   0.000  1.00  0.00           H  
CONECT    1    2    3    4    5
END
`;
    
    // Create and configure the PDB loader
    let loader;
    if (typeof THREE !== 'undefined' && THREE.PDBLoader) {
        loader = new THREE.PDBLoader();
    } else if (typeof window !== 'undefined' && window.PDBLoader) {
        // If we manually attached PDBLoader to the window
        loader = new window.PDBLoader();
    } else if (typeof PDBLoader !== 'undefined') {
        loader = new PDBLoader();
    } else {
        console.error('PDBLoader not found. Make sure it is loaded first.');
        return root;
    }
    
    const pdbBlob = new Blob([pdbData], { type: 'text/plain' });
    const pdbUrl = URL.createObjectURL(pdbBlob);

    // Load and process the PDB data
    loader.load(pdbUrl, function (pdb) {
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
        for (let i = 0; i < positions.count; i++) {
            position.x = positions.getX(i);
            position.y = positions.getY(i);
            position.z = positions.getZ(i);

            color.r = colors.getX(i);
            color.g = colors.getY(i);
            color.b = colors.getZ(i);

            const material = new THREE.MeshPhongMaterial({ color: color });
            const object = new THREE.Mesh(sphereGeometry, material);
            object.position.copy(position);
            object.position.multiplyScalar(20);
            object.scale.multiplyScalar(7.5);
            root.add(object);
        }

        // Add bonds
        positions = geometryBonds.getAttribute('position');
        const start = new THREE.Vector3();
        const end = new THREE.Vector3();

        for (let i = 0; i < positions.count; i += 2) {
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
                new THREE.MeshPhongMaterial({ color: 0xffffff })
            );
            object.position.copy(start);
            object.position.lerp(end, 0.5);
            object.scale.set(5, 5, start.distanceTo(end));
            object.lookAt(end);
            root.add(object);
        }

        // Clean up
        URL.revokeObjectURL(pdbUrl);
    });

    // Return the root group for external control
    return root;
}

// Execute the function to create the visualization
createMoleculeVisualization(THREE, scene);
