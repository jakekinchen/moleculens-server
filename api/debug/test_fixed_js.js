// GeometryAgent LLM-generated code
function createMoleculeVisualization(THREE, scene) {
    // Create a group for the molecule
    const root = new THREE.Group();
    scene.add(root);

    // Convert SDF -> PDB in Python, embed it here
    const pdbData = `COMPND    702
HETATM    1  O1  UNL     1       4.124   0.452   0.000  1.00  0.00           O  
HETATM    2  C1  UNL     1       3.022  -0.390   0.000  1.00  0.00           C  
HETATM    3  C2  UNL     1       1.733   0.337   0.000  1.00  0.00           C  
HETATM    4  H1  UNL     1       2.516  -1.387   0.000  1.00  0.00           H  
HETATM    5  H2  UNL     1       3.667  -1.286   0.000  1.00  0.00           H  
HETATM    6  H3  UNL     1       2.188   1.333   0.000  1.00  0.00           H  
HETATM    7  H4  UNL     1       0.968   1.138   0.000  1.00  0.00           H  
HETATM    8  H5  UNL     1       0.887  -0.347   0.000  1.00  0.00           H  
HETATM    9  H6  UNL     1       4.957  -0.049   0.000  1.00  0.00           H  
CONECT    1    2    9
CONECT    2    3    4    5
CONECT    3    6    7    8
END`;
    
    // Create and configure the PDB loader
    // First check if it's available globally (for React/embedded environments)
    let loader;
    if (typeof window !== 'undefined' && window.PDBLoader) {
        loader = new window.PDBLoader();
    } else if (typeof PDBLoader !== 'undefined') {
        // For module environments where PDBLoader is imported directly
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
            object.position.multiplyScalar(75);
            object.scale.multiplyScalar(25);
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

            start.multiplyScalar(75);
            end.multiplyScalar(75);

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