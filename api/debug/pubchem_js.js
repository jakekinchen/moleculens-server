
function createMoleculeVisualization(THREE, scene, options = {}) {
    console.log('createMoleculeVisualization');
    // Configuration options with defaults
    const config = {
        enableAnnotations: true,  // Toggle atomic annotations
        scaleFactor: 0.25,       // Scale factor to control molecule size (reduced from 0.6)
        camera: null,           // Camera instance (optional)
        controls: null,         // Controls instance (optional)
        ...options
    };
    
    // Create a group for the molecule
    const root = new THREE.Group();
    scene.add(root);
    
    // Store labels in a separate group for easier toggling
    const labelsGroup = new THREE.Group();
    root.add(labelsGroup);
    
    // Set a public property to allow external toggling of annotations
    root.enableAnnotations = config.enableAnnotations;

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
            .atom-label {
                text-shadow: -1px 1px 1px rgb(0,0,0);
                margin-left: 5px;
                font-size: 14px;
                color: white;
                pointer-events: none;
            }
        `;
        document.head.appendChild(style);
    }

    // Add molecule label if not already present
    if (!document.getElementById('molecule-label')) {
        const container = document.querySelector('#container');
        if (container) {
            const labelContainer = document.createElement('div');
            labelContainer.innerHTML = `<div id="molecule-label">[2.2.2]cryptand</div>`;
            container.appendChild(labelContainer.firstChild);
        }
    }

    // Set up CSS2D renderer for atom labels if it doesn't exist yet
    if (config.enableAnnotations && !window.labelRenderer && typeof THREE.CSS2DRenderer !== 'undefined') {
        window.labelRenderer = new THREE.CSS2DRenderer();
        const updateLabelRendererSize = () => {
            const container = document.querySelector('#container');
            if (container) {
                const rect = container.getBoundingClientRect();
                window.labelRenderer.setSize(rect.width, rect.height);
            }
        };
        updateLabelRendererSize();
        window.labelRenderer.domElement.style.position = 'absolute';
        window.labelRenderer.domElement.style.top = '0px';
        window.labelRenderer.domElement.style.pointerEvents = 'none';
        const container = document.querySelector('#container');
        if (container) {
            container.appendChild(window.labelRenderer.domElement);
        } else {
            document.body.appendChild(window.labelRenderer.domElement);
        }
        
        // Add resize listener for labelRenderer if not already present
        if (!window.labelRendererResizeListener) {
            window.labelRendererResizeListener = true;
            window.addEventListener('resize', updateLabelRendererSize);
        }
    }

    // Convert SDF -> PDB in Python, embed it here
    const pdbData = `COMPND    72801
HETATM    1  O1  UNL     1       7.006   0.546   0.000  1.00  0.00           O  
HETATM    2  O2  UNL     1       5.493  -2.009   0.000  1.00  0.00           O  
HETATM    3  O3  UNL     1      10.708   1.362   0.000  1.00  0.00           O  
HETATM    4  O4  UNL     1       2.754   2.381   0.000  1.00  0.00           O  
HETATM    5  O5  UNL     1       8.391  -3.475   0.000  1.00  0.00           O  
HETATM    6  O6  UNL     1       0.361  -1.557   0.000  1.00  0.00           O  
HETATM    7  N1  UNL     1       6.513   4.101   0.000  1.00  0.00           N  
HETATM    8  N2  UNL     1       3.931  -4.827   0.000  1.00  0.00           N  
HETATM    9  C1  UNL     1       6.264   2.739   0.000  1.00  0.00           C  
HETATM   10  C2  UNL     1       8.323   4.252   0.000  1.00  0.00           C  
HETATM   11  C3  UNL     1       4.857   4.747   0.000  1.00  0.00           C  
HETATM   12  C4  UNL     1       4.578  -3.434   0.000  1.00  0.00           C  
HETATM   13  C5  UNL     1       5.552  -5.761   0.000  1.00  0.00           C  
HETATM   14  C6  UNL     1       1.975  -4.543   0.000  1.00  0.00           C  
HETATM   15  C7  UNL     1       7.249   1.780   0.000  1.00  0.00           C  
HETATM   16  C8  UNL     1       4.216  -2.042   0.000  1.00  0.00           C  
HETATM   17  C9  UNL     1       9.802   2.909   0.000  1.00  0.00           C  
HETATM   18  C10 UNL     1       3.053   3.986   0.000  1.00  0.00           C  
HETATM   19  C11 UNL     1       7.479  -5.075   0.000  1.00  0.00           C  
HETATM   20  C12 UNL     1       0.504  -3.269   0.000  1.00  0.00           C  
HETATM   21  C13 UNL     1       5.821   0.135   0.000  1.00  0.00           C  
HETATM   22  C14 UNL     1       6.429  -1.158   0.000  1.00  0.00           C  
HETATM   23  C15 UNL     1      10.381  -0.445   0.000  1.00  0.00           C  
HETATM   24  C16 UNL     1       1.387   1.520   0.000  1.00  0.00           C  
HETATM   25  C17 UNL     1       9.759  -2.317   0.000  1.00  0.00           C  
HETATM   26  C18 UNL     1       1.347  -0.308   0.000  1.00  0.00           C  
HETATM   27  H1  UNL     1       5.234   2.664   0.000  1.00  0.00           H  
HETATM   28  H2  UNL     1       5.570   1.997   0.000  1.00  0.00           H  
HETATM   29  H3  UNL     1       8.068   5.327   0.000  1.00  0.00           H  
HETATM   30  H4  UNL     1       9.147   5.023   0.000  1.00  0.00           H  
HETATM   31  H5  UNL     1       4.688   5.804   0.000  1.00  0.00           H  
HETATM   32  H6  UNL     1       4.497   3.798   0.000  1.00  0.00           H  
HETATM   33  H7  UNL     1       3.695  -3.630   0.000  1.00  0.00           H  
HETATM   34  H8  UNL     1       5.578  -3.555   0.000  1.00  0.00           H  
HETATM   35  H9  UNL     1       5.786  -4.852   0.000  1.00  0.00           H  
HETATM   36  H10 UNL     1       5.587  -6.822   0.000  1.00  0.00           H  
HETATM   37  H11 UNL     1       1.831  -3.418   0.000  1.00  0.00           H  
HETATM   38  H12 UNL     1       2.440  -3.670   0.000  1.00  0.00           H  
HETATM   39  H13 UNL     1       8.156   1.289   0.000  1.00  0.00           H  
HETATM   40  H14 UNL     1       7.926   2.488   0.000  1.00  0.00           H  
HETATM   41  H15 UNL     1       3.192  -2.230   0.000  1.00  0.00           H  
HETATM   42  H16 UNL     1       3.705  -1.119   0.000  1.00  0.00           H  
HETATM   43  H17 UNL     1      10.724   3.454   0.000  1.00  0.00           H  
HETATM   44  H18 UNL     1       9.351   2.030   0.000  1.00  0.00           H  
HETATM   45  H19 UNL     1       2.650   5.049   0.000  1.00  0.00           H  
HETATM   46  H20 UNL     1       1.963   4.159   0.000  1.00  0.00           H  
HETATM   47  H21 UNL     1       7.071  -4.107   0.000  1.00  0.00           H  
HETATM   48  H22 UNL     1       8.296  -5.781   0.000  1.00  0.00           H  
HETATM   49  H23 UNL     1      -0.608  -3.195   0.000  1.00  0.00           H  
HETATM   50  H24 UNL     1      -0.139  -4.198   0.000  1.00  0.00           H  
HETATM   51  H25 UNL     1       4.890   0.705   0.000  1.00  0.00           H  
HETATM   52  H26 UNL     1       4.959  -0.361   0.000  1.00  0.00           H  
HETATM   53  H27 UNL     1       7.038  -2.031   0.000  1.00  0.00           H  
HETATM   54  H28 UNL     1       7.505  -1.056   0.000  1.00  0.00           H  
HETATM   55  H29 UNL     1      11.442  -0.550   0.000  1.00  0.00           H  
HETATM   56  H30 UNL     1       9.392  -0.067   0.000  1.00  0.00           H  
HETATM   57  H31 UNL     1       0.831   2.471   0.000  1.00  0.00           H  
HETATM   58  H32 UNL     1       0.273   1.492   0.000  1.00  0.00           H  
HETATM   59  H33 UNL     1      10.531  -3.061   0.000  1.00  0.00           H  
HETATM   60  H34 UNL     1       8.882  -1.805   0.000  1.00  0.00           H  
HETATM   61  H35 UNL     1       2.009  -1.114   0.000  1.00  0.00           H  
HETATM   62  H36 UNL     1       2.425  -0.092   0.000  1.00  0.00           H  
CONECT    1   15   21
CONECT    2   16   22
CONECT    3   17   23
CONECT    4   18   24
CONECT    5   19   25
CONECT    6   20   26
CONECT    7    9   10   11
CONECT    8   12   13   14
CONECT    9   15   27   28
CONECT   10   17   29   30
CONECT   11   18   31   32
CONECT   12   16   33   34
CONECT   13   19   35   36
CONECT   14   20   37   38
CONECT   15   39   40
CONECT   16   41   42
CONECT   17   43   44
CONECT   18   45   46
CONECT   19   47   48
CONECT   20   49   50
CONECT   21   22   51   52
CONECT   22   53   54
CONECT   23   25   55   56
CONECT   24   26   57   58
CONECT   25   59   60
CONECT   26   61   62
END
`;
    
    // Create and configure the PDB loader
    let loader;
    if (typeof THREE.PDBLoader !== 'undefined') {
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

        // Add atoms and their labels
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
            object.position.multiplyScalar(1.5 * config.scaleFactor);
            object.scale.multiplyScalar(0.75 * config.scaleFactor);
            root.add(object);
            
            // Create atom annotation using CSS2DObject if available
            if (config.enableAnnotations && typeof THREE.CSS2DObject !== 'undefined') {
                const atom = json.atoms[i];
                const atomSymbol = atom ? (atom[4] || '') : '';
                
                if (atomSymbol) {
                    const text = document.createElement('div');
                    text.className = 'atom-label';
                    text.textContent = atomSymbol;
                    text.style.color = `rgb(${Math.round(color.r*255)},${Math.round(color.g*255)},${Math.round(color.b*255)})`;
                    
                    // Create CSS2DObject and attach it directly to the scene (not labelsGroup)
                    // This ensures it's not affected by group transformations
                    const label = new THREE.CSS2DObject(text);
                    label.position.copy(object.position);
                    scene.add(label);
                    
                    // Add reference to the label in the labelsGroup array for toggling
                    if (!labelsGroup.labels) labelsGroup.labels = [];
                    labelsGroup.labels.push(label);
                }
            }
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

            start.multiplyScalar(1.5 * config.scaleFactor);
            end.multiplyScalar(1.5 * config.scaleFactor);

            const object = new THREE.Mesh(
                boxGeometry,
                new THREE.MeshPhongMaterial({ color: 0xffffff })
            );
            object.position.copy(start);
            object.position.lerp(end, 0.5);
            object.scale.set(0.25 * config.scaleFactor, 0.25 * config.scaleFactor, start.distanceTo(end));
            object.lookAt(end);
            root.add(object);
        }

        // Clean up
        URL.revokeObjectURL(pdbUrl);
        
        // Set initial visibility based on config
        labelsGroup.visible = config.enableAnnotations;

        // Fit camera to the molecule after loading
        root.fitCameraToMolecule();
    });
    
    // Add a method to toggle annotations visibility
    root.toggleAnnotations = function(enable) {
        if (typeof enable === 'boolean') {
            root.enableAnnotations = enable;
        } else {
            root.enableAnnotations = !root.enableAnnotations;
        }
        
        // Toggle visibility of each label in the labels array
        if (labelsGroup.labels && Array.isArray(labelsGroup.labels)) {
            labelsGroup.labels.forEach(label => {
                label.visible = root.enableAnnotations;
            });
        }
        
        return root.enableAnnotations;
    };

    // Add method to fit camera to molecule
    root.fitCameraToMolecule = function() {
        const box = new THREE.Box3();
        root.children.forEach(child => {
            if (!(child instanceof THREE.Light)) {
                box.expandByObject(child);
            }
        });
        
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
    };
    
    // Return the root group for external control
    return root;
}

// Function to make sure CSS2DRenderer is included in render loop
function setupAnnotationRenderer(renderer, scene, camera) {
    if (!renderer || !scene || !camera) {
        console.error('setupAnnotationRenderer requires renderer, scene, and camera parameters');
        return;
    }
    const originalRender = renderer.render.bind(renderer);
    renderer.render = function(scene, camera) {
        originalRender(scene, camera);
        if (window.labelRenderer) {
            window.labelRenderer.render(scene, camera);
        }
    };
}

createMoleculeVisualization(THREE, scene);
