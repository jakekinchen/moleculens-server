
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
            labelContainer.innerHTML = `<div id="molecule-label">Tetraphenylporphyrin</div>`;
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
    const pdbData = `COMPND    86280046
HETATM    1  N1  UNL     1       4.604   1.473   0.000  1.00  0.00           N  
HETATM    2  N2  UNL     1       7.909  -1.812   0.000  1.00  0.00           N  
HETATM    3  N3  UNL     1       4.442  -1.405   0.000  1.00  0.00           N  
HETATM    4  N4  UNL     1       7.758   1.092   0.000  1.00  0.00           N  
HETATM    5  C1  UNL     1       3.176   1.372   0.000  1.00  0.00           C  
HETATM    6  C2  UNL     1       4.981   2.790   0.000  1.00  0.00           C  
HETATM    7  C3  UNL     1       7.438  -3.090   0.000  1.00  0.00           C  
HETATM    8  C4  UNL     1       9.253  -1.690   0.000  1.00  0.00           C  
HETATM    9  C5  UNL     1       2.270  -0.014   0.000  1.00  0.00           C  
HETATM   10  C6  UNL     1       6.469   3.440   0.000  1.00  0.00           C  
HETATM   11  C7  UNL     1       5.791  -3.730   0.000  1.00  0.00           C  
HETATM   12  C8  UNL     1      10.100  -0.283   0.000  1.00  0.00           C  
HETATM   13  C9  UNL     1       3.102  -1.127   0.000  1.00  0.00           C  
HETATM   14  C10 UNL     1       7.448   2.496   0.000  1.00  0.00           C  
HETATM   15  C11 UNL     1       4.707  -2.787   0.000  1.00  0.00           C  
HETATM   16  C12 UNL     1       9.221   0.792   0.000  1.00  0.00           C  
HETATM   17  C13 UNL     1       2.736   2.637   0.000  1.00  0.00           C  
HETATM   18  C14 UNL     1       3.841   3.492   0.000  1.00  0.00           C  
HETATM   19  C15 UNL     1       8.652  -3.718   0.000  1.00  0.00           C  
HETATM   20  C16 UNL     1       9.781  -2.909   0.000  1.00  0.00           C  
HETATM   21  C17 UNL     1       0.785  -0.073   0.000  1.00  0.00           C  
HETATM   22  C18 UNL     1       6.681   5.136   0.000  1.00  0.00           C  
HETATM   23  C19 UNL     1       5.593  -5.300   0.000  1.00  0.00           C  
HETATM   24  C20 UNL     1      11.870  -0.326   0.000  1.00  0.00           C  
HETATM   25  C21 UNL     1       2.452  -2.369   0.000  1.00  0.00           C  
HETATM   26  C22 UNL     1       8.720   3.055   0.000  1.00  0.00           C  
HETATM   27  C23 UNL     1       3.378  -3.310   0.000  1.00  0.00           C  
HETATM   28  C24 UNL     1       9.626   2.105   0.000  1.00  0.00           C  
HETATM   29  C25 UNL     1       0.003   1.199   0.000  1.00  0.00           C  
HETATM   30  C26 UNL     1       7.847   5.891   0.000  1.00  0.00           C  
HETATM   31  C27 UNL     1      -0.167  -1.271   0.000  1.00  0.00           C  
HETATM   32  C28 UNL     1       5.630   6.127   0.000  1.00  0.00           C  
HETATM   33  C29 UNL     1       4.490  -6.072   0.000  1.00  0.00           C  
HETATM   34  C30 UNL     1      12.953  -1.105   0.000  1.00  0.00           C  
HETATM   35  C31 UNL     1       6.749  -6.135   0.000  1.00  0.00           C  
HETATM   36  C32 UNL     1      12.341   1.015   0.000  1.00  0.00           C  
HETATM   37  C33 UNL     1      -1.515   1.196   0.000  1.00  0.00           C  
HETATM   38  C34 UNL     1       8.093   7.266   0.000  1.00  0.00           C  
HETATM   39  C35 UNL     1      -1.688  -1.304   0.000  1.00  0.00           C  
HETATM   40  C36 UNL     1       5.760   7.543   0.000  1.00  0.00           C  
HETATM   41  C37 UNL     1       4.499  -7.477   0.000  1.00  0.00           C  
HETATM   42  C38 UNL     1      14.371  -0.402   0.000  1.00  0.00           C  
HETATM   43  C39 UNL     1       6.782  -7.595   0.000  1.00  0.00           C  
HETATM   44  C40 UNL     1      13.498   1.770   0.000  1.00  0.00           C  
HETATM   45  C41 UNL     1      -2.337  -0.021   0.000  1.00  0.00           C  
HETATM   46  C42 UNL     1       6.984   8.101   0.000  1.00  0.00           C  
HETATM   47  C43 UNL     1       5.668  -8.290   0.000  1.00  0.00           C  
HETATM   48  C44 UNL     1      14.630   1.022   0.000  1.00  0.00           C  
HETATM   49  H1  UNL     1       5.240   0.692   0.000  1.00  0.00           H  
HETATM   50  H2  UNL     1       7.346  -0.975   0.000  1.00  0.00           H  
HETATM   51  H3  UNL     1       1.793   3.139   0.000  1.00  0.00           H  
HETATM   52  H4  UNL     1       3.667   4.513   0.000  1.00  0.00           H  
HETATM   53  H5  UNL     1       8.926  -4.768   0.000  1.00  0.00           H  
HETATM   54  H6  UNL     1      10.839  -2.906   0.000  1.00  0.00           H  
HETATM   55  H7  UNL     1       1.474  -2.769   0.000  1.00  0.00           H  
HETATM   56  H8  UNL     1       9.241   3.966   0.000  1.00  0.00           H  
HETATM   57  H9  UNL     1       3.007  -4.325   0.000  1.00  0.00           H  
HETATM   58  H10 UNL     1      10.534   2.703   0.000  1.00  0.00           H  
HETATM   59  H11 UNL     1       0.493   2.114   0.000  1.00  0.00           H  
HETATM   60  H12 UNL     1       8.838   5.466   0.000  1.00  0.00           H  
HETATM   61  H13 UNL     1      -0.679  -0.443   0.000  1.00  0.00           H  
HETATM   62  H14 UNL     1       4.622   5.820   0.000  1.00  0.00           H  
HETATM   63  H15 UNL     1       3.458  -5.793   0.000  1.00  0.00           H  
HETATM   64  H16 UNL     1      13.142  -2.119   0.000  1.00  0.00           H  
HETATM   65  H17 UNL     1       7.711  -5.664   0.000  1.00  0.00           H  
HETATM   66  H18 UNL     1      11.456   1.602   0.000  1.00  0.00           H  
HETATM   67  H19 UNL     1      -2.093   2.133   0.000  1.00  0.00           H  
HETATM   68  H20 UNL     1       9.072   7.538   0.000  1.00  0.00           H  
HETATM   69  H21 UNL     1      -2.159  -2.253   0.000  1.00  0.00           H  
HETATM   70  H22 UNL     1       4.851   8.155   0.000  1.00  0.00           H  
HETATM   71  H23 UNL     1       3.495  -7.914   0.000  1.00  0.00           H  
HETATM   72  H24 UNL     1      15.189  -1.055   0.000  1.00  0.00           H  
HETATM   73  H25 UNL     1       7.729  -8.111   0.000  1.00  0.00           H  
HETATM   74  H26 UNL     1      13.471   2.898   0.000  1.00  0.00           H  
HETATM   75  H27 UNL     1      -3.358   0.056   0.000  1.00  0.00           H  
HETATM   76  H28 UNL     1       7.168   9.167   0.000  1.00  0.00           H  
HETATM   77  H29 UNL     1       5.597  -9.419   0.000  1.00  0.00           H  
HETATM   78  H30 UNL     1      15.569   1.455   0.000  1.00  0.00           H  
CONECT    1    5    6   49
CONECT    2    7    8   50
CONECT    3   13   13   15
CONECT    4   14   14   16
CONECT    5    9    9   17
CONECT    6   10   10   18
CONECT    7   11   19   19
CONECT    8   12   20   20
CONECT    9   13   21
CONECT   10   14   22
CONECT   11   15   15   23
CONECT   12   16   16   24
CONECT   13   25
CONECT   14   26
CONECT   15   27
CONECT   16   28
CONECT   17   18   18   51
CONECT   18   52
CONECT   19   20   53
CONECT   20   54
CONECT   21   29   29   31
CONECT   22   30   30   32
CONECT   23   33   33   35
CONECT   24   34   34   36
CONECT   25   27   27   55
CONECT   26   28   28   56
CONECT   27   57
CONECT   28   58
CONECT   29   37   59
CONECT   30   38   60
CONECT   31   39   39   61
CONECT   32   40   40   62
CONECT   33   41   63
CONECT   34   42   64
CONECT   35   43   43   65
CONECT   36   44   44   66
CONECT   37   45   45   67
CONECT   38   46   46   68
CONECT   39   45   69
CONECT   40   46   70
CONECT   41   47   47   71
CONECT   42   48   48   72
CONECT   43   47   73
CONECT   44   48   74
CONECT   45   75
CONECT   46   76
CONECT   47   77
CONECT   48   78
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
