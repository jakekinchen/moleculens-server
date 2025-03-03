
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
            labelContainer.innerHTML = `<div id="molecule-label">(C60-Ih)[5,6]fullerene</div>`;
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
    const pdbData = `COMPND    123591
HETATM    1  C1  UNL     1       2.677  -0.465   0.000  1.00  0.00           C  
HETATM    2  C2  UNL     1       2.053   1.120   0.000  1.00  0.00           C  
HETATM    3  C3  UNL     1       8.968  -2.043   0.000  1.00  0.00           C  
HETATM    4  C4  UNL     1       6.852  -1.754   0.000  1.00  0.00           C  
HETATM    5  C5  UNL     1       2.445   3.363   0.000  1.00  0.00           C  
HETATM    6  C6  UNL     1       7.468   4.674   0.000  1.00  0.00           C  
HETATM    7  C7  UNL     1       4.697  -2.840   0.000  1.00  0.00           C  
HETATM    8  C8  UNL     1       7.448  -4.129   0.000  1.00  0.00           C  
HETATM    9  C9  UNL     1       0.589   2.878   0.000  1.00  0.00           C  
HETATM   10  C10 UNL     1       7.339  -0.968   0.000  1.00  0.00           C  
HETATM   11  C11 UNL     1       7.059   1.009   0.000  1.00  0.00           C  
HETATM   12  C12 UNL     1       8.675  -0.029   0.000  1.00  0.00           C  
HETATM   13  C13 UNL     1       4.715  -1.397   0.000  1.00  0.00           C  
HETATM   14  C14 UNL     1       4.911   0.392   0.000  1.00  0.00           C  
HETATM   15  C15 UNL     1       6.648   3.088   0.000  1.00  0.00           C  
HETATM   16  C16 UNL     1       1.854   0.142   0.000  1.00  0.00           C  
HETATM   17  C17 UNL     1       5.616  -2.383   0.000  1.00  0.00           C  
HETATM   18  C18 UNL     1       2.063  -4.122   0.000  1.00  0.00           C  
HETATM   19  C19 UNL     1       6.145  -3.981   0.000  1.00  0.00           C  
HETATM   20  C20 UNL     1       8.577   4.049   0.000  1.00  0.00           C  
HETATM   21  C21 UNL     1       5.875  -5.150   0.000  1.00  0.00           C  
HETATM   22  C22 UNL     1       7.054   2.385   0.000  1.00  0.00           C  
HETATM   23  C23 UNL     1       9.341   1.786   0.000  1.00  0.00           C  
HETATM   24  C24 UNL     1       6.243  -0.205   0.000  1.00  0.00           C  
HETATM   25  C25 UNL     1       7.086  -3.109   0.000  1.00  0.00           C  
HETATM   26  C26 UNL     1       7.201  -4.714   0.000  1.00  0.00           C  
HETATM   27  C27 UNL     1       4.602  -5.685   0.000  1.00  0.00           C  
HETATM   28  C28 UNL     1       3.790  -4.737   0.000  1.00  0.00           C  
HETATM   29  C29 UNL     1       3.774   3.854   0.000  1.00  0.00           C  
HETATM   30  C30 UNL     1       5.591   1.626   0.000  1.00  0.00           C  
HETATM   31  C31 UNL     1       8.049   3.492   0.000  1.00  0.00           C  
HETATM   32  C32 UNL     1       3.267  -3.205   0.000  1.00  0.00           C  
HETATM   33  C33 UNL     1       8.845  -2.621   0.000  1.00  0.00           C  
HETATM   34  C34 UNL     1       1.729  -2.309   0.000  1.00  0.00           C  
HETATM   35  C35 UNL     1       0.480  -0.822   0.000  1.00  0.00           C  
HETATM   36  C36 UNL     1       0.670  -0.504   0.000  1.00  0.00           C  
HETATM   37  C37 UNL     1       4.520   2.469   0.000  1.00  0.00           C  
HETATM   38  C38 UNL     1       5.585   5.101   0.000  1.00  0.00           C  
HETATM   39  C39 UNL     1       5.014  -4.590   0.000  1.00  0.00           C  
HETATM   40  C40 UNL     1       7.696  -1.129   0.000  1.00  0.00           C  
HETATM   41  C41 UNL     1       1.325  -3.121   0.000  1.00  0.00           C  
HETATM   42  C42 UNL     1       3.586   0.951   0.000  1.00  0.00           C  
HETATM   43  C43 UNL     1       3.464  -1.623   0.000  1.00  0.00           C  
HETATM   44  C44 UNL     1       5.616  -3.657   0.000  1.00  0.00           C  
HETATM   45  C45 UNL     1       7.595   1.062   0.000  1.00  0.00           C  
HETATM   46  C46 UNL     1       7.159   4.560   0.000  1.00  0.00           C  
HETATM   47  C47 UNL     1       2.652   4.522   0.000  1.00  0.00           C  
HETATM   48  C48 UNL     1      -0.034  -1.800   0.000  1.00  0.00           C  
HETATM   49  C49 UNL     1       5.083  -1.180   0.000  1.00  0.00           C  
HETATM   50  C50 UNL     1       9.997  -1.403   0.000  1.00  0.00           C  
HETATM   51  C51 UNL     1       3.500   2.632   0.000  1.00  0.00           C  
HETATM   52  C52 UNL     1       3.054  -0.296   0.000  1.00  0.00           C  
HETATM   53  C53 UNL     1       0.722   1.222   0.000  1.00  0.00           C  
HETATM   54  C54 UNL     1       9.916   0.543   0.000  1.00  0.00           C  
HETATM   55  C55 UNL     1       5.532   4.201   0.000  1.00  0.00           C  
HETATM   56  C56 UNL     1       8.919   1.993   0.000  1.00  0.00           C  
HETATM   57  C57 UNL     1       1.964   2.606   0.000  1.00  0.00           C  
HETATM   58  C58 UNL     1       4.581   5.211   0.000  1.00  0.00           C  
HETATM   59  C59 UNL     1       0.341   1.367   0.000  1.00  0.00           C  
HETATM   60  C60 UNL     1       5.542   3.680   0.000  1.00  0.00           C  
CONECT    1   13   13   16   36
CONECT    2   16   51   51   53
CONECT    3   12   12   25   50
CONECT    4   10   19   40   40
CONECT    5    9   29   29   47
CONECT    6   20   31   31   55
CONECT    7   17   25   25   28
CONECT    8   19   21   21   33
CONECT    9   53   59   59
CONECT   10   11   11   49
CONECT   11   14   22
CONECT   12   23   24
CONECT   13   32   44
CONECT   14   37   37   52
CONECT   15   37   38   46   46
CONECT   16   52   52
CONECT   17   24   24   43
CONECT   18   28   34   41   41
CONECT   19   26   26
CONECT   20   46   56   56
CONECT   21   25   27
CONECT   22   45   45   46
CONECT   23   31   54   54
CONECT   24   30
CONECT   26   39   44
CONECT   27   28   28   39
CONECT   29   38   51
CONECT   30   42   42   60
CONECT   31   60
CONECT   32   39   39   41
CONECT   33   40   50   50
CONECT   34   35   43   43
CONECT   35   48   48   59
CONECT   36   48   53   53
CONECT   37   51
CONECT   38   55   55
CONECT   40   45
CONECT   41   48
CONECT   42   43   57
CONECT   44   49   49
CONECT   45   56
CONECT   47   57   57   58
CONECT   49   52
CONECT   50   54
CONECT   54   56
CONECT   55   58
CONECT   57   59
CONECT   58   60   60
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
