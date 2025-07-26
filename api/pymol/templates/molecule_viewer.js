function createMoleculeVisualization(THREE, scene, options = {}) {
    console.log('createMoleculeVisualization');

    // Configuration options with defaults
    const config = {
        enableAnnotations: {{ enable_annotations }},  // Toggle atomic annotations
        scaleFactor: 0.25,       // Scale factor to control molecule size
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
        style.textContent = `{{ molecule_styles }}`;
        document.head.appendChild(style);
    }

    // Add molecule label if not already present
    if (!document.getElementById('molecule-label')) {
        const container = document.querySelector('#container');
        if (container) {
            const labelContainer = document.createElement('div');
            labelContainer.innerHTML = `<div id="molecule-label">{{ molecule_name }}</div>`;
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
    const pdbData = `{{ pdb_data }}`;

    // Create and configure the PDB loader
    let loader;
    if (typeof THREE.PDBLoader !== 'undefined') {
        loader = new THREE.PDBLoader();
    } else if (typeof window !== 'undefined' && window.PDBLoader) {
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

                    // Create CSS2DObject and attach it directly to the scene
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

createMoleculeVisualization(THREE, scene, {
  camera,
  controls
});
