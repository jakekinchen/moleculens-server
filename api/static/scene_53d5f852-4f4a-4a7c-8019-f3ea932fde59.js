// Scientific visualization: Understanding Cyclohexane Chair Conformations
// Auto-generated code

// Global variables
let scene, camera, renderer, controls, clock;
let isPlaying = true;
const TOTAL_DURATION = 120; // 2 minutes in seconds


// Initialize the scene
function init() {
    // Create scene
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x111122);
    
    // Create camera
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    camera.position.z = 10;
    
    // Create renderer
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setPixelRatio(window.devicePixelRatio);
    document.body.appendChild(renderer.domElement);
    
    // Create lighting
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);
    
    const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
    directionalLight.position.set(1, 1, 1).normalize();
    scene.add(directionalLight);
    
    // Add orbit controls
    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.25;
    
    // Initialize clock
    clock = new THREE.Clock();
    
    // Create geometry
    createGeometry();
    
    // Handle window resize
    window.addEventListener('resize', onWindowResize);
    
    // Set up UI controls
    setupControls();
}

// Create all geometry in the scene
function createGeometry() {
    // Geometry created by the GeometryAgent
// GeometryAgent LLM-generated code
(function() {
  // Materials
  const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x777777, flatShading: true });
  const axialMaterial = new THREE.MeshPhongMaterial({ color: 0x0000ff });
  const equatorialMaterial = new THREE.MeshPhongMaterial({ color: 0x00ff00 });
  const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, flatShading: true });

  // Geometries for atoms and markers
  const carbonGeometry = new THREE.SphereGeometry(0.5, 32, 32);
  const markerGeometry = new THREE.SphereGeometry(0.2, 16, 16);

  // Create the cyclohexane group (molecule)
  const cyclohexane = new THREE.Group();
  window.cyclohexane = cyclohexane;
  scene.add(cyclohexane);

  // Arrays to store carbon atoms and bond meshes
  const carbons = [];
  const bonds = [];
  
  // Positions for flat and chair conformations, for 6 carbon atoms.
  // Flat conformation: carbons arranged in a perfect hexagon (radius = 5 units)
  // Chair conformation: slightly distorted ring: alternate carbons lie at different radial distances and with an offset in z
  const flatPositions = [];
  const chairPositions = [];
  const numCarbons = 6;
  const twoPi = Math.PI * 2;
  for (let i = 0; i < numCarbons; i++) {
    const angle = (i / numCarbons) * twoPi;
    // Flat (planar) positions: on circle of radius 5 in x-y plane, z = 0.
    const flatX = 5 * Math.cos(angle);
    const flatY = 5 * Math.sin(angle);
    flatPositions.push(new THREE.Vector3(flatX, flatY, 0));
    // Chair conformation: alternate carbons with slightly different radii and z offsets.
    // For even index: radius 4.5 and z = +1; for odd: radius 5.5 and z = -1.
    const r = (i % 2 === 0) ? 4.5 : 5.5;
    const chairX = r * Math.cos(angle);
    const chairY = r * Math.sin(angle);
    const chairZ = (i % 2 === 0) ? 1 : -1;
    chairPositions.push(new THREE.Vector3(chairX, chairY, chairZ));
  }

  // Create carbon atoms with labeling markers
  for (let i = 0; i < numCarbons; i++) {
    // Create a sphere for a carbon atom. Start at flat conformation.
    const carbon = new THREE.Mesh(carbonGeometry, carbonMaterial);
    carbon.position.copy(flatPositions[i]);
    // Store both conformations so we can morph later.
    carbon.userData = {
      flatPos: flatPositions[i].clone(),
      chairPos: chairPositions[i].clone()
    };

    // Create axial marker.
    // In chair conformation each carbon carries one axial label,
    // offset vertically; even-index carbons have axial marker up, odd-index down.
    const axialMarker = new THREE.Mesh(markerGeometry, axialMaterial);
    const axialOffset = (i % 2 === 0) ? 0.8 : -0.8;
    axialMarker.position.set(0, 0, axialOffset);
    // Create equatorial marker.
    // Place it in the x-y plane, offset outward from the carbon center.
    const angle = (i / numCarbons) * twoPi;
    const ex = Math.cos(angle) * 1.0;
    const ey = Math.sin(angle) * 1.0;
    const equatorialMarker = new THREE.Mesh(markerGeometry, equatorialMaterial);
    equatorialMarker.position.set(ex, ey, 0);

    // Add markers as children of the carbon atom for easy transformation.
    carbon.add(axialMarker);
    carbon.add(equatorialMarker);

    // Add the carbon to the cyclohexane group and store in array.
    cyclohexane.add(carbon);
    carbons.push(carbon);
  }

  // Function to create a bond (cylinder) between two points.
  function createBond(start, end) {
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();
    const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
    const bond = new THREE.Mesh(bondGeometry, bondMaterial);
    // Position: midpoint between start and end.
    bond.position.copy(start).lerp(end, 0.5);
    // Orient the cylinder so that it connects start to end.
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
    return bond;
  }

  // Create bonds connecting neighboring carbons in the ring (and closing the ring)
  for (let i = 0; i < numCarbons; i++) {
    const nextIndex = (i + 1) % numCarbons;
    // Create a bond based on current positions (flat conformation initially)
    const bond = createBond(carbons[i].position, carbons[nextIndex].position);
    // Store endpoints indices so that we can update bonds during morph animations.
    bond.userData = { a: i, b: nextIndex };
    cyclohexane.add(bond);
    bonds.push(bond);
  }

  // Function to update bond geometry based on new positions.
  function updateBond(bond, start, end) {
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();
    bond.position.copy(start).lerp(end, 0.5);
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
    bond.geometry.dispose();
    bond.geometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
  }

  // Animation update function:
  // This function should be called in your render loop with the elapsed time.
  // It smoothly morphs the cyclohexane molecule between its flat and chair conformations
  // and continuously rotates the molecule.
  window.updateCyclohexaneMolecule = function(time) {
    // t oscillates between 0 and 1 over time (morph factor)
    const t = 0.5 * (1 + Math.sin(time));
    // Update each carbon atom by interpolating between flat and chair position.
    for (let i = 0; i < numCarbons; i++) {
      const carbon = carbons[i];
      const flatPos = carbon.userData.flatPos;
      const chairPos = carbon.userData.chairPos;
      // Linear interpolation: pos = flat*(1-t) + chair*(t)
      carbon.position.set(
        flatPos.x * (1 - t) + chairPos.x * t,
        flatPos.y * (1 - t) + chairPos.y * t,
        flatPos.z * (1 - t) + chairPos.z * t
      );
      // (Markers are children so they move with the carbon)
    }
    // Update bonds based on current positions of their endpoints.
    for (let i = 0; i < bonds.length; i++) {
      const bond = bonds[i];
      const a = bond.userData.a;
      const b = bond.userData.b;
      updateBond(bond, carbons[a].position, carbons[b].position);
    }
    // Apply a continuous rotation to the entire molecule about its vertical (z) axis.
    cyclohexane.rotation.z += 0.005;
  };

  // Optionally, if you have an animation loop elsewhere, call:
  // function animate(time) {
  //   requestAnimationFrame(animate);
  //   window.updateCyclohexaneMolecule(time * 0.001);
  //   renderer.render(scene, camera);
  // }
  // animate();
})();

// GeometryAgent LLM-generated code
// Create a group for the energy diagram overlay
window.energyDiagramOverlay = new THREE.Group();
const overlay = window.energyDiagramOverlay;

// (Optional) Set custom userData for animation parameters (fade-in, etc.)
overlay.userData = {
  animation: {
    type: "fade-in",
    duration: 2000,  // in ms
    description: "Gradual fade-in overlay with dynamic energy level indicators synchronized with the molecule's flip"
  }
};

// Create a semi-transparent, glossy background plane (8 units wide x 6 units high)
const planeGeometry = new THREE.PlaneGeometry(8, 6);
const planeMaterial = new THREE.MeshPhongMaterial({
  color: 0x222222,
  transparent: true,
  opacity: 0.5,
  shininess: 100,
});
const backgroundPlane = new THREE.Mesh(planeGeometry, planeMaterial);
// Slightly offset the plane so that overlay elements appear on top
backgroundPlane.position.set(0, 0, 0);
overlay.add(backgroundPlane);

// Create energy level indicators as thin horizontal bars
const barWidth = 6;
const barHeight = 0.1;
const barDepth = 0.1;
const barGeometry = new THREE.BoxGeometry(barWidth, barHeight, barDepth);
const barMaterial = new THREE.MeshPhongMaterial({ color: 0xffcc00, shininess: 80 });

// High energy level bar (upper level)
const highEnergyBar = new THREE.Mesh(barGeometry, barMaterial);
highEnergyBar.position.set(0, 1, 0.1); // positioned on the background plane
overlay.add(highEnergyBar);

// Low energy level bar (lower level) – highlighting the favorable chair conformation
const lowEnergyBar = new THREE.Mesh(barGeometry, barMaterial);
lowEnergyBar.position.set(0, -1, 0.1);
overlay.add(lowEnergyBar);

// Create a bright arrow to indicate the energy difference (from high to low energy)
const arrowDir = new THREE.Vector3(0, -1, 0);
const arrowOrigin = new THREE.Vector3(0, 1, 0.2); // starting at the high energy level (slightly in front)
const arrowLength = 2;  // distance between energy levels
const arrowColor = 0xffffff;
const headLength = 0.4;
const headWidth = 0.2;
const energyArrow = new THREE.ArrowHelper(arrowDir, arrowOrigin, arrowLength, arrowColor, headLength, headWidth);
overlay.add(energyArrow);

// Create an annotation for the lower energy state using a canvas texture
function createTextSprite(message, parameters) {
  parameters = parameters || {};
  const fontface = parameters.fontface || "Arial";
  const fontsize = parameters.fontsize || 24;
  const borderThickness = parameters.borderThickness || 4;
  const borderColor = parameters.borderColor || { r: 0, g: 0, b: 0, a: 1.0 };
  const backgroundColor = parameters.backgroundColor || { r: 255, g: 255, b: 255, a: 1.0 };

  const canvas = document.createElement("canvas");
  const context = canvas.getContext("2d");
  context.font = "Bold " + fontsize + "px " + fontface;
  // get size data (height depends only on font size)
  const metrics = context.measureText(message);
  const textWidth = metrics.width;
  canvas.width = textWidth + borderThickness * 2;
  canvas.height = fontsize * 1.4 + borderThickness * 2;

  // need to reset font since we changed canvas size
  context.font = "Bold " + fontsize + "px " + fontface;
  // background color
  context.fillStyle = "rgba(" + backgroundColor.r + "," + backgroundColor.g + "," + backgroundColor.b + "," + backgroundColor.a + ")";
  context.fillRect(0, 0, canvas.width, canvas.height);
  // border
  context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + "," + borderColor.b + "," + borderColor.a + ")";
  context.lineWidth = borderThickness;
  context.strokeRect(0, 0, canvas.width, canvas.height);
  // text color
  context.fillStyle = "rgba(0, 0, 0, 1.0)";
  context.textAlign = "center";
  context.textBaseline = "middle";
  context.fillText(message, canvas.width / 2, canvas.height / 2);

  const texture = new THREE.CanvasTexture(canvas);
  texture.needsUpdate = true;

  const spriteMaterial = new THREE.SpriteMaterial({ map: texture, transparent: true });
  const sprite = new THREE.Sprite(spriteMaterial);
  // scale sprite based on canvas dimensions
  const scaleFactor = 0.01;  // adjust as needed for scene scale
  sprite.scale.set(canvas.width * scaleFactor, canvas.height * scaleFactor, 1);

  return sprite;
}

const annotation = createTextSprite("Lower Energy (Chair)", {
  fontsize: 32,
  borderThickness: 2,
  backgroundColor: { r: 255, g: 255, b: 255, a: 0.8 },
  borderColor: { r: 0, g: 0, b: 0, a: 1.0 },
});

// Position the annotation near the low energy bar, offset to the left and above it
annotation.position.set(-2.5, -0.8, 0.3);
overlay.add(annotation);

// (Optional) Store overlay elements for animation purposes
overlay.userData.elements = {
  backgroundPlane: backgroundPlane,
  highEnergyBar: highEnergyBar,
  lowEnergyBar: lowEnergyBar,
  energyArrow: energyArrow,
  annotation: annotation,
};

// Finally, add the energy diagram overlay group to the scene
scene.add(overlay);
}

// Animation loop
function animate() {
    if (!isPlaying) {
        requestAnimationFrame(animate);
        renderer.render(scene, camera);
        controls.update();
        return;
    }
    
    requestAnimationFrame(animate);
    
    // Animation created by the AnimationAgent
// Get elapsed time and delta time
const elapsedTime = clock.getElapsedTime();
const deltaTime = clock.getDelta();

// ────────────────────────────────────────────────────────────
// At 00:00 - Introduction: Show cyclohexane 6-membered ring (0s - 20s)
if (elapsedTime < 20) {
  // Find the cyclohexane molecule
  const molecule = scene.getObjectByName("Cyclohexane_Molecule");
  if (molecule) {
    // Slowly rotate the molecule to emphasize its cyclic structure
    molecule.rotation.y += 0.5 * deltaTime;
    
    // (Optional) If the molecule uses material opacity, gradually fade it in
    molecule.traverse(child => {
      if (child.isMesh && child.material && child.material.opacity !== undefined) {
        child.material.opacity = Math.min(1, elapsedTime / 5);
        child.material.transparent = child.material.opacity < 1;
      }
    });
  }
  // Position the camera to get an overall view of the molecule
  if (camera) {
    // Set a gentle circular orbit around the origin
    const radius = 8;
    const angle = elapsedTime * 0.2;
    camera.position.x = Math.cos(angle) * radius;
    camera.position.z = Math.sin(angle) * radius;
    camera.position.y = 3; // slightly above for better view
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "Cyclohexane: 6-membered ring")
}

// ────────────────────────────────────────────────────────────
// At 00:20 - Morphing from flat representation to chair conformation (20s - 40s)
else if (elapsedTime >= 20 && elapsedTime < 40) {
  const progress = (elapsedTime - 20) / 20; // Progress from 0 to 1
  const molecule = scene.getObjectByName("Cyclohexane_Molecule");
  if (molecule) {
    // Assume the molecule has a property 'morphProgress' that drives its conformation state
    molecule.userData.morphProgress = progress;
    // For visual effect, a slight rotation adjustment during morphing
    molecule.rotation.y += 0.3 * deltaTime;
  }
  // Keep the camera in a similar orbit for continuity
  if (camera) {
    const radius = 8;
    const angle = elapsedTime * 0.2;
    camera.position.x = Math.cos(angle) * radius;
    camera.position.z = Math.sin(angle) * radius;
    camera.position.y = 3;
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "Chair shape: optimal structure")
}

// ────────────────────────────────────────────────────────────
// At 00:40 - Focus on carbons with labels (40s - 60s)
else if (elapsedTime >= 40 && elapsedTime < 60) {
  const progress = (elapsedTime - 40) / 20; // Progress 0 to 1 for label appearance
  const molecule = scene.getObjectByName("Cyclohexane_Molecule");
  if (molecule) {
    // Continue a subtle rotation for smooth viewing
    molecule.rotation.y += 0.4 * deltaTime;
    
    // Traverse molecule to find individual carbon labels
    molecule.traverse(child => {
      // Assume that each carbon’s label mesh has a name including "Label"
      if (child.isMesh && child.name.includes("Label") && child.material && child.material.opacity !== undefined) {
        // Gradually reveal the label with axial/equatorial info
        child.material.opacity = Math.min(1, progress);
        child.material.transparent = child.material.opacity < 1;
        // Optionally, scale the label in for a subtle pop-in effect
        const scale = 0.5 + 0.5 * progress;
        child.scale.set(scale, scale, scale);
      }
    });
  }
  // Adjust camera closer to the molecule to inspect details
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 2, 4);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "Axial vs Equatorial")
}

// ────────────────────────────────────────────────────────────
// At 01:00 - Chair flip dynamics (60s - 90s)
else if (elapsedTime >= 60 && elapsedTime < 90) {
  const progress = (elapsedTime - 60) / 30; // 0 to 1 progress during chair flip
  const molecule = scene.getObjectByName("Cyclohexane_Molecule");
  if (molecule) {
    // Perform a smooth flip along the X-axis (simulate chair flip)
    // This rotates from 0 to Math.PI (180°) gradually
    molecule.rotation.x = progress * Math.PI;
    
    // Optionally highlight key atoms or bonds (assume they have a property to modify emissive color)
    molecule.traverse(child => {
      if (child.isMesh && child.material) {
        // Briefly pulse emissive intensity during the flip
        child.material.emissive = new THREE.Color(0x333333);
        child.material.emissiveIntensity = 0.5 + 0.5 * Math.sin(progress * Math.PI);
      }
    });
  }
  // Slowly adjust the camera to follow the flip dynamics
  if (camera) {
    const angle = progress * Math.PI * 2;
    const radius = 6;
    camera.position.x = Math.cos(angle) * radius;
    camera.position.z = Math.sin(angle) * radius;
    camera.position.y = 3 + progress; // slight elevation change
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "Chair flip dynamics")
}

// ────────────────────────────────────────────────────────────
// At 01:30 - Energy diagram overlay appears (90s - 110s)
else if (elapsedTime >= 90 && elapsedTime < 110) {
  const progress = (elapsedTime - 90) / 20; // 0 to 1 progress for the overlay's appearance
  // Continue rotating the molecule
  const molecule = scene.getObjectByName("Cyclohexane_Molecule");
  if (molecule) {
    molecule.rotation.y += 0.3 * deltaTime;
  }
  // Animate the Energy Diagram Overlay
  const energyOverlay = scene.getObjectByName("Energy_Diagram_Overlay");
  if (energyOverlay) {
    // Fade in the overlay (assume it uses material.opacity)
    energyOverlay.traverse(child => {
      if (child.isMesh && child.material && child.material.opacity !== undefined) {
        child.material.opacity = Math.min(1, progress);
        child.material.transparent = child.material.opacity < 1;
      }
    });
    // Move the overlay into position (e.g., right side of the molecule)
    energyOverlay.position.lerp(new THREE.Vector3(4, 2, 0), 0.05);
  }
  // Place the camera to see both the molecule and the overlay
  if (camera) {
    const targetPosition = new THREE.Vector3(2, 3, 8);
    camera.position.lerp(targetPosition, 0.02);
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "Energy & strain explained")
}

// ────────────────────────────────────────────────────────────
// At 01:50 - Split-screen display showing both chair conformations (110s - 120s)
else if (elapsedTime >= 110 && elapsedTime < 120) {
  const progress = (elapsedTime - 110) / 10; // fast transition in the final 10 seconds
  const molecule = scene.getObjectByName("Cyclohexane_Molecule");
  if (molecule) {
    // Create (or retrieve) a clone of the molecule for the alternate conformation
    let moleculeClone = scene.getObjectByName("Cyclohexane_Clone");
    if (!moleculeClone) {
      moleculeClone = molecule.clone();
      moleculeClone.name = "Cyclohexane_Clone";
      // Offset the clone slightly initially; its final position will be set below.
      moleculeClone.position.x = 2;
      scene.add(moleculeClone);
    }
    // Animate repositioning to display side-by-side conformations
    // Main molecule moves to the left; clone moves to the right.
    molecule.position.x = THREE.MathUtils.lerp(0, -2, progress);
    moleculeClone.position.x = THREE.MathUtils.lerp(2, 2, progress); // fixed on the right side
    
    // Optionally, add a subtle rotation to each for better viewing
    molecule.rotation.y += 0.2 * deltaTime;
    moleculeClone.rotation.y += 0.2 * deltaTime;
  }
  // Slowly revolve the camera around the split-screen display for summary view
  if (camera) {
    const angle = elapsedTime * 0.1;
    const radius = 7;
    camera.position.x = Math.cos(angle) * radius;
    camera.position.z = Math.sin(angle) * radius;
    camera.position.y = 3;
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "Conformation impacts reactivity")
}
    
    // Update controls
    controls.update();
    
    // Update captions and timeline
    updateUI();
    
    // Render the scene
    renderer.render(scene, camera);
}

// Update captions and UI based on current time
function updateUI() {
    const elapsedTime = clock.getElapsedTime();
    
    // Update progress bar
    const progressBar = document.getElementById('progress-bar');
    const progress = Math.min(elapsedTime / TOTAL_DURATION, 1);
    progressBar.style.width = progress * 100 + '%';
    
    // Update captions
    const captions = document.querySelectorAll('.caption');
    captions.forEach(caption => {
        const timeStr = caption.getAttribute('data-time');
        const [min, sec] = timeStr.split(':').map(Number);
        const timeInSeconds = min * 60 + sec;
        
        // Show caption if we're within 5 seconds of its timecode
        if (elapsedTime >= timeInSeconds && elapsedTime < timeInSeconds + 5) {
            caption.style.display = 'block';
        } else {
            caption.style.display = 'none';
        }
    });
}

// Window resize handler
function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

// Set up UI controls
function setupControls() {
    const playPauseButton = document.getElementById('play-pause');
    const resetButton = document.getElementById('reset');
    
    playPauseButton.addEventListener('click', () => {
        isPlaying = !isPlaying;
        playPauseButton.textContent = isPlaying ? 'Pause' : 'Play';
        if (isPlaying) {
            clock.start();
        } else {
            clock.stop();
        }
    });
    
    resetButton.addEventListener('click', () => {
        clock = new THREE.Clock();
        isPlaying = true;
        playPauseButton.textContent = 'Pause';
    });
}

// Initialize and start animation
init();
animate();
