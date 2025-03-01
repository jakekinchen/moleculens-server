// Scientific visualization: Cyclohexane Chair Conformations Explored
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
// Create a group for the scene title text and caption overlays
window.sceneTitleTextGroup = new THREE.Group();

// Load a font (using the built-in helvetiker font in Three.js)
const fontLoader = new THREE.FontLoader();
fontLoader.load('fonts/helvetiker_regular.typeface.json', function(font) {

  // Create the text geometry for the main title
  const titleGeometry = new THREE.TextGeometry("Cyclohexane Chair Conformations Explored", {
    font: font,
    size: 1,
    height: 0.2,
    curveSegments: 12,
    bevelEnabled: true,
    bevelThickness: 0.05,
    bevelSize: 0.05,
    bevelOffset: 0,
    bevelSegments: 5
  });
  // Center the title text so it is positioned at the origin
  titleGeometry.center();

  // Create a material with a soft white color and a glow-like emissive effect.
  // Start fully transparent for the fade-in effect.
  const titleMaterial = new THREE.MeshPhongMaterial({
    color: 0xf5f5f5,
    emissive: 0xffffff,
    emissiveIntensity: 0.5,
    transparent: true,
    opacity: 0.0
  });

  const titleMesh = new THREE.Mesh(titleGeometry, titleMaterial);
  // Position the title mesh at the center of the 3D space
  titleMesh.position.set(0, 0.5, 0);
  window.sceneTitleTextGroup.add(titleMesh);

  // Create a caption overlay geometry to guide the viewer
  const captionGeometry = new THREE.TextGeometry("Explore the molecular journey", {
    font: font,
    size: 0.3,
    height: 0.05,
    curveSegments: 10,
    bevelEnabled: true,
    bevelThickness: 0.02,
    bevelSize: 0.02,
    bevelOffset: 0,
    bevelSegments: 3
  });
  // Center the caption text and position it just below the main title
  captionGeometry.center();

  // Create a similar material for the caption with initial transparency
  const captionMaterial = new THREE.MeshPhongMaterial({
    color: 0xf5f5f5,
    emissive: 0xffffff,
    emissiveIntensity: 0.5,
    transparent: true,
    opacity: 0.0
  });

  const captionMesh = new THREE.Mesh(captionGeometry, captionMaterial);
  captionMesh.position.set(0, -0.8, 0);
  window.sceneTitleTextGroup.add(captionMesh);

  // Add a simple fade-in transition for both title and caption.
  // This uses a simple setInterval to update the opacity smoothly.
  const fadeDuration = 2000; // duration in milliseconds
  const steps = 60;
  const stepTime = fadeDuration / steps;
  let opacity = 0;
  const fadeInterval = setInterval(function(){
    opacity += 1/steps;
    if(opacity >= 1) {
      opacity = 1;
      clearInterval(fadeInterval);
    }
    titleMaterial.opacity = opacity;
    captionMaterial.opacity = opacity;
  }, stepTime);

  // Finally, add the group to the scene (positioned centered in 3D space)
  scene.add(window.sceneTitleTextGroup);
});

// GeometryAgent LLM-generated code
// Create a parent group for the cyclohexane molecule and store it globally
const cyclohexane = new THREE.Group();
window.cyclohexaneMolecule = cyclohexane;

// Materials (using semi-transparent, glossy finishes)
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0xCCCCCC, shininess: 100, opacity: 0.8, transparent: true });
const hydrogenAxialMaterial = new THREE.MeshPhongMaterial({ color: 0x0000ff, shininess: 100, opacity: 0.8, transparent: true });
const hydrogenEquatorialMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 100, opacity: 0.8, transparent: true });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, shininess: 80 });
const lineMaterial = new THREE.LineBasicMaterial({ color: 0x888888 });

// Basic geometries for atoms
const carbonGeometry = new THREE.SphereGeometry(0.3, 32, 32);
const hydrogenGeometry = new THREE.SphereGeometry(0.2, 32, 32);

// Utility: Create a bond (cylinder) between two points
function createBond(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);
  // Position bond at the midpoint
  const midPoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  bond.position.copy(midPoint);
  // Align the cylinder with the bond direction (default is along Y)
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize());
  return bond;
}

// Create the initial stylized ring outline (flat hexagon)
const ringOutlineGeometry = new THREE.BufferGeometry();
const outlineVertices = [];
const flatRadius = 2;
for (let i = 0; i < 6; i++) {
  const angle = i * Math.PI / 3;
  const x = flatRadius * Math.cos(angle);
  const z = flatRadius * Math.sin(angle);
  outlineVertices.push(x, 0, z);
}
ringOutlineGeometry.setAttribute('position', new THREE.Float32BufferAttribute(outlineVertices, 3));
const ringOutline = new THREE.LineLoop(ringOutlineGeometry, lineMaterial);

// Create the detailed chair conformation model group
const chairGroup = new THREE.Group();
const numCarbons = 6;
const carbonPositions = []; // will store the computed positions for carbon atoms

// For the chair conformation, position the carbons in a distorted hexagon:
for (let i = 0; i < numCarbons; i++) {
  const angle = i * Math.PI / 3;
  // Use the same x/z as the flat ring but offset Y to simulate a chair:
  const x = flatRadius * Math.cos(angle);
  const z = flatRadius * Math.sin(angle);
  // Alternate y positions: even-indexed carbons raised, odd-indexed lowered
  const y = (i % 2 === 0) ? 0.5 : -0.5;
  const pos = new THREE.Vector3(x, y, z);
  carbonPositions.push(pos);
  
  // Create the carbon atom
  const carbonMesh = new THREE.Mesh(carbonGeometry, carbonMaterial);
  carbonMesh.position.copy(pos);
  chairGroup.add(carbonMesh);
}

// Create bonds between adjacent carbon atoms (cyclic)
for (let i = 0; i < numCarbons; i++) {
  const nextIndex = (i + 1) % numCarbons;
  const bond = createBond(carbonPositions[i], carbonPositions[nextIndex]);
  chairGroup.add(bond);
}

// Add hydrogens to each carbon, color-coded for axial (blue) and equatorial (red)
for (let i = 0; i < numCarbons; i++) {
  const carbonPos = carbonPositions[i];
  
  // Axial hydrogen: along Y axis (pointing outward based on the carbon's position)
  const axialOffset = new THREE.Vector3(0, (carbonPos.y >= 0 ? 0.8 : -0.8), 0);
  const axialPos = carbonPos.clone().add(axialOffset);
  const axialHydrogen = new THREE.Mesh(hydrogenGeometry, hydrogenAxialMaterial);
  axialHydrogen.position.copy(axialPos);
  chairGroup.add(axialHydrogen);
  const bondAxial = createBond(carbonPos, axialPos);
  chairGroup.add(bondAxial);
  
  // Equatorial hydrogen: radially outward in the xz-plane
  const equatorialDir = new THREE.Vector3(carbonPos.x, 0, carbonPos.z).normalize();
  const equatorialOffset = equatorialDir.multiplyScalar(0.8);
  const equatorialPos = carbonPos.clone().add(equatorialOffset);
  const equatorialHydrogen = new THREE.Mesh(hydrogenGeometry, hydrogenEquatorialMaterial);
  equatorialHydrogen.position.copy(equatorialPos);
  chairGroup.add(equatorialHydrogen);
  const bondEquatorial = createBond(carbonPos, equatorialPos);
  chairGroup.add(bondEquatorial);
}

// Add both the initial ring outline and the detailed chair conformation to the parent cyclohexane group
cyclohexane.add(ringOutline);
cyclohexane.add(chairGroup);

// The cyclohexane object is now fully constructed with dynamic conformation elements,
// and it will appear in the scene as part of "Cyclohexane Chair Conformations Explored".
// External animation routines can access window.cyclohexaneMolecule, for example to:
// - Slowly rotate the entire object (e.g., cyclohexane.rotation.y += delta)
// - Morph between the flat ring outline and the chair conformation (using tweening or morph targets)
// - Trigger the ring-flip sequence with motion trails and annotation overlays

// Finally, add the cyclohexane molecule to the scene
scene.add(cyclohexane);
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
// Get elapsed time and deltaTime for smooth animations
const elapsedTime = clock.getElapsedTime();
const deltaTime = clock.getDelta();

// ----------------------------
// At 00:00 - Scene Title Emerges and Ring Outline Materializes
// ----------------------------
if (window.Scene_Title_Text) {
  // Fade in the scene title over the first 3 seconds
  const titleFade = Math.min(elapsedTime / 3, 1);
  if (window.Scene_Title_Text.material) {
    window.Scene_Title_Text.material.opacity = titleFade;
  }
  // (Optional) If the text has a “textMesh” property, update it here.
}

if (window.Cyclohexane_Molecule_Model) {
  // For the first 5 seconds, have the model (ring outline) scale from 0 to full size
  if (elapsedTime < 5) {
    const scaleVal = elapsedTime / 5;
    window.Cyclohexane_Molecule_Model.scale.set(scaleVal, scaleVal, scaleVal);
  } else {
    window.Cyclohexane_Molecule_Model.scale.set(1, 1, 1);
  }
}

// ----------------------------
// At 00:20 - Displaying Cyclohexane's 3D Structure: Slow rotation and camera orbit
// ----------------------------
if (elapsedTime > 20 && elapsedTime <= 40 && window.Cyclohexane_Molecule_Model) {
  // Slowly rotate the molecule on its Y-axis
  window.Cyclohexane_Molecule_Model.rotation.y += 0.2 * deltaTime;

  // Camera orbit: Assume an orbit around the molecule center.
  if (window.camera) {
    const orbitRadius = 10;
    // Calculate orbit angle based on time after 20 seconds
    const orbitAngle = (elapsedTime - 20) * 0.5;
    window.camera.position.x = Math.cos(orbitAngle) * orbitRadius;
    window.camera.position.z = Math.sin(orbitAngle) * orbitRadius;
    window.camera.lookAt(window.Cyclohexane_Molecule_Model.position);
  }
}

// ----------------------------
// At 00:40 - Highlighting Chair with Axial/Equatorial: Transition to chair conformation
// ----------------------------
if (elapsedTime > 40 && elapsedTime <= 60 && window.Cyclohexane_Molecule_Model) {
  // Compute progress of the chair transition (0 to 1 over 20 seconds)
  const transitionProgress = Math.min((elapsedTime - 40) / 20, 1);
  // You might be using a morph target or some custom uniform to distort the flat ring into a chair:
  window.Cyclohexane_Molecule_Model.userData.chairTransition = transitionProgress;

  // Optionally, adjust a slight tilt to simulate the three-dimensional distortion
  window.Cyclohexane_Molecule_Model.rotation.x = transitionProgress * 0.5;

  // Fade in color-coded labels (assumed to be children with names including 'axial' or 'equatorial')
  window.Cyclohexane_Molecule_Model.traverse(function(child) {
    if (child.name && (child.name.toLowerCase().includes('axial') || child.name.toLowerCase().includes('equatorial'))) {
      if (child.material) {
        child.material.opacity = transitionProgress; // Fade in with the same progress
      }
    }
  });
}

// ----------------------------
// At 01:00 - Axial vs Equatorial Orientations: Zoom in and continue rotation
// ----------------------------
if (elapsedTime > 60 && elapsedTime <= 80 && window.Cyclohexane_Molecule_Model) {
  // Continue a gentle rotation to emphasize different bonds
  window.Cyclohexane_Molecule_Model.rotation.y += 0.3 * deltaTime;

  // Zoom in the camera gradually to focus on the chair conformation
  if (window.camera) {
    // Define starting and target camera positions
    const startDistance = 8;
    const targetDistance = 4;
    const zoomProgress = Math.min((elapsedTime - 60) / 20, 1);
    const currentDistance = startDistance + (targetDistance - startDistance) * zoomProgress;

    // Assume the camera orbits on the XZ-plane around the model center; preserve angle from previous orbit if desired
    const orbitAngle = (elapsedTime - 20) * 0.5; // reuse previous angle calculation
    window.camera.position.x = Math.cos(orbitAngle) * currentDistance;
    window.camera.position.z = Math.sin(orbitAngle) * currentDistance;
    window.camera.lookAt(window.Cyclohexane_Molecule_Model.position);
  }

  // (Optional) Show on-screen annotations regarding steric interactions and energy profiles here.
}

// ----------------------------
// At 01:20 - Visualizing a Dynamic Ring-Flip: Animate the inversion (ring flip)
// ----------------------------
if (elapsedTime > 80 && elapsedTime <= 100 && window.Cyclohexane_Molecule_Model) {
  // Calculate flip progress; 0 at start (1:20) to 1 at end (1:40)
  const flipProgress = Math.min((elapsedTime - 80) / 20, 1);
  // Animate the flip: smoothly rotate the model around the X-axis to simulate inversion
  window.Cyclohexane_Molecule_Model.rotation.x = flipProgress * Math.PI;

  // Create a motion trail effect by modulating a custom property (to be used by a shader or separate trail mesh)
  window.Cyclohexane_Molecule_Model.userData.trailOpacity = Math.abs(Math.sin(flipProgress * Math.PI));
  // (Optional) If a dedicated motion trail object exists (e.g., window.RingFlipTrail), update its opacity or animate it here.
}

// ----------------------------
// At 01:40 - Cyclohexane Conformer Energy Balance: Side-by-side conformations and summary text
// ----------------------------
if (elapsedTime > 100 && window.Cyclohexane_Molecule_Model) {
  // Simulate side-by-side display by shifting the model position horizontally
  // Here, we shift the current model left and assume a mirrored display for the post-flip conformation.
  window.Cyclohexane_Molecule_Model.position.x = -1;

  // Optionally, if the model has separate sub-components for pre- and post-flip conformers,
  // adjust their positions so both appear side-by-side.
}

if (elapsedTime > 100 && window.Scene_Title_Text) {
  // Update the summary caption text (assumed to be updated via some external text update mechanism)
  // Fade in the summary text "Cyclohexane Conformer Energy Balance" over 3 seconds
  const summaryFade = Math.min((elapsedTime - 100) / 3, 1);
  if (window.Scene_Title_Text.material) {
    window.Scene_Title_Text.material.opacity = summaryFade;
  }
}

// ----------------------------
// Final Fade-Out: Begin slowly fading out the entire scene after the summary is shown
// ----------------------------
if (elapsedTime > 105) {
  // Compute a global fade factor over 5 seconds past 105 seconds (fade from 105 to 110 seconds)
  const fadeFactor = 1 - Math.min((elapsedTime - 105) / 5, 1);

  // Fade out the title text
  if (window.Scene_Title_Text && window.Scene_Title_Text.material) {
    window.Scene_Title_Text.material.opacity *= fadeFactor;
  }

  // Fade out the cyclohexane model
  if (window.Cyclohexane_Molecule_Model) {
    window.Cyclohexane_Molecule_Model.traverse(function(child) {
      if (child.material && typeof child.material.opacity === 'number') {
        child.material.opacity *= fadeFactor;
      }
    });
  }

  // Optionally, fade out the camera’s background or any other scene-wide elements here.
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
