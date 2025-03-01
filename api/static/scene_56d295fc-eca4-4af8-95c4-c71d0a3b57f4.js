// Scientific visualization: Exploring Cyclohexane Chair Conformations
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
// Materials for the cyclohexane model
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, shininess: 100 });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x666666, shininess: 100 });
const axialColor = 0xffaa00;      // bright color for axial substituents
const equatorialColor = 0x00aaff; // cool tone for equatorial substituents

// Geometry for a carbon atom (a sphere)
const carbonGeometry = new THREE.SphereGeometry(0.5, 32, 32);

// Create a group for the cyclohexane molecule model
const cyclohexane = new THREE.Group();
window.cyclohexane = cyclohexane;

// Arrays to store carbon meshes, their flat and chair positions,
// bonds and substituent arrow helpers for later animation updates.
const carbonMeshes = [];
const flatPositions = [];
const chairPositions = [];
const bonds = [];
const bondIndices = []; // keep track of which carbons each bond connects
const axialArrows = [];
const equatorialArrows = [];

// Utility function to create a bond (a cylinder) between two points.
// This function creates a cylinder oriented along the Y direction,
// then positions and rotates it to connect the two given points.
function createBond(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);
  
  // Position bond at the midpoint.
  const midPoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  bond.position.copy(midPoint);
  
  // Align the bond with the direction vector.
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
  
  // Save the original length so it can update smoothly later.
  bond.userData.originalLength = length;
  return bond;
}

// Number of carbons in cyclohexane
const nCarbons = 6;
const radius = 3; // used for the initial flat ring

// Calculate positions for the flat and chair conformations.
// Flat: carbons arranged in a perfect regular hexagon (all z = 0).
// Chair: use the same (x,y) but alternate z-positions; even-index carbons get z = +1, odd get z = -1.
for (let i = 0; i < nCarbons; i++) {
  const angle = (i / nCarbons) * Math.PI * 2;
  // Flat conformation
  const flatPos = new THREE.Vector3(
    radius * Math.cos(angle),
    radius * Math.sin(angle),
    0
  );
  flatPositions.push(flatPos);
  
  // Chair conformation: alternating z values for a basic chair shape.
  const chairPos = new THREE.Vector3(
    radius * Math.cos(angle),
    radius * Math.sin(angle),
    (i % 2 === 0 ? 1 : -1)
  );
  chairPositions.push(chairPos);
  
  // Create the carbon atom mesh at flat position initially.
  const carbon = new THREE.Mesh(carbonGeometry, carbonMaterial);
  carbon.position.copy(flatPos);
  carbonMeshes.push(carbon);
  cyclohexane.add(carbon);
  
  // Create substituent arrows for each carbon.
  // Axial arrow: points vertically (up if even, down if odd)
  const axialDir = new THREE.Vector3(0, 0, (i % 2 === 0 ? 1 : -1)).normalize();
  const axialArrow = new THREE.ArrowHelper(axialDir, carbon.position.clone(), 1, axialColor);
  axialArrows.push(axialArrow);
  cyclohexane.add(axialArrow);
  
  // Equatorial arrow: points perpendicular to the radial direction.
  // Compute the radial direction (from center to carbon, projected on XY plane).
  const radial = flatPos.clone().setZ(0).normalize();
  const equatorialDir = new THREE.Vector3(-radial.y, radial.x, 0).normalize();
  const equatorialArrow = new THREE.ArrowHelper(equatorialDir, carbon.position.clone(), 1, equatorialColor);
  equatorialArrows.push(equatorialArrow);
  cyclohexane.add(equatorialArrow);
}

// Create bonds between adjacent carbon atoms (closing the ring)
for (let i = 0; i < nCarbons; i++) {
  const j = (i + 1) % nCarbons;
  const bond = createBond(flatPositions[i], flatPositions[j]);
  bonds.push(bond);
  bondIndices.push({ i: i, j: j });
  cyclohexane.add(bond);
}

// Finally add the cyclohexane molecule to the scene
scene.add(cyclohexane);

// Animation update function to smoothly morph between flat and chair conformations,
// update bonds and substituent arrow orientations.
window.updateCyclohexane = function(time) {
  // time is assumed in seconds; create a morph progress that oscillates between 0 and 1
  // This parameter can be used to drive the ring-flip and morph transitions.
  const morphProgress = 0.5 + 0.5 * Math.sin(time);
  
  // Update each carbon position (linear interpolation between flat and chair positions)
  for (let i = 0; i < nCarbons; i++) {
    const newPos = new THREE.Vector3().lerpVectors(flatPositions[i], chairPositions[i], morphProgress);
    carbonMeshes[i].position.copy(newPos);
    
    // Update axial arrow: always vertical, direction determined by the chair position sign.
    const axialDir = new THREE.Vector3(0, 0, (i % 2 === 0 ? 1 : -1)).normalize();
    axialArrows[i].position.copy(newPos);
    axialArrows[i].setDirection(axialDir);
    
    // Update equatorial arrow: recalc the radial direction based on the current (x, y).
    const radial = newPos.clone().setZ(0);
    if (radial.lengthSq() > 0.0001) {
      radial.normalize();
    } else {
      radial.set(1, 0, 0);
    }
    const equatorialDir = new THREE.Vector3(-radial.y, radial.x, 0).normalize();
    equatorialArrows[i].position.copy(newPos);
    equatorialArrows[i].setDirection(equatorialDir);
  }
  
  // Update bonds: reposition and reorient bonds based on the updated carbon positions.
  for (let k = 0; k < bonds.length; k++) {
    const indices = bondIndices[k];
    const posA = carbonMeshes[indices.i].position;
    const posB = carbonMeshes[indices.j].position;
    
    // Compute new midpoint and direction
    const newDir = new THREE.Vector3().subVectors(posB, posA);
    const newLength = newDir.length();
    newDir.normalize();
    const midPoint = new THREE.Vector3().addVectors(posA, posB).multiplyScalar(0.5);
    
    // Update bond's position and quaternion.
    const bond = bonds[k];
    bond.position.copy(midPoint);
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), newDir);
    
    // Update scale to reflect new bond length.
    const originalLength = bond.userData.originalLength;
    bond.scale.set(1, newLength / originalLength, 1);
  }
}; 

// (Optional) To see continuous rotation of the entire cyclohexane group, you can add an animation,
// for example, by updating cyclohexane.rotation.y in your render loop.
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
const elapsedTime = clock.getElapsedTime();
const delta = clock.getDelta();

// ------------------------------
// At 00:00 - Initial rotating flat molecule and title overlay
if (window.Cyclohexane_Molecule_Model) {
  // Rotate the flat cyclohexane model slowly around the y-axis
  if (elapsedTime < 20) {
    window.Cyclohexane_Molecule_Model.rotation.y += 0.5 * delta;
  }
}
if (window.TitleOverlay) {
  // Show title overlay "Cyclohexane Chair Conformations" briefly (fade out after ~2 sec)
  if (elapsedTime < 2) {
    window.TitleOverlay.visible = true;
    window.TitleOverlay.material.opacity = 1;
  } else if (elapsedTime < 4) {
    // Fade out the title overlay
    window.TitleOverlay.material.opacity = Math.max(0, 1 - ((elapsedTime - 2) / 2));
  } else {
    window.TitleOverlay.visible = false;
  }
}
// (Optional) Caption text overlay "Introducing cyclohexane" can be managed similarly


// ------------------------------
// At 00:20 - Zoom in on cyclohexane ring and display atom labels
if (elapsedTime >= 20 && elapsedTime < 40) {
  // Smoothly interpolate camera position for zoom-in effect
  if (window.camera) {
    // Assume initial camera z-position is 10 and target is 5
    let t = Math.min((elapsedTime - 20) / 20, 1);
    window.camera.position.z = 10 - 5 * t;
    // Optionally adjust camera lookAt if needed (e.g., window.camera.lookAt(target))
  }
  // Show atom labels and bond markers by gradually increasing their opacity
  if (window.AtomLabels) {
    // Assume AtomLabels is an object with material.opacity property
    window.AtomLabels.visible = true;
    window.AtomLabels.material.opacity = Math.min((elapsedTime - 20) / 5, 1);
  }
  // (Optional) A caption "Cyclohexane ring structure" can be triggered here.
}


// ------------------------------
// At 00:40 - Transition to chair conformation view
if (elapsedTime >= 40 && elapsedTime < 60) {
  if (window.Cyclohexane_Molecule_Model) {
    let t = (elapsedTime - 40) / 20; // normalized progress 0 -> 1
    // If the model has morph targets for chair conformation, use them:
    if (window.Cyclohexane_Molecule_Model.morphTargetInfluences && 
        window.Cyclohexane_Molecule_Model.morphTargetInfluences.length > 0) {
      window.Cyclohexane_Molecule_Model.morphTargetInfluences[0] = t;
    }
    // Additionally, tilt the model to give a three-dimensional perspective
    window.Cyclohexane_Molecule_Model.rotation.x = t * (Math.PI / 4);
  }
  // (Optional) Caption "Chair conformation view" can be displayed here.
}


// ------------------------------
// At 01:00 - Highlight axial vs. equatorial substituent positions
if (elapsedTime >= 60 && elapsedTime < 80) {
  // Ensure directional arrows for axial and equatorial positions are visible
  if (window.AxialArrow) {
    window.AxialArrow.visible = true;
    // Create a subtle pulsing effect
    const pulse = 1 + 0.1 * Math.sin(elapsedTime * 5);
    window.AxialArrow.scale.set(pulse, pulse, pulse);
  }
  if (window.EquatorialArrow) {
    window.EquatorialArrow.visible = true;
    const pulse = 1 + 0.1 * Math.sin(elapsedTime * 5);
    window.EquatorialArrow.scale.set(pulse, pulse, pulse);
  }
  // Optionally, color-code specific atoms on the model if the material supports it.
  // (Optional) Display caption "Axial vs. Equatorial".
}


// ------------------------------
// At 01:20 - Dynamic ring-flip animation sequence
if (elapsedTime >= 80 && elapsedTime < 100) {
  if (window.Cyclohexane_Molecule_Model) {
    let t = (elapsedTime - 80) / 20; // normalized progress for the ring flip
    // Animate a flip: gradually rotate the model's x-axis from current tilt to a flipped conformation.
    // For example, transition from a ~45° tilt to a full 180° flip (adjust values as needed)
    window.Cyclohexane_Molecule_Model.rotation.x = (Math.PI / 4) + t * (Math.PI - (Math.PI / 4));
  }
  // Reverse or rotate the axial/equatorial arrows to indicate the interchange:
  if (window.AxialArrow) {
    window.AxialArrow.rotation.z = t * Math.PI;
  }
  if (window.EquatorialArrow) {
    window.EquatorialArrow.rotation.z = t * Math.PI;
  }
  // (Optional) Caption "Dynamic ring-flip" can remain displayed during this sequence.
}


// ------------------------------
// At 01:40 - Final scene: Display both chair conformations and summary overlay
if (elapsedTime >= 100) {
  // Clone the cyclohexane model to display the alternative chair conformation side by side.
  if (window.Cyclohexane_Molecule_Model && !window.SecondChair) {
    window.SecondChair = window.Cyclohexane_Molecule_Model.clone();
    // Offset the clone's position to the right to create a side-by-side view
    window.SecondChair.position.x = window.Cyclohexane_Molecule_Model.position.x + 3;
    // Add to the scene if available
    if (window.scene) {
      window.scene.add(window.SecondChair);
    }
  }
  // Display minimal energy diagrams and comparative labels.
  if (window.EnergyDiagram) {
    window.EnergyDiagram.visible = true;
    // Fade in the energy diagram over 2 seconds starting at 100 sec elapsed time.
    window.EnergyDiagram.material.opacity = Math.min((elapsedTime - 100) / 2, 1);
  }
  // After 110 sec, fade out the entire models (both chair conformations) for a concluding recapitulation.
  if (elapsedTime > 110) {
    const fade = Math.max(0, 1 - (elapsedTime - 110) / 5);
    // Fade out the main model
    window.Cyclohexane_Molecule_Model.traverse((child) => {
      if (child.material && child.material.opacity !== undefined) {
        child.material.opacity = fade;
      }
    });
    // Fade out the second chair conformation clone
    if (window.SecondChair) {
      window.SecondChair.traverse((child) => {
        if (child.material && child.material.opacity !== undefined) {
          child.material.opacity = fade;
        }
      });
    }
  }
  // (Optional) Display caption "Conformational summary" during this final phase.
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
