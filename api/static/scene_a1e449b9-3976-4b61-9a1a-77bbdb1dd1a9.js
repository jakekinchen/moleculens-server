// Scientific visualization: Decoding Ethene's Double Bond: A 3D Exploration
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
// Materials for atoms, bonds, and electron density overlay
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x555555, shininess: 100 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 100 });
const singleBondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, shininess: 100 });
const doubleBondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, emissive: 0xffff00, shininess: 100 });
const electronDensityMaterial = new THREE.MeshPhongMaterial({
  color: 0x0000ff,
  transparent: true,
  opacity: 0.3,
  side: THREE.DoubleSide
});

// Base sphere geometry (unit sphere) for atoms
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Helper function to create bonds (cylinders)
// Optionally accepts an offset vector to displace the bond (useful for double bond components).
function createBond(start, end, material, radius = 0.08, offset) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  const bondGeometry = new THREE.CylinderGeometry(radius, radius, length, 16);
  const bond = new THREE.Mesh(bondGeometry, material);

  // Position the bond at the midpoint between start and end.
  bond.position.copy(start).addScaledVector(direction, 0.5);

  // Align the cylinder: default cylinder is along the Y-axis.
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());

  // Apply additional offset if provided.
  if (offset) {
    bond.position.add(offset);
  }
  return bond;
}

// Create a group for the ethene molecule
const etheneMolecule = new THREE.Group();
window.etheneMolecule = etheneMolecule; // store globally

// Define positions for atoms (values chosen for molecular visualization)
// Carbon atoms: placed along the X-axis.
const carbon1Pos = new THREE.Vector3(-0.75, 0, 0);
const carbon2Pos = new THREE.Vector3(0.75, 0, 0);

// Hydrogen atoms: ethene (CH2=CH2) has two hydrogens per carbon.
// They are arranged approximately 120° apart in the plane (xy-plane).
const h1Pos = new THREE.Vector3(-1.5, 0.9, 0);
const h2Pos = new THREE.Vector3(-1.5, -0.9, 0);
const h3Pos = new THREE.Vector3(1.5, 0.9, 0);
const h4Pos = new THREE.Vector3(1.5, -0.9, 0);

// Create carbon atoms with a polished spherical appearance.
const carbon1 = new THREE.Mesh(atomGeometry, carbonMaterial);
carbon1.position.copy(carbon1Pos);
carbon1.scale.set(0.5, 0.5, 0.5); // scale for carbon

const carbon2 = new THREE.Mesh(atomGeometry, carbonMaterial);
carbon2.position.copy(carbon2Pos);
carbon2.scale.set(0.5, 0.5, 0.5);

// Create hydrogen atoms.
const hydrogen1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen1.position.copy(h1Pos);
hydrogen1.scale.set(0.3, 0.3, 0.3);

const hydrogen2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen2.position.copy(h2Pos);
hydrogen2.scale.set(0.3, 0.3, 0.3);

const hydrogen3 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen3.position.copy(h3Pos);
hydrogen3.scale.set(0.3, 0.3, 0.3);

const hydrogen4 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen4.position.copy(h4Pos);
hydrogen4.scale.set(0.3, 0.3, 0.3);

// Create bonds between atoms

// Carbon to Hydrogen bonds (single bonds)
const bondC1H1 = createBond(carbon1Pos, h1Pos, singleBondMaterial);
const bondC1H2 = createBond(carbon1Pos, h2Pos, singleBondMaterial);
const bondC2H3 = createBond(carbon2Pos, h3Pos, singleBondMaterial);
const bondC2H4 = createBond(carbon2Pos, h4Pos, singleBondMaterial);

// Carbon-Carbon double bond:
// For the double bond, create two parallel cylinders.
// The first (sigma) bond: centered along the connection.
const ccBondSigma = createBond(carbon1Pos, carbon2Pos, doubleBondMaterial, 0.1);
// The second (pi) bond: offset slightly out-of-plane (along Z) to simulate the electron cloud above and below.
const ccBondPi = createBond(carbon1Pos, carbon2Pos, doubleBondMaterial, 0.1, new THREE.Vector3(0, 0, 0.12));

// Group all atoms and bonds together.
etheneMolecule.add(carbon1, carbon2, hydrogen1, hydrogen2, hydrogen3, hydrogen4,
  bondC1H1, bondC1H2, bondC2H3, bondC2H4,
  ccBondSigma, ccBondPi);

// Create an electron density overlay covering the molecule.
// This overlay (a semi-transparent circle with contour lines) visualizes the electron cloud.
// The circle is centered at the midpoint of the carbon-carbon bond.
const overlayCenter = new THREE.Vector3().addVectors(carbon1Pos, carbon2Pos).multiplyScalar(0.5);
const electronOverlayGeometry = new THREE.CircleGeometry(2, 64);
const electronOverlay = new THREE.Mesh(electronOverlayGeometry, electronDensityMaterial);
electronOverlay.position.copy(overlayCenter);
electronOverlay.rotation.x = -Math.PI / 2;  // lie in the xy-plane
// Slightly lower it along z to avoid z-fighting.
electronOverlay.position.z = -0.05;

// Create contour lines for the electron density overlay using LineLoop.
const contoursMaterial = new THREE.LineBasicMaterial({ color: 0x0000ff, transparent: true, opacity: 0.6 });
const contoursGeometry = new THREE.BufferGeometry().setFromPoints(electronOverlayGeometry.vertices || []);
const contourLines = new THREE.LineLoop(contoursGeometry, contoursMaterial);
contourLines.rotation.x = -Math.PI / 2;
contourLines.position.copy(overlayCenter);
contourLines.position.z = -0.04;

// Add the electron density overlay and contour lines to the molecule group.
etheneMolecule.add(electronOverlay, contourLines);

// Attach a simple continuous rotation animation to the molecule group.
// This animation rotates the molecule slowly around its Y-axis,
// revealing its planar geometry and emphasizing bonding details.
etheneMolecule.userData.animate = function(deltaTime) {
  this.rotation.y += deltaTime * 0.2;
};

// Finally, add the ethene molecule group to the global scene.
scene.add(etheneMolecule);
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
// Get elapsed time and deltaTime in seconds
const elapsedTime = clock.getElapsedTime();
const delta = clock.getDelta();

// ----- TIMECODE 00:00 - Introducing Ethene Molecule ----- 
if (window.Ethene_Molecule_Animation) {
  // Rotate molecule gradually to reveal planar geometry
  if (elapsedTime < 20) {
    window.Ethene_Molecule_Animation.rotation.y += 0.5 * delta;
  }
}

// ----- TIMECODE 00:20 - Highlighting the Double Bond -----
if (elapsedTime >= 20 && elapsedTime < 40) {
  // Camera zooms in on the carbon-carbon bond
  if (window.camera) {
    // Interpolate camera position.z from 10 (default) to 5 (close-up)
    window.camera.position.z = THREE.MathUtils.lerp(10, 5, (elapsedTime - 20) / 20);
  }
  // Highlight the double bond glow effect
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.doubleBond) {
    // Gradually increase emissive intensity for highlight effect
    window.Ethene_Molecule_Animation.doubleBond.material.emissiveIntensity = THREE.MathUtils.lerp(0, 2, (elapsedTime - 20) / 20);
  }
  // Optionally trigger split view elements (if defined)
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.splitView) {
    window.Ethene_Molecule_Animation.splitView.visible = true;
    window.Ethene_Molecule_Animation.splitView.material.opacity = THREE.MathUtils.lerp(0, 1, (elapsedTime - 20) / 20);
  }
}

// ----- TIMECODE 00:40 - Sigma & Pi Overlap -----
if (elapsedTime >= 40 && elapsedTime < 60) {
  let progress = Math.min((elapsedTime - 40) / 20, 1);
  
  // Animate sigma bond representation (head-on sp2 overlap)
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.sigmaBond) {
    window.Ethene_Molecule_Animation.sigmaBond.material.opacity = progress;
    // Optionally add a slight scaling effect to emphasize formation
    window.Ethene_Molecule_Animation.sigmaBond.scale.setScalar(THREE.MathUtils.lerp(0.8, 1, progress));
  }
  
  // Animate pi bond representation (lobes above and below)
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.piBond) {
    window.Ethene_Molecule_Animation.piBond.material.opacity = progress;
    // The electron clouds may slide into position; simulate with slight vertical translation
    window.Ethene_Molecule_Animation.piBond.position.y = THREE.MathUtils.lerp(0, 0.5, progress);
  }
  
  // Optionally animate arrows/lines indicating orbital overlaps if available
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.orbitalArrows) {
    window.Ethene_Molecule_Animation.orbitalArrows.material.opacity = progress;
  }
}

// ----- TIMECODE 01:00 - Mapping Electron Density -----
if (elapsedTime >= 60 && elapsedTime < 80) {
  let progress = Math.min((elapsedTime - 60) / 20, 1);
  
  // Animate the electron density map with contour lines and shading
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.electronDensity) {
    window.Ethene_Molecule_Animation.electronDensity.material.opacity = progress;
    // Optionally adjust transparency to emphasize semi-transparent shading
    window.Ethene_Molecule_Animation.electronDensity.material.transparent = true;
    window.Ethene_Molecule_Animation.electronDensity.material.opacity = progress;
  }
  
  // Optionally adjust the camera slightly to focus on the electron density overlay
  if (window.camera) {
    // A small pan to highlight electron cloud distribution
    window.camera.position.x = THREE.MathUtils.lerp(0, 1, progress);
  }
}

// ----- TIMECODE 01:20 - Bond Rigidity & Reactivity -----
if (elapsedTime >= 80 && elapsedTime < 100) {
  let progress = Math.min((elapsedTime - 80) / 20, 1);
  
  // Animate visual markers and labels showing bond rigidity
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.rigidityMarker) {
    window.Ethene_Molecule_Animation.rigidityMarker.material.opacity = progress;
    // Scale marker up for emphasis
    window.Ethene_Molecule_Animation.rigidityMarker.scale.setScalar(THREE.MathUtils.lerp(0.5, 1, progress));
  }
  
  // Optionally, stop or slow down rotation to indicate rigidity on the C-C axis
  if (window.Ethene_Molecule_Animation) {
    window.Ethene_Molecule_Animation.rotation.y += 0.1 * delta; // slower rotation to emphasize fixed structure
  }
  
  // Additional reactive site labels may fade in (if defined)
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.reactivityLabels) {
    window.Ethene_Molecule_Animation.reactivityLabels.material.opacity = progress;
  }
}

// ----- TIMECODE 01:40 - Recap of Key Features -----
if (elapsedTime >= 100) {
  // Gradually zoom out the camera back to the original view
  if (window.camera) {
    window.camera.position.z = THREE.MathUtils.lerp(window.camera.position.z, 10, 0.02);
    // Optionally, gently re-center the camera
    window.camera.position.x = THREE.MathUtils.lerp(window.camera.position.x, 0, 0.02);
    window.camera.position.y = THREE.MathUtils.lerp(window.camera.position.y, 0, 0.02);
  }
  
  // Slowly rotate the entire ethene molecule for a comprehensive view
  if (window.Ethene_Molecule_Animation) {
    window.Ethene_Molecule_Animation.rotation.y += 0.2 * delta;
  }
  
  // Optionally, bring back any overview annotations or summaries (if defined)
  if (window.Ethene_Molecule_Animation && window.Ethene_Molecule_Animation.overviewAnnotations) {
    window.Ethene_Molecule_Animation.overviewAnnotations.material.opacity = THREE.MathUtils.lerp(window.Ethene_Molecule_Animation.overviewAnnotations.material.opacity || 0, 1, 0.02);
  }
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
