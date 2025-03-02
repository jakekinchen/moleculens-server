// Scientific visualization: 3D Exploration of Acetic Acid Geometry
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
// Atom and bond materials
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x808080, shininess: 80 }); // grey, semi-glossy
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 80 });  // red, semi-glossy
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 80 }); // white, semi-glossy
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999 }); // matte for bonds

// Dashed overlay material for bond connectivity and angle indicators
const dashedMaterial = new THREE.LineDashedMaterial({ color: 0xffff00, dashSize: 0.15, gapSize: 0.1 });

// Atom geometry (base sphere)
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create helper functions for bonds
function createBond(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();

  // Create a cylinder that connects the two atoms
  const bondGeom = new THREE.CylinderGeometry(0.1, 0.1, length, 12);
  const bond = new THREE.Mesh(bondGeom, bondMaterial);

  // Position bond: move to midpoint between start and end
  bond.position.copy(start).lerp(end, 0.5);

  // Align the cylinder to connect start and end.
  bond.quaternion.setFromUnitVectors(
    new THREE.Vector3(0, 1, 0), 
    direction.clone().normalize()
  );

  return bond;
}

function createDashedBond(start, end) {
  const points = [];
  points.push(start);
  points.push(end);
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  // computeLineDistances() is required for dashed lines
  geometry.computeLineDistances();
  const line = new THREE.Line(geometry, dashedMaterial);
  return line;
}

// Create a group to hold the whole acetic acid molecule
const aceticAcid = new THREE.Group();
window.aceticAcid = aceticAcid;

// ---------------------
// Define atomic positions
// ---------------------

// CH3 group (methyl)
const C1_pos = new THREE.Vector3(0, 0, 0);
const H1_pos = new THREE.Vector3(-0.7, 0.7, 0);
const H2_pos = new THREE.Vector3(-0.7, -0.7, 0);
const H3_pos = new THREE.Vector3(0, 0, 1);

// Carboxyl group (COOH)
// Position the carboxyl carbon a bit to the right of the CH3 carbon:
const C2_pos = new THREE.Vector3(2, 0, 0);
// Double bonded oxygen positioned upward (planar arrangement)
const O_dbl_pos = new THREE.Vector3(3.2, 0.8, 0);
// Hydroxyl oxygen positioned downward (with planarity)
const O_OH_pos = new THREE.Vector3(3.2, -0.8, 0);
// Hydrogen on the hydroxyl group
const HO_pos = new THREE.Vector3(3.2, -1.8, 0);

// ---------------------
// Create atom meshes
// ---------------------

// Methyl carbon (C1)
const C1 = new THREE.Mesh(atomGeometry, carbonMaterial);
C1.position.copy(C1_pos);
C1.scale.set(0.5, 0.5, 0.5);

// Methyl hydrogens
const H1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
H1.position.copy(H1_pos);
H1.scale.set(0.3, 0.3, 0.3);

const H2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
H2.position.copy(H2_pos);
H2.scale.set(0.3, 0.3, 0.3);

const H3 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
H3.position.copy(H3_pos);
H3.scale.set(0.3, 0.3, 0.3);

// Carboxyl carbon (C2)
const C2 = new THREE.Mesh(atomGeometry, carbonMaterial);
C2.position.copy(C2_pos);
C2.scale.set(0.5, 0.5, 0.5);

// Oxygen atoms
const O_dbl = new THREE.Mesh(atomGeometry, oxygenMaterial);
O_dbl.position.copy(O_dbl_pos);
O_dbl.scale.set(0.55, 0.55, 0.55);

const O_OH = new THREE.Mesh(atomGeometry, oxygenMaterial);
O_OH.position.copy(O_OH_pos);
O_OH.scale.set(0.55, 0.55, 0.55);

// Hydroxyl hydrogen
const HO = new THREE.Mesh(atomGeometry, hydrogenMaterial);
HO.position.copy(HO_pos);
HO.scale.set(0.3, 0.3, 0.3);

// ---------------------
// Create bonds (cylinders and dashed overlays)
// ---------------------
const bonds = [];

// Bond between CH3 carbon (C1) and carboxyl carbon (C2)
bonds.push(createBond(C1_pos, C2_pos));
bonds.push(createDashedBond(C1_pos, C2_pos));

// Bonds for CH3 hydrogens
bonds.push(createBond(C1_pos, H1_pos));
bonds.push(createDashedBond(C1_pos, H1_pos));

bonds.push(createBond(C1_pos, H2_pos));
bonds.push(createDashedBond(C1_pos, H2_pos));

bonds.push(createBond(C1_pos, H3_pos));
bonds.push(createDashedBond(C1_pos, H3_pos));

// Bonds in the carboxyl group
bonds.push(createBond(C2_pos, O_dbl_pos));
bonds.push(createDashedBond(C2_pos, O_dbl_pos));

bonds.push(createBond(C2_pos, O_OH_pos));
bonds.push(createDashedBond(C2_pos, O_OH_pos));

bonds.push(createBond(O_OH_pos, HO_pos));
bonds.push(createDashedBond(O_OH_pos, HO_pos));

// ---------------------
// Assemble molecule: add atoms and bonds to the acetic acid group
// ---------------------
aceticAcid.add(C1, H1, H2, H3, C2, O_dbl, O_OH, HO);
bonds.forEach(bond => aceticAcid.add(bond));

// ---------------------
// Optional: Create overlay helper arrows to indicate bond orientations
// (These arrow helpers simulate dynamic overlays and label connectivity)
// ---------------------
function addArrowOverlay(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start).normalize();
  const length = new THREE.Vector3().subVectors(end, start).length();
  const arrow = new THREE.ArrowHelper(direction, start, length * 0.5, 0xffff00, 0.2, 0.1);
  aceticAcid.add(arrow);
}
// For demonstration, add arrow overlays on key bonds:
addArrowOverlay(C1_pos, C2_pos);
addArrowOverlay(C2_pos, O_dbl_pos);
addArrowOverlay(C2_pos, O_OH_pos);

// ---------------------
// Apply subtle rotations for functional group distinctions
// (This sets the initial orientation; animations may override these transforms)
const methylGroup = new THREE.Group();
methylGroup.add(C1, H1, H2, H3);
// Slight rotation to emphasize tetrahedral geometry of the methyl group
methylGroup.rotation.y = 0.2;

const carboxylGroup = new THREE.Group();
carboxylGroup.add(C2, O_dbl, O_OH, HO);
// Slight independent rotation for the carboxyl group
carboxylGroup.rotation.z = -0.1;

// Remove the directly added C1 and C2 from aceticAcid and reassemble group
aceticAcid.remove(C1, H1, H2, H3, C2, O_dbl, O_OH, HO);
aceticAcid.add(methylGroup, carboxylGroup);
// Bonds and overlay arrows remain at the aceticAcid level

// ---------------------
// Final scaling of the entire molecule according to scene context (approx. 10 units diameter)
aceticAcid.scale.set(2, 2, 2);

// Add the acetic acid molecule to the global scene
scene.add(aceticAcid);

// Optional: Set initial gentle rotation animation (to be updated in your render loop)
// Example: window.aceticAcid.rotation.y += 0.005; (Place this in your animation loop)
//
// Optional: Additional animated overlays such as dynamic electron density contours can be added
// using custom shaders or texture animations as separate objects and then appended to window.aceticAcid.
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

// Find the acetic acid molecule object via scene.getObjectByName()
const aceticAcid = scene.getObjectByName("aceticAcid");

// ────────────────────────────────────────────────────────────
// At 00:00 - Introducing Acetic Acid
if (elapsedTime >= 0 && elapsedTime < 20) {
  if (aceticAcid) {
    // If the geometry provides an update function, call it with phase "intro"
    if (typeof window.updateAceticAcid === 'function') {
      window.updateAceticAcid(elapsedTime, 'intro');
    } else {
      // Fallback: gentle rotation
      aceticAcid.rotation.y += 0.2 * deltaTime;
    }
  }
  
  // Camera: slow zoom in toward the molecule and center it
  if (camera) {
    const targetPos = new THREE.Vector3(0, 0, 10);
    camera.position.lerp(targetPos, 0.05);
    camera.lookAt(aceticAcid ? aceticAcid.position : new THREE.Vector3(0, 0, 0));
  }
}

// ────────────────────────────────────────────────────────────
// At 00:20 - Atomic Structure & Bonds
else if (elapsedTime >= 20 && elapsedTime < 40) {
  if (aceticAcid) {
    if (typeof window.updateAceticAcid === 'function') {
      window.updateAceticAcid(elapsedTime, 'atomic');
    } else {
      // Fallback: moderate rotation and highlight atoms via emissive pulse
      aceticAcid.rotation.y += 0.15 * deltaTime;
      aceticAcid.traverse(child => {
        if (child.isMesh && child.name.toLowerCase().includes("atom")) {
          child.material.emissive = new THREE.Color(0x660000);
          child.material.emissiveIntensity = 0.5 + 0.5 * Math.abs(Math.sin(elapsedTime * 2));
        }
        // Optionally, you can add dashed-line or arrow effects if available
      });
    }
  }
  
  // Camera: subtle shift to focus on atomic details
  if (camera) {
    const targetPos = new THREE.Vector3(2, 0, 10);
    camera.position.lerp(targetPos, 0.05);
    camera.lookAt(aceticAcid ? aceticAcid.position : new THREE.Vector3(0, 0, 0));
  }
}

// ────────────────────────────────────────────────────────────
// At 00:40 - Methyl vs. Carboxyl Groups
else if (elapsedTime >= 40 && elapsedTime < 60) {
  if (aceticAcid) {
    if (typeof window.updateAceticAcid === 'function') {
      window.updateAceticAcid(elapsedTime, 'groups');
    } else {
      // Fallback: slight rotation; if functional groups exist as children, animate them separately
      aceticAcid.rotation.y += 0.1 * deltaTime;
      aceticAcid.traverse(child => {
        if (child.isMesh) {
          const nameLC = child.name.toLowerCase();
          if (nameLC.includes("methyl")) {
            // Separate methyl group slightly using sinusoidal motion
            child.position.x = Math.sin(elapsedTime * 0.5) * 0.2;
          } else if (nameLC.includes("carboxyl")) {
            // Separate carboxyl group in opposite phase
            child.position.x = Math.cos(elapsedTime * 0.5) * 0.2;
          }
        }
      });
    }
  }
  
  // Camera: adjust viewpoint to better showcase group separation
  if (camera) {
    const targetPos = new THREE.Vector3(0, 1, 8);
    camera.position.lerp(targetPos, 0.05);
    camera.lookAt(aceticAcid ? aceticAcid.position : new THREE.Vector3(0, 0, 0));
  }
}

// ────────────────────────────────────────────────────────────
// At 01:00 - Bond Angles & Geometry
else if (elapsedTime >= 60 && elapsedTime < 90) {
  if (aceticAcid) {
    if (typeof window.updateAceticAcid === 'function') {
      window.updateAceticAcid(elapsedTime, 'bondAngles');
    } else {
      // Fallback: slow rotation and fade in any bond-angle overlays
      aceticAcid.rotation.y += 0.05 * deltaTime;
      aceticAcid.traverse(child => {
        if (child.isMesh && child.name.toLowerCase().includes("angle")) {
          // Compute fade-in progress over 5 seconds
          const opacity = Math.min(1, (elapsedTime - 60) / 5);
          if(child.material) {
            child.material.opacity = opacity;
            child.material.transparent = true;
          }
        }
      });
    }
  }
  
  // Camera: pan along the molecular backbone for clearer angle visualization
  if (camera) {
    const targetPos = new THREE.Vector3(0, 0, 7);
    camera.position.lerp(targetPos, 0.05);
    camera.lookAt(aceticAcid ? aceticAcid.position : new THREE.Vector3(0, 0, 0));
  }
}

// ────────────────────────────────────────────────────────────
// At 01:30 - Electron Density Mapping
else if (elapsedTime >= 90 && elapsedTime < 110) {
  if (aceticAcid) {
    if (typeof window.updateAceticAcid === 'function') {
      window.updateAceticAcid(elapsedTime, 'orbitals');
    } else {
      // Fallback: slight additional rotation and animate electron density contours
      aceticAcid.rotation.y += 0.05 * deltaTime;
      aceticAcid.traverse(child => {
        if (child.isMesh && child.name.toLowerCase().includes("electron")) {
          // Fade in electron cloud overlays over 5 seconds
          const opacity = Math.min(1, (elapsedTime - 90) / 5);
          if (child.material) {
            child.material.opacity = opacity;
            child.material.transparent = true;
          }
        }
      });
    }
  }
}

// ────────────────────────────────────────────────────────────
// At 01:50 - Summary: 3D Molecular Shape
else if (elapsedTime >= 110) {
  if (aceticAcid) {
    if (typeof window.updateAceticAcid === 'function') {
      window.updateAceticAcid(elapsedTime, 'summary');
    } else {
      // Fallback: slow rotation and gradual fade-out of highlighted overlays
      aceticAcid.rotation.y += 0.02 * deltaTime;
      aceticAcid.traverse(child => {
        if (child.isMesh && child.material && child.material.opacity !== undefined) {
          // Fade out elements smoothly over 10 seconds
          child.material.opacity = Math.max(0, 1 - (elapsedTime - 110) / 10);
          child.material.transparent = true;
        }
      });
    }
  }
  
  // Camera: pull back to provide an overall view as highlighted elements disappear
  if (camera) {
    const targetPos = new THREE.Vector3(0, 0, 12);
    camera.position.lerp(targetPos, 0.05);
    camera.lookAt(aceticAcid ? aceticAcid.position : new THREE.Vector3(0, 0, 0));
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
