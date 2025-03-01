// Scientific visualization: Exploring Acetic Acid Geometry
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
// Materials for atoms and bonds
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x333333, shininess: 100 });
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 100 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 100 });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x99ccff, emissive: 0x222244, shininess: 150 });

// Atom geometry (spheres)
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create a group for the acetic acid molecule
const aceticMolecule = new THREE.Group();
window.aceticMolecule = aceticMolecule;

// Positions for acetic acid (CH3COOH)
// CH3 (methyl) group
const c1Position = new THREE.Vector3(0, 0, 0);
const h1Position = new THREE.Vector3(-0.8, 0.8, 0);
const h2Position = new THREE.Vector3(-0.8, -0.8, 0);
const h3Position = new THREE.Vector3(0, 0, -1.0);

// Carboxyl group (COOH)
const c2Position = new THREE.Vector3(1.5, 0, 0);
const o1Position = new THREE.Vector3(2.8, 0.6, 0);  // Double-bonded oxygen
const o2Position = new THREE.Vector3(2.8, -0.6, 0); // Single-bonded oxygen
const h4Position = new THREE.Vector3(3.6, -0.6, 0); // Hydrogen bonded to O2

// Create atoms
// Methyl carbon
const c1 = new THREE.Mesh(atomGeometry, carbonMaterial);
c1.position.copy(c1Position);
c1.scale.set(0.5, 0.5, 0.5);

// Methyl hydrogens
const h1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h1.position.copy(h1Position);
h1.scale.set(0.3, 0.3, 0.3);

const h2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h2.position.copy(h2Position);
h2.scale.set(0.3, 0.3, 0.3);

const h3 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h3.position.copy(h3Position);
h3.scale.set(0.3, 0.3, 0.3);

// Carboxyl carbon
const c2 = new THREE.Mesh(atomGeometry, carbonMaterial);
c2.position.copy(c2Position);
c2.scale.set(0.5, 0.5, 0.5);

// Carboxyl oxygens
const o1 = new THREE.Mesh(atomGeometry, oxygenMaterial);
o1.position.copy(o1Position);
o1.scale.set(0.55, 0.55, 0.55);

const o2 = new THREE.Mesh(atomGeometry, oxygenMaterial);
o2.position.copy(o2Position);
o2.scale.set(0.55, 0.55, 0.55);

// Hydrogen bonded to O2
const h4 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h4.position.copy(h4Position);
h4.scale.set(0.3, 0.3, 0.3);

// Function to create bonds between two points
function createBond(start, end) {
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();

    const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
    const bond = new THREE.Mesh(bondGeometry, bondMaterial);
    
    // Position bond at midpoint
    bond.position.copy(start).lerp(end, 0.5);
    
    // Align cylinder with the direction vector
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
    
    return bond;
}

// Create bonds for the molecule
const bonds = [];

// Bonds within methyl group
bonds.push(createBond(c1.position, h1.position));
bonds.push(createBond(c1.position, h2.position));
bonds.push(createBond(c1.position, h3.position));

// Bond between methyl carbon and carboxyl carbon
bonds.push(createBond(c1.position, c2.position));

// Bonds in carboxyl group
bonds.push(createBond(c2.position, o1.position)); // double-bonded oxygen (displayed as a bond)
bonds.push(createBond(c2.position, o2.position)); // single-bonded oxygen
bonds.push(createBond(o2.position, h4.position)); // O2-H bond

// Add atoms to the acetic molecule group
aceticMolecule.add(c1, h1, h2, h3, c2, o1, o2, h4);

// Add bonds to the molecule group
bonds.forEach(bond => aceticMolecule.add(bond));

// Optional: Add an arrow helper to highlight the planar carboxyl group (COOH)
// Arrow pointing from carboxyl carbon (c2) towards the center between oxygens
const carboxylCenter = new THREE.Vector3().addVectors(o1.position, o2.position).multiplyScalar(0.5);
const arrowDir = new THREE.Vector3().subVectors(carboxylCenter, c2.position).normalize();
const arrowLength = 1.0;
const arrowColor = 0xffff00;
const arrowHelper = new THREE.ArrowHelper(arrowDir, c2.position, arrowLength, arrowColor, 0.2, 0.1);
aceticMolecule.add(arrowHelper);

// Optional: Add bond angle marker (using a simple line for demonstration)
// This line connects the two oxygens to visually emphasize the COOH group's planar nature.
const markerMaterial = new THREE.LineBasicMaterial({ color: 0xffffff });
const markerGeometry = new THREE.BufferGeometry().setFromPoints([o1.position, o2.position]);
const bondAngleMarker = new THREE.Line(markerGeometry, markerMaterial);
aceticMolecule.add(bondAngleMarker);

// Add the completed acetic acid molecule to the scene
scene.add(aceticMolecule);

// GeometryAgent LLM-generated code
// Hybridization_Comparison_Overlay

// Create a parent group for the overlay and store it globally
const hybridOverlay = new THREE.Group();
window.Hybridization_Comparison_Overlay = hybridOverlay;

// ---------- LEFT PANEL: Carboxyl Group (sp2 hybridization) ----------

// Create a canvas texture for a blue-green gradient background
const canvas = document.createElement('canvas');
canvas.width = 256;
canvas.height = 256;
const ctx = canvas.getContext('2d');
const gradient = ctx.createLinearGradient(0, 0, canvas.width, canvas.height);
gradient.addColorStop(0, '#0000ff'); // blue
gradient.addColorStop(1, '#00ff00'); // green
ctx.fillStyle = gradient;
ctx.fillRect(0, 0, canvas.width, canvas.height);
const gradientTexture = new THREE.CanvasTexture(canvas);

// Create material using the gradient texture
const electronMaterial = new THREE.MeshBasicMaterial({ 
  map: gradientTexture, 
  side: THREE.DoubleSide,
  transparent: true,
  opacity: 0.8
});

// Create a plane for the left panel background
const panelWidth = 10, panelHeight = 20;
const leftPanelGeometry = new THREE.PlaneGeometry(panelWidth, panelHeight);
const leftPanelMesh = new THREE.Mesh(leftPanelGeometry, electronMaterial);
// Position the left panel to the left half of the scene
leftPanelMesh.position.set(-panelWidth * 0.5, 0, 0);
 
// Create a left panel group and add the background
const leftPanelGroup = new THREE.Group();
leftPanelGroup.add(leftPanelMesh);

// Add animated molecular orbitals overlay (semicircular torus shapes)
const orbitalMaterial = new THREE.MeshPhongMaterial({ 
  color: 0xffffff, 
  emissive: 0x0000ff, 
  transparent: true, 
  opacity: 0.5 
});
for (let i = 0; i < 3; i++) {
  // Create a torus to simulate an orbital path
  const torus = new THREE.Mesh(
    new THREE.TorusGeometry(3 - i * 0.5, 0.1, 16, 100),
    orbitalMaterial
  );
  torus.rotation.x = Math.PI / 2;
  // Slightly offset each orbital for a layered look
  torus.position.z = i * 0.2;
  leftPanelGroup.add(torus);
}

// Store left panel group for potential interactive animation updates
window.leftHybridPanel = leftPanelGroup;

// ---------- RIGHT PANEL: Methyl Group (sp3 hybridization) ----------

// Define warm-tone material with transparency for tetrahedral geometry
const tetraMaterial = new THREE.MeshPhongMaterial({
  color: 0xcc6600,
  transparent: true,
  opacity: 0.6
});

// Create a group for the right panel
const rightPanelGroup = new THREE.Group();

// Create a background plane for the right panel with a warm transparent tone
const rightBgMaterial = new THREE.MeshBasicMaterial({
  color: 0xffcc99,
  side: THREE.DoubleSide,
  transparent: true,
  opacity: 0.3
});
const rightPanelGeometry = new THREE.PlaneGeometry(panelWidth, panelHeight);
const rightPanelMesh = new THREE.Mesh(rightPanelGeometry, rightBgMaterial);
rightPanelMesh.position.set(panelWidth * 0.5, 0, 0);
rightPanelGroup.add(rightPanelMesh);

// Create a subgroup for the tetrahedral (methyl) molecule overlay
const tetraGroup = new THREE.Group();

// Central atom (carbon-like)
const centerAtomGeometry = new THREE.SphereGeometry(0.6, 32, 32);
const centerAtom = new THREE.Mesh(centerAtomGeometry, tetraMaterial);
tetraGroup.add(centerAtom);

// Create peripheral atoms (simulating the hydrogen atoms in a tetrahedral arrangement)
// Tetrahedral vertex directions (normalized)
const directions = [
  new THREE.Vector3( 1,  1,  1).normalize(),
  new THREE.Vector3(-1, -1,  1).normalize(),
  new THREE.Vector3(-1,  1, -1).normalize(),
  new THREE.Vector3( 1, -1, -1).normalize()
];

const peripheralAtoms = [];
const peripheralAtomGeometry = new THREE.SphereGeometry(0.4, 32, 32);
for (let i = 0; i < directions.length; i++) {
  const atom = new THREE.Mesh(peripheralAtomGeometry, tetraMaterial);
  // Position atoms at a fixed bond length from the center (e.g., 2.0 units)
  atom.position.copy(directions[i]).multiplyScalar(2.0);
  peripheralAtoms.push(atom);
  tetraGroup.add(atom);
}

// Create bright yellow overlay lines for bond angle visualizations
const overlayLineMaterial = new THREE.MeshBasicMaterial({ color: 0xffff00 });

// Function to create a bond line between two points
function createBondLine(start, end) {
  const bondDir = new THREE.Vector3().subVectors(end, start);
  const bondLength = bondDir.length();
  const bondGeometry = new THREE.CylinderGeometry(0.05, 0.05, bondLength, 8);
  const bond = new THREE.Mesh(bondGeometry, overlayLineMaterial);
  // Position bond at midpoint
  bond.position.copy(start).lerp(end, 0.5);
  // Orient the cylinder along the bond direction
  bond.quaternion.setFromUnitVectors(
      new THREE.Vector3(0, 1, 0),
      bondDir.normalize()
  );
  return bond;
}

// Create bonds from the center to each peripheral atom
for (let i = 0; i < peripheralAtoms.length; i++) {
  const bond = createBondLine(centerAtom.position, peripheralAtoms[i].position);
  tetraGroup.add(bond);
}

// Add the tetrahedral molecule overlay to the right panel group
rightPanelGroup.add(tetraGroup);

// Position the tetraGroup slightly forward for visual clarity
tetraGroup.position.set(panelWidth * 0.5, 0, 0.5);

// Store right panel group for potential interactive animation updates
window.rightHybridPanel = rightPanelGroup;

// ---------- Assemble Overall Overlay ----------

// Position left and right panels appropriately in the scene (split-screen)
// Assuming the scene's center is 0, offset panels along the x-axis.
const offsetX = panelWidth; // adjust spacing as needed
leftPanelGroup.position.set(-offsetX * 0.5, 0, 0);
rightPanelGroup.position.set(offsetX * 0.5, 0, 0);

hybridOverlay.add(leftPanelGroup);
hybridOverlay.add(rightPanelGroup);

// Save the overlay group globally and add to the scene
window.Hybridization_Comparison_Overlay = hybridOverlay;
scene.add(hybridOverlay);

// Optional: Store animation parameters in userData for smooth split-screen appearing and interactive highlighting
hybridOverlay.userData = {
  transitionEffect: 'smooth split-screen appearance',
  interactiveElements: 'highlight overlay on user focus'
};

// Note: Animation/update routines should be added in the render loop to animate orbitals and interactive effects.
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
// Get elapsed time in seconds
const elapsedTime = clock.getElapsedTime();

// ------------------------------
// 00:00 - 00:20: Initial Overview & Orbiting Rotation
// ------------------------------
if (elapsedTime < 20) {
  // Rotate the acetic acid molecule and add a subtle glow to its bonds.
  if (window.Acetic_Acid_Molecule_Model) {
    window.Acetic_Acid_Molecule_Model.rotation.y += 0.5 * deltaTime;
    // Traverse children to modulate emissive intensity for a glowing effect.
    window.Acetic_Acid_Molecule_Model.traverse(child => {
      if (child.material && child.material.emissive) {
        // Create a gentle pulsating glow effect on bonds.
        child.material.emissiveIntensity = 1 + 0.3 * Math.sin(elapsedTime * 2);
      }
    });
  }
  // Orbit the camera slowly around the molecule.
  if (window.camera) {
    const radius = 10;
    const orbitSpeed = 0.5; // radians per second
    const angle = elapsedTime * orbitSpeed;
    window.camera.position.x = radius * Math.cos(angle);
    window.camera.position.z = radius * Math.sin(angle);
    window.camera.position.y = 4; // Keep a fixed elevation
    window.camera.lookAt(0, 0, 0);
  }
}

// ------------------------------
// 00:20 - 00:40: Zoom In on the Carboxyl Group (COOH)
// ------------------------------
if (elapsedTime >= 20 && elapsedTime < 40) {
  const t = (elapsedTime - 20) / 20; // Normalized 0 to 1 over 20 seconds.
  if (window.camera) {
    // Define a hypothetical target position that zooms into the carboxyl group.
    // (Assume (5, 2, 5) approximates the desired focus area.)
    const startPos = window.camera.position.clone();
    const targetPos = new THREE.Vector3(5, 2, 5);
    window.camera.position.lerp(targetPos, t);
    // Set the focus to the carboxyl group area.
    window.camera.lookAt(new THREE.Vector3(4.5, 1.8, 4.5));
  }
  // Animate the appearance of labels and arrows highlighting key atoms.
  if (window.Acetic_Acid_Molecule_Model) {
    const carboxylHighlight = window.Acetic_Acid_Molecule_Model.getObjectByName('CarboxylHighlight');
    if (carboxylHighlight) {
      carboxylHighlight.visible = true;
      // Fade in the highlight overlay.
      if (carboxylHighlight.material) {
        carboxylHighlight.material.opacity = t;
      }
    }
  }
}

// ------------------------------
// 00:40 - 01:00: Side View Emphasizing Carboxyl Planarity
// ------------------------------
if (elapsedTime >= 40 && elapsedTime < 60) {
  const t = (elapsedTime - 40) / 20; // Normalize between 0 and 1.
  if (window.camera) {
    // Transition to a side view position.
    const currentPos = window.camera.position.clone();
    const sideViewPos = new THREE.Vector3(10, 2, 0); // Hypothetical side view location.
    window.camera.position.lerp(sideViewPos, t);
    window.camera.lookAt(new THREE.Vector3(4.5, 1.8, 4.5)); // Focus remains on the carboxyl group.
  }
  // Emphasize key bond elements with animated bond lines and angle markers.
  if (window.Acetic_Acid_Molecule_Model) {
    const doubleBond = window.Acetic_Acid_Molecule_Model.getObjectByName('DoubleBondHighlight');
    const angleMarker = window.Acetic_Acid_Molecule_Model.getObjectByName('AngleMarker');
    if (doubleBond) {
      doubleBond.visible = true;
      if (doubleBond.material) {
        doubleBond.material.opacity = t;
      }
    }
    if (angleMarker) {
      angleMarker.visible = true;
      if (angleMarker.material) {
        angleMarker.material.opacity = t;
      }
    }
  }
}

// ------------------------------
// 01:00 - 01:20: Explore Central Carbon Hybridization (sp2 vs sp3)
// ------------------------------
if (elapsedTime >= 60 && elapsedTime < 80) {
  const t = (elapsedTime - 60) / 20;
  if (window.camera) {
    // Move the camera to an intermediate position highlighting both the planar part (sp2)
    // and the tetrahedral methyl group (sp3). (Assume (7,3,7) works well.)
    const hybridViewPos = new THREE.Vector3(7, 3, 7);
    window.camera.position.lerp(hybridViewPos, t);
    window.camera.lookAt(new THREE.Vector3(5, 2, 5));
  }
  // Animate overlays or annotations contrasting sp2 and sp3 hybridization.
  if (window.Acetic_Acid_Molecule_Model) {
    const sp2Overlay = window.Acetic_Acid_Molecule_Model.getObjectByName('sp2Overlay');
    const sp3Overlay = window.Acetic_Acid_Molecule_Model.getObjectByName('sp3Overlay');
    if (sp2Overlay) {
      sp2Overlay.visible = true;
      if (sp2Overlay.material) {
        sp2Overlay.material.opacity = t;
      }
    }
    if (sp3Overlay) {
      sp3Overlay.visible = true;
      if (sp3Overlay.material) {
        sp3Overlay.material.opacity = t;
      }
    }
  }
}

// ------------------------------
// 01:20 - 01:40: Split-Screen Hybridization Comparison Overlay
// ------------------------------
if (elapsedTime >= 80 && elapsedTime < 100) {
  const t = (elapsedTime - 80) / 20;
  if (window.Hybridization_Comparison_Overlay) {
    window.Hybridization_Comparison_Overlay.visible = true;
    // Fade in the split-screen overlay by adjusting its opacity.
    if (window.Hybridization_Comparison_Overlay.material) {
      window.Hybridization_Comparison_Overlay.material.opacity = t;
    }
    // Optional: add subtle rotation to the overlay for an interactive feel.
    window.Hybridization_Comparison_Overlay.rotation.y += 0.3 * deltaTime;
  }
}

// ------------------------------
// 01:40 and Beyond: Final Summary and Pull-Back View
// ------------------------------
if (elapsedTime >= 100) {
  if (window.camera) {
    // Gradually pull the camera back to provide a comprehensive overview.
    const finalPos = new THREE.Vector3(10, 10, 10); // Final viewing position.
    window.camera.position.lerp(finalPos, 0.02); // Use a small factor for a slow pull-back.
    window.camera.lookAt(new THREE.Vector3(4.5, 1.8, 4.5));
  }
  // Animate the reappearance of the initial model with annotated bond angles.
  if (window.Acetic_Acid_Molecule_Model) {
    const bondAnnotations = window.Acetic_Acid_Molecule_Model.getObjectByName('BondAngleAnnotations');
    if (bondAnnotations) {
      bondAnnotations.visible = true;
      if (bondAnnotations.material) {
        // Slowly increase opacity for clarity.
        bondAnnotations.material.opacity = Math.min(1, bondAnnotations.material.opacity + 0.5 * deltaTime);
      }
    }
  }
  // Fade out the hybridization overlay if it is still visible.
  if (window.Hybridization_Comparison_Overlay && window.Hybridization_Comparison_Overlay.material) {
    window.Hybridization_Comparison_Overlay.material.opacity = Math.max(0, window.Hybridization_Comparison_Overlay.material.opacity - 0.5 * deltaTime);
    if (window.Hybridization_Comparison_Overlay.material.opacity <= 0.01) {
      window.Hybridization_Comparison_Overlay.visible = false;
    }
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
