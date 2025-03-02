// Scientific visualization: Visualizing Methanol's 3D Structure
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
// Atom materials
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x808080 }); // gray
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, emissive: 0xff6666 }); // red with gentle glow
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff }); // white
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, opacity: 0.8, transparent: true });
const electronPairMaterial = new THREE.MeshPhongMaterial({ color: 0x66ccff, opacity: 0.6, transparent: true });

// Atom geometries (spheres)
const carbonGeometry = new THREE.SphereGeometry(0.5, 32, 32);
const oxygenGeometry = new THREE.SphereGeometry(0.55, 32, 32);
const hydrogenGeometry = new THREE.SphereGeometry(0.3, 32, 32);
const electronPairGeometry = new THREE.SphereGeometry(0.15, 16, 16);

// Create a group for the Methanol molecule
const methanol = new THREE.Group();
window.MethanolMolecule = methanol;

// Define bond lengths (in arbitrary units)
const r_CO = 1.0;  // Carbon-Oxygen bond length
const r_CH = 0.8;  // Carbon-Hydrogen bond length

// Place central Carbon atom at origin
const carbon = new THREE.Mesh(carbonGeometry, carbonMaterial);
carbon.position.set(0, 0, 0);
methanol.add(carbon);

// Position Oxygen along +Z direction from carbon
const oxygen = new THREE.Mesh(oxygenGeometry, oxygenMaterial);
oxygen.position.set(0, 0, r_CO);
methanol.add(oxygen);

// Tetrahedral geometry for CH3 hydrogens
// Standard tetrahedral vertices (with one vertex along +Z already taken by oxygen)
// Using normalized vectors for the remaining three vertices: 
// v2 = (  (2√2)/3,    0,   -1/3 )
// v3 = ( -√2/3,   √(6)/3,  -1/3 )
// v4 = ( -√2/3,  -√(6)/3,  -1/3 )
const hPositions = [];

const v2 = new THREE.Vector3((2 * Math.sqrt(2)) / 3, 0, -1 / 3).normalize().multiplyScalar(r_CH);
const v3 = new THREE.Vector3(-Math.sqrt(2) / 3, Math.sqrt(6) / 3, -1 / 3).normalize().multiplyScalar(r_CH);
const v4 = new THREE.Vector3(-Math.sqrt(2) / 3, -Math.sqrt(6) / 3, -1 / 3).normalize().multiplyScalar(r_CH);

hPositions.push(v2, v3, v4);

// Create and add the three Hydrogen atoms (CH3 group)
const hydrogens = [];
hPositions.forEach(pos => {
  const hAtom = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);
  hAtom.position.copy(pos);
  hydrogens.push(hAtom);
  methanol.add(hAtom);
});

// Create a helper function to generate bonds between two points
function createBond(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  // Cylinder geometry aligned along y-axis
  const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 12);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);
  // Position bond at midpoint
  bond.position.copy(start).lerp(end, 0.5);
  // Align the bond with the direction vector
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize());
  return bond;
}

// Create bonds between Carbon and each Hydrogen (CH3 bonds)
const bonds = [];
hydrogens.forEach(hAtom => {
  bonds.push(createBond(carbon.position, hAtom.position));
});

// Create bond between Carbon and Oxygen
bonds.push(createBond(carbon.position, oxygen.position));

// Add all bonds to the group
bonds.forEach(bond => methanol.add(bond));

// Create representations for oxygen's two lone electron pairs.
// These are positioned in the plane perpendicular to the C-O bond.
// Here, we offset them in the X direction from the oxygen.
const lpOffset = 0.3;
const lp1 = new THREE.Mesh(electronPairGeometry, electronPairMaterial);
const lp2 = new THREE.Mesh(electronPairGeometry, electronPairMaterial);
// Position lone pairs relative to oxygen; placed in oxygen's local XY-plane
lp1.position.set( lpOffset, 0.15, oxygen.position.z );
lp2.position.set(-lpOffset, 0.15, oxygen.position.z );
methanol.add(lp1, lp2);

// Create a semi-transparent guide overlay (wireframe sphere) to accentuate spatial orientation
const overlayGeometry = new THREE.SphereGeometry(1.4, 32, 32);
const overlayMaterial = new THREE.MeshBasicMaterial({
  color: 0xaaaaaa,
  wireframe: true,
  opacity: 0.3,
  transparent: true
});
const guideOverlay = new THREE.Mesh(overlayGeometry, overlayMaterial);
guideOverlay.position.set(0, 0, 0);
methanol.add(guideOverlay);

// (Optional) Attach basic animation data in userData for external animation handlers
methanol.userData = {
  rotationSpeed: 0.005,  // For smooth 360° orbit animation
  pulsateBonds: true,    // Flag to indicate bonds should pulsate (electron density effect)
  zoomFocus: new THREE.Vector3(0, 0, 0), // Center for zooming on the CH3 group (could be carbon's position)
  annotationText: "CH3OH\nTetrahedral carbon center\nPolar hydroxyl group"
};

// Add the completed methanol molecule group to the scene
scene.add(methanol);
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

// Retrieve the methanol molecule object
const methanol = scene.getObjectByName("MethanolMolecule");

// ────────────────────────────────────────────────────────────
// At 00:00 - 00:20: Reveal Methanol Molecule with labels and formula
if (elapsedTime < 20) {
  if (methanol) {
    // If a provided update function exists, use it (with "intro" mode)
    if (typeof window.updateMethanolMolecule === 'function') {
      window.updateMethanolMolecule(elapsedTime, { mode: "intro" });
    } else {
      // Fallback: slowly reveal the molecule by fading in its opacity
      methanol.traverse(child => {
        if (child.isMesh && child.material && child.material.opacity !== undefined) {
          child.material.opacity = Math.min(1, elapsedTime / 20);
          child.material.transparent = (child.material.opacity < 1);
        }
      });
      // Subtly rotate for a dynamic entrance
      methanol.rotation.y += 0.1 * deltaTime;
    }
  }
  // Optional: Fade in ambient light by adjusting scene background brightness if it's a color
  if (scene.background && scene.background.isColor) {
    // transitioning from black to a dark bluish hue
    const intensity = Math.min(1, elapsedTime / 20);
    scene.background.setRGB(0.02 * intensity, 0.02 * intensity, 0.08 * intensity);
  }
  // Position camera for introduction (initial view)
  if (camera) {
    const introPos = new THREE.Vector3(0, 0, 8);
    camera.position.lerp(introPos, 0.05);
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "Methanol: CH3OH structure" appears at bottom via UI overlay)
}

// ────────────────────────────────────────────────────────────
// At 00:20 - 00:40: Zoom in on the central carbon and display tetrahedral bond angles
else if (elapsedTime >= 20 && elapsedTime < 40) {
  if (methanol) {
    if (typeof window.updateMethanolMolecule === 'function') {
      window.updateMethanolMolecule(elapsedTime, { mode: "zoomCarbon" });
    } else {
      // Fallback: slowly move camera closer towards the molecule center (assumed near carbon)
      if (camera) {
        const zoomTarget = new THREE.Vector3(0, 0, 5);
        camera.position.lerp(zoomTarget, 0.05);
        camera.lookAt(0, 0, 0);
      }
      // Optionally highlight bonds by modulating their material opacity
      methanol.traverse(child => {
        if (child.isMesh && child.name.toLowerCase().includes("bond")) {
          child.material.opacity = 0.5 + 0.5 * Math.abs(Math.sin((elapsedTime - 20) * 3));
          child.material.transparent = true;
        }
      });
    }
  }
  // (Caption: "Tetrahedral carbon center")
}

// ────────────────────────────────────────────────────────────
// At 00:40 - 01:00: Focus on the methyl (CH3) group and animate hydrogens with pulsating bonds
else if (elapsedTime >= 40 && elapsedTime < 60) {
  if (methanol) {
    if (typeof window.updateMethanolMolecule === 'function') {
      window.updateMethanolMolecule(elapsedTime, { mode: "animateCH3" });
    } else {
      // Fallback: animate the three hydrogens (assumed to have names containing "Hydrogen" but not "OH")
      methanol.traverse(child => {
        if (child.isMesh && child.name.toLowerCase().includes("hydrogen") &&
           !child.name.toLowerCase().includes("oh")) {
          const pulsate = 1 + 0.1 * Math.sin((elapsedTime - 40) * 5);
          child.scale.set(pulsate, pulsate, pulsate);
        }
        // Also modulate bond thickness or emissive intensity as a simple pulsation effect
        if (child.isMesh && child.name.toLowerCase().includes("bond") &&
           child.name.toLowerCase().includes("ch3")) {
          child.material.emissive = new THREE.Color(0xffffff);
          child.material.emissiveIntensity = 0.5 + 0.5 * Math.abs(Math.sin((elapsedTime - 40) * 5));
        }
      });
    }
  }
  // (Caption: "Visualizing the CH3 group")
}

// ────────────────────────────────────────────────────────────
// At 01:00 - 01:20: Transition to hydroxyl (-OH) group; highlight oxygen glow & lone pairs
else if (elapsedTime >= 60 && elapsedTime < 80) {
  if (methanol) {
    if (typeof window.updateMethanolMolecule === 'function') {
      window.updateMethanolMolecule(elapsedTime, { mode: "hydroxyl" });
    } else {
      // Fallback: search for the oxygen atom in the hydroxyl group by name and apply a glowing effect
      methanol.traverse(child => {
        if (child.isMesh && child.name.toLowerCase().includes("oxygen")) {
          // Set a glowing emissive color
          child.material.emissive = new THREE.Color(0xff0000);
          child.material.emissiveIntensity = 1 + 0.3 * Math.sin((elapsedTime - 60) * 4);
        }
        // If there is a bond connecting the oxygen, add pulsation to indicate polarity
        if (child.isMesh && child.name.toLowerCase().includes("oh_bond")) {
          child.material.emissive = new THREE.Color(0xffaaaa);
          child.material.emissiveIntensity = 0.5 + 0.5 * Math.abs(Math.sin((elapsedTime - 60) * 4));
        }
      });
    }
  }
  // (Caption: "Introducing the -OH group")
}

// ────────────────────────────────────────────────────────────
// At 01:20 - 01:40: Begin smooth 360° rotation of the entire methanol molecule
else if (elapsedTime >= 80 && elapsedTime < 100) {
  if (methanol) {
    if (typeof window.updateMethanolMolecule === 'function') {
      window.updateMethanolMolecule(elapsedTime, { mode: "rotate" });
    } else {
      // Fallback: rotate the molecule slowly for 360° over 20 seconds
      const rotationSpeed = (2 * Math.PI) / 20; // full rotation in 20 seconds
      methanol.rotation.y += rotationSpeed * deltaTime;
    }
  }
  // Orbit the camera around the molecule (if camera is available)
  if (camera) {
    const orbitRadius = 5;
    const angle = ((elapsedTime - 80) / 20) * (2 * Math.PI);
    camera.position.x = orbitRadius * Math.cos(angle);
    camera.position.z = orbitRadius * Math.sin(angle);
    camera.position.y = 2; // slight elevation for a better view
    camera.lookAt(0, 0, 0);
  }
  // (Caption: "3D rotation reveals shape")
}

// ────────────────────────────────────────────────────────────
// At 01:40 - 02:00: Freeze rotation and overlay scientific annotations that gradually fade out
else if (elapsedTime >= 100 && elapsedTime < 120) {
  if (methanol) {
    // Freeze the rotation by not updating it further (or use an update function mode "freeze")
    if (typeof window.updateMethanolMolecule === 'function') {
      window.updateMethanolMolecule(elapsedTime, { mode: "freeze" });
    }
  }
  // Keep camera static at its last orbit position
  if (camera) {
    // Optionally, slowly adjust camera to a final predetermined position for a neat framing
    const finalCamPos = new THREE.Vector3(4, 2, 4);
    camera.position.lerp(finalCamPos, 0.02);
    camera.lookAt(0, 0, 0);
  }
  // Fading out any on-screen annotations or overlays handled by a UI layer
  if (typeof window.updateAnnotationsOpacity === "function") {
    // Assume this function gradually fades out annotations
    window.updateAnnotationsOpacity(1 - ((elapsedTime - 100) / 20));
  }
  // (Caption: "Summary: bonds & polarity" fades out)
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
