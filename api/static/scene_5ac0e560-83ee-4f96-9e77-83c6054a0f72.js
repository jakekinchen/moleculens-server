// Scientific visualization: Exploring Tetrahedral Carbon Geometry
// Auto-generated code

// Global variables
// Use existing globals if defined, otherwise create new ones
var scene = window.scene || null;
var camera = window.camera || null;
var renderer = window.renderer || null;
var controls = window.controls || null;
var clock = window.clock || null;
var isPlaying = typeof window.isPlaying !== 'undefined' ? window.isPlaying : true;
// Set animation duration to 2 minutes
var ANIMATION_DURATION = 120;


// Initialize the scene
function init() {
    // Create scene if it doesn't exist
    if (!scene) {
        scene = new THREE.Scene();
        scene.background = new THREE.Color(0x111122);
    }
    
    // Create camera if it doesn't exist
    if (!camera) {
        camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
        camera.position.z = 10;
    }
    
    // Look for an existing canvas or create renderer
    if (!renderer) {
        const canvas = document.getElementById('scene-canvas');
        if (canvas) {
            renderer = new THREE.WebGLRenderer({ 
                canvas: canvas,
                antialias: true 
            });
            renderer.setSize(window.innerWidth, window.innerHeight);
            renderer.setPixelRatio(window.devicePixelRatio);
        } else {
            console.warn('No canvas element found with id "scene-canvas"');
            return; // Exit init if canvas not found
        }
    }
    
    // Create lighting
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);
    
    const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
    directionalLight.position.set(1, 1, 1).normalize();
    scene.add(directionalLight);
    
    // Add orbit controls if they don't exist
    if (!controls && renderer) {
        controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.25;
    }
    
    // Initialize clock if it doesn't exist
    if (!clock && typeof THREE !== 'undefined') {
        clock = new THREE.Clock();
    }
    
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
// Create a group for the tetrahedral carbon molecule
const tetrahedralCarbon = new THREE.Group();
window.tetrahedralCarbon = tetrahedralCarbon;

// ---------- Central Nucleus ----------
// The nucleus is a glowing dark gray sphere with subtle blue emissive highlights
const nucleusGeometry = new THREE.SphereGeometry(1, 32, 32); // Diameter = 2 units (radius = 1)
const nucleusMaterial = new THREE.MeshPhongMaterial({ 
    color: 0x333333, 
    emissive: 0x001133, 
    shininess: 100 
});
const nucleus = new THREE.Mesh(nucleusGeometry, nucleusMaterial);
tetrahedralCarbon.add(nucleus);

// ---------- Tetrahedral Directions & Constants ----------
const bondLength = 5; // as provided
// Define four tetrahedral directions (normalized)
const directions = [
    new THREE.Vector3( 1,  1,  1).normalize(),
    new THREE.Vector3( 1, -1, -1).normalize(),
    new THREE.Vector3(-1,  1, -1).normalize(),
    new THREE.Vector3(-1, -1,  1).normalize(),
];

// Pastel colors for substituents (small spheres)
const pastelColors = [0xffb6c1, 0xbfefff, 0xffdab9, 0xe6e6fa];

// Bond material: metallic silver look with a hint for animated vibration (animation to be handled externally)
const bondMaterial = new THREE.MeshPhongMaterial({ 
    color: 0xc0c0c0, 
    shininess: 80 
});

// Cylinder geometry creation helper function for bonds
function createBond(start, end) {
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();
    // Cylinder aligned along Y; adjust its radius as needed for a bond appearance
    const bondGeometry = new THREE.CylinderGeometry(0.15, 0.15, length, 16);
    const bond = new THREE.Mesh(bondGeometry, bondMaterial);
    
    // Position the bond in the middle between start and end
    bond.position.copy(start).lerp(end, 0.5);
    
    // Rotate the bond so that its Y-axis aligns with the direction vector
    bond.quaternion.setFromUnitVectors(
        new THREE.Vector3(0, 1, 0),
        direction.clone().normalize()
    );
    
    // Mark this mesh for vibration animation
    bond.userData.vibrate = true;
    
    return bond;
}

// ---------- Bonds & Substituents ----------
// For each tetrahedral direction, create a bond extending from the nucleus center and a substituent sphere at the end.
for (let i = 0; i < directions.length; i++) {
    // Calculate the end point of the bond (and the position of the substituent)
    const endPoint = directions[i].clone().multiplyScalar(bondLength);
    
    // Create and add the bond cylinder
    const bond = createBond(new THREE.Vector3(0, 0, 0), endPoint);
    tetrahedralCarbon.add(bond);
    
    // Create a small substituent sphere at the end of the bond.
    const substituentGeometry = new THREE.SphereGeometry(0.4, 16, 16);
    const substituentMaterial = new THREE.MeshPhongMaterial({ 
        color: pastelColors[i], 
        shininess: 10 
    });
    const substituent = new THREE.Mesh(substituentGeometry, substituentMaterial);
    substituent.position.copy(endPoint);
    tetrahedralCarbon.add(substituent);
}

// ---------- Outlined Tetrahedron Structure ----------
// Create an outlined tetrahedron (with vertices matching the bond endpoints).
// THREE.TetrahedronGeometry generates a tetrahedron with a given radius (here, bondLength).
const tetrahedronGeometry = new THREE.TetrahedronGeometry(bondLength, 0);
const tetraEdges = new THREE.EdgesGeometry(tetrahedronGeometry);
const tetraOutlineMaterial = new THREE.LineBasicMaterial({
    color: 0x000000, // black outline for clarity
    linewidth: 2
});
const tetraOutline = new THREE.LineSegments(tetraEdges, tetraOutlineMaterial);
tetrahedralCarbon.add(tetraOutline);

// ---------- Semi-transparent Animated Electron Cloud ----------
// Create an encompassing electron cloud with radius 6 units.
const electronCloudGeometry = new THREE.SphereGeometry(6, 32, 32);
const electronCloudMaterial = new THREE.MeshPhongMaterial({
    color: 0x00ffff, // translucent cyan
    transparent: true,
    opacity: 0.4,
    side: THREE.DoubleSide
});
const electronCloud = new THREE.Mesh(electronCloudGeometry, electronCloudMaterial);
// Flag this mesh for pulsation animation (to be implemented in the render loop)
electronCloud.userData.pulsate = true;
tetrahedralCarbon.add(electronCloud);
window.electronCloud = electronCloud;

// ---------- Final Touches ----------
// Optionally, mark the entire group for slow rotation animation (to be handled in the render loop)
tetrahedralCarbon.userData.rotateSlowly = true;

// Add the complete tetrahedral carbon molecule to the global scene
scene.add(tetrahedralCarbon);

// GeometryAgent LLM-generated code
// Custom curve to create an arc (used for bond angle markers)
// WARNING: In modern versions of THREE.js, THREE.Curve is an ES6 class and must be instantiated with 'new'.
// Always use 'new ArcCurve(...)' to create instances of this class.
class ArcCurve extends THREE.Curve {
  constructor(radius, startAngle, endAngle) {
    super();
    this.radius = radius;
    this.startAngle = startAngle;
    this.endAngle = endAngle;
  }
  
  getPoint(t) {
    const angle = this.startAngle + t * (this.endAngle - this.startAngle);
    return new THREE.Vector3(
      Math.cos(angle) * this.radius,
      Math.sin(angle) * this.radius,
      0
    );
  }
}

// Create a group for the Bond_Angle_and_Lattice_Overlay
const bondAngleLatticeOverlay = new THREE.Group();
window.Bond_Angle_and_Lattice_Overlay = bondAngleLatticeOverlay;

// ------------------------
// 1. Create Angle Markers
// ------------------------
// Define parameters
const markerRadius = 2; // proportional to scene scale
const markerThickness = 0.05; // tube radius for the marker
const angleInRad = 109.5 * Math.PI / 180; // 109.5 degrees in radians

// Material for the angle markers: bright white with subtle neon accents,
// glossy and emissive; overlay opacity set to 50%
const angleMarkerMaterial = new THREE.MeshPhongMaterial({
  color: 0xffffff,
  emissive: 0x00ffff,
  shininess: 100,
  transparent: true,
  opacity: 0.5
});

// Create a TubeGeometry for an arc curve representing the bond angle marker
const arcCurve = new ArcCurve(markerRadius, 0, angleInRad);
const tubeGeometry = new THREE.TubeGeometry(arcCurve, 32, markerThickness, 8, false);
const markerMesh1 = new THREE.Mesh(tubeGeometry, angleMarkerMaterial);

// For better spatial indication, clone and rotate the marker (displaying similar angles in different orientations)
const markerMesh2 = markerMesh1.clone();
markerMesh2.rotation.set(THREE.Math.degToRad(45), 0, 0);

const markerMesh3 = markerMesh1.clone();
markerMesh3.rotation.set(0, THREE.Math.degToRad(45), 0);

// Add the markers to the overlay group
bondAngleLatticeOverlay.add(markerMesh1, markerMesh2, markerMesh3);

// ------------------------
// 2. Create Lattice Background
// ------------------------
// Create a canvas to generate a smooth gradient texture for the lattice background
const canvas = document.createElement('canvas');
canvas.width = 256;
canvas.height = 256;
const context = canvas.getContext('2d');

// Create a vertical gradient - soft light blue gradients
const gradient = context.createLinearGradient(0, 0, 0, 256);
gradient.addColorStop(0, '#BFEFFF'); // light, soft accent
gradient.addColorStop(1, '#ADD8E6'); // light blue
context.fillStyle = gradient;
context.fillRect(0, 0, 256, 256);

// Create a texture from the canvas
const gradientTexture = new THREE.CanvasTexture(canvas);
gradientTexture.minFilter = THREE.LinearFilter;

// Material for the lattice background: translucent with smooth gradient and 50% opacity
const latticeMaterial = new THREE.MeshBasicMaterial({
  map: gradientTexture,
  side: THREE.DoubleSide,
  transparent: true,
  opacity: 0.5
});

// Create a large plane that covers the full scene background
const latticeGeometry = new THREE.PlaneGeometry(100, 100, 10, 10);
const latticeMesh = new THREE.Mesh(latticeGeometry, latticeMaterial);
// Position the lattice slightly behind the angle markers
latticeMesh.position.set(0, 0, -10);

// Optional: mark lattice for later animation (e.g., fade in/out) through userData
latticeMesh.userData.animate = true;

// Add the lattice background to the overlay group
bondAngleLatticeOverlay.add(latticeMesh);

// ------------------------
// 3. Finalize Overlay Object
// ------------------------
// Optionally, add userData to the group for animation tuning (e.g., synchronized with scene rotation)
bondAngleLatticeOverlay.userData = {
  animation: "synchronized with scene rotation; fade-in and fade-out during key moments"
};

// Add the Bond_Angle_and_Lattice_Overlay group to the scene
scene.add(bondAngleLatticeOverlay);
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
    
    // Get the adjusted time with offset applied
    const adjustedTime = updateUI();
    
    // We'll make this time variable available to the animation code
    window.animationTime = adjustedTime;
    
    // Animation created by the AnimationAgent
// ──────────────────────────────────────────────
// Get elapsed time and deltaTime from global clock/control
const elapsedTime = window.animationTime || clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Helper: fade object's material opacity over a duration (applied per mesh)
function applyFade(object, progress) {
  object.traverse(child => {
    if (child.isMesh && child.material) {
      child.material.transparent = true;
      child.material.opacity = THREE.MathUtils.clamp(progress, 0, 1);
    }
  });
}

// ──────────────────────────────────────────────
// At 00:00 - Introducing tetrahedral carbon
if (elapsedTime >= 0 && elapsedTime < 20) {
  const tetrahedralCarbon = scene.getObjectByName("tetrahedralCarbon");
  if (tetrahedralCarbon) {
    // Fade in over the first 2 seconds
    const fadeProgress = (elapsedTime < 2) ? elapsedTime / 2 : 1;
    applyFade(tetrahedralCarbon, fadeProgress);

    // Use the provided update function if available; else, perform simple rotation
    if (typeof window.updateTetrahedralCarbon === "function") {
      window.updateTetrahedralCarbon(elapsedTime);
    } else {
      tetrahedralCarbon.rotation.y += 0.1 * deltaTime;
    }
  }

  // Ensure electron cloud and bond overlay are hidden
  const electronCloud = scene.getObjectByName("electronCloud");
  if (electronCloud) {
    electronCloud.visible = false;
  }
  const bondOverlay = scene.getObjectByName("Bond_Angle_and_Lattice_Overlay");
  if (bondOverlay) {
    bondOverlay.visible = false;
  }

  // Position camera for initial view (gradually moving to z = 10)
  if (camera) {
    camera.position.lerp(new THREE.Vector3(0, 0, 10), 0.1);
    camera.lookAt(0, 0, 0);
  }
  
  // (Caption: "Introducing tetrahedral carbon")
}

// ──────────────────────────────────────────────
// At 00:20 - Sp³ hybridization explained
else if (elapsedTime >= 20 && elapsedTime < 40) {
  const tetrahedralCarbon = scene.getObjectByName("tetrahedralCarbon");
  if (tetrahedralCarbon) {
    // If available, update with parameter to highlight tetrahedral (sp³ hybridization) features
    if (typeof window.updateTetrahedralCarbon === "function") {
      window.updateTetrahedralCarbon(elapsedTime, { highlightTetrahedron: true });
    } else {
      tetrahedralCarbon.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Keep electronCloud and bond overlay hidden in this phase
  const electronCloud = scene.getObjectByName("electronCloud");
  if (electronCloud) {
    electronCloud.visible = false;
  }
  const bondOverlay = scene.getObjectByName("Bond_Angle_and_Lattice_Overlay");
  if (bondOverlay) {
    bondOverlay.visible = false;
  }
  
  // Slight camera zoom (move closer to the molecule; target z = 8)
  if (camera) {
    camera.position.lerp(new THREE.Vector3(0, 0, 8), 0.05);
    camera.lookAt(0, 0, 0);
  }
  
  // (Caption: "Sp³ hybridization explained")
}

// ──────────────────────────────────────────────
// At 00:40 - Visualizing directional bonds
else if (elapsedTime >= 40 && elapsedTime < 60) {
  const tetrahedralCarbon = scene.getObjectByName("tetrahedralCarbon");
  if (tetrahedralCarbon) {
    // Trigger bond-specific animations (e.g., vibrating rods) via update function if provided
    if (typeof window.updateTetrahedralCarbon === "function") {
      window.updateTetrahedralCarbon(elapsedTime, { animateBonds: true });
    } else {
      tetrahedralCarbon.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Maintain camera steady on the model
  if (camera) {
    camera.lookAt(0, 0, 0);
  }
  
  // (Caption: "Visualizing directional bonds")
}

// ──────────────────────────────────────────────
// At 01:00 - Bond angles: dynamic overlay and camera rotation
else if (elapsedTime >= 60 && elapsedTime < 80) {
  // Ensure Bond_Angle_and_Lattice_Overlay is visible and fade it in over 2 seconds
  const bondOverlay = scene.getObjectByName("Bond_Angle_and_Lattice_Overlay");
  if (bondOverlay) {
    bondOverlay.visible = true;
    const overlayFadeProgress = (elapsedTime - 60 < 2) ? (elapsedTime - 60) / 2 : 1;
    applyFade(bondOverlay, overlayFadeProgress);
    if (typeof window.updateBondOverlay === "function") {
      window.updateBondOverlay(elapsedTime);
    }
  }
  
  // Continue updating tetrahedralCarbon normally
  const tetrahedralCarbon = scene.getObjectByName("tetrahedralCarbon");
  if (tetrahedralCarbon) {
    if (typeof window.updateTetrahedralCarbon === "function") {
      window.updateTetrahedralCarbon(elapsedTime);
    } else {
      tetrahedralCarbon.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Animate camera rotation around the molecule
  if (camera) {
    // Compute an angle so that the camera starts orbiting from 01:00 onwards
    const angle = (elapsedTime - 60) * (Math.PI / 20); // Slow orbit: completes ~2π in ~40 sec if continued
    const radius = 8;
    const targetX = radius * Math.sin(angle);
    const targetZ = radius * Math.cos(angle);
    camera.position.lerp(new THREE.Vector3(targetX, 0, targetZ), 0.1);
    camera.lookAt(0, 0, 0);
  }
  
  // (Caption: "Bond angles: 109.5° each")
}

// ──────────────────────────────────────────────
// At 01:20 - Electron cloud visualization
else if (elapsedTime >= 80 && elapsedTime < 100) {
  // Fade in the electron cloud over 2 seconds
  const electronCloud = scene.getObjectByName("electronCloud");
  if (electronCloud) {
    electronCloud.visible = true;
    const cloudFadeProgress = (elapsedTime - 80 < 2) ? (elapsedTime - 80) / 2 : 1;
    applyFade(electronCloud, cloudFadeProgress);
    if (typeof window.updateElectronCloud === "function") {
      window.updateElectronCloud(elapsedTime);
    }
  }
  
  // Keep bond overlay visible at full opacity
  const bondOverlay = scene.getObjectByName("Bond_Angle_and_Lattice_Overlay");
  if (bondOverlay) {
    bondOverlay.visible = true;
  }
  
  // Continue updating tetrahedralCarbon
  const tetrahedralCarbon = scene.getObjectByName("tetrahedralCarbon");
  if (tetrahedralCarbon) {
    if (typeof window.updateTetrahedralCarbon === "function") {
      window.updateTetrahedralCarbon(elapsedTime);
    } else {
      tetrahedralCarbon.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Continue smooth camera orbit
  if (camera) {
    const angle = (elapsedTime - 60) * (Math.PI / 20);
    const radius = 8;
    const targetX = radius * Math.sin(angle);
    const targetZ = radius * Math.cos(angle);
    camera.position.lerp(new THREE.Vector3(targetX, 0, targetZ), 0.1);
    camera.lookAt(0, 0, 0);
  }
  
  // (Caption: "Electron cloud visualization")
}

// ──────────────────────────────────────────────
// At 01:40 - Recap and molecular significance
else if (elapsedTime >= 100 && elapsedTime < 120) {
  // Continue the full 360° rotation recap using updated parameters
  const tetrahedralCarbon = scene.getObjectByName("tetrahedralCarbon");
  if (tetrahedralCarbon) {
    if (typeof window.updateTetrahedralCarbon === "function") {
      window.updateTetrahedralCarbon(elapsedTime, { recap: true });
    } else {
      tetrahedralCarbon.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Briefly highlight bond overlay and electron cloud by keeping full opacity
  const bondOverlay = scene.getObjectByName("Bond_Angle_and_Lattice_Overlay");
  if (bondOverlay) {
    bondOverlay.traverse(child => {
      if (child.isMesh && child.material) { child.material.opacity = 1; }
    });
  }
  const electronCloud = scene.getObjectByName("electronCloud");
  if (electronCloud) {
    electronCloud.traverse(child => {
      if (child.isMesh && child.material) { child.material.opacity = 1; }
    });
  }
  
  // Accelerate camera to complete a full 360° orbit during this phase
  if (camera) {
    // Calculate a recap angle that completes one full rotation over 20 sec
    const recapAngle = ((elapsedTime - 100) / 20) * (2 * Math.PI);
    const radius = 8;
    const targetX = radius * Math.sin(recapAngle);
    const targetZ = radius * Math.cos(recapAngle);
    camera.position.lerp(new THREE.Vector3(targetX, 0, targetZ), 0.1);
    camera.lookAt(0, 0, 0);
  }
  
  // (Caption: "Recap and molecular significance")
}
    
    // Update controls
    controls.update();
    
    // Render the scene
    renderer.render(scene, camera);
}

// Update captions and UI based on current time
function updateUI() {
    // Add the time offset to the elapsed time
    const elapsedTime = clock.getElapsedTime() + (window.timeOffset || 0);
    
    // Constrain elapsed time to be between 0 and animation duration
    const constrainedTime = Math.max(0, Math.min(elapsedTime, ANIMATION_DURATION));
    
    // Update progress bar
    const progressBar = document.getElementById('progress-bar');
    const progress = constrainedTime / ANIMATION_DURATION;
    progressBar.style.width = progress * 100 + '%';
    
    // Update captions
    const captions = document.querySelectorAll('.caption');
    captions.forEach(caption => {
        const timeStr = caption.getAttribute('data-time');
        const [min, sec] = timeStr.split(':').map(Number);
        const timeInSeconds = min * 60 + sec;
        
        // Show caption if we're within 5 seconds of its timecode
        if (constrainedTime >= timeInSeconds && constrainedTime < timeInSeconds + 5) {
            caption.style.display = 'block';
        } else {
            caption.style.display = 'none';
        }
    });
    
    return constrainedTime; // Return the constrained time for use in animations
}

// Window resize handler
function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

// Set up UI controls
function setupControls() {
    // Store button references
    var playPauseButton = document.getElementById('play-pause');
    var resetButton = document.getElementById('reset');
    var rewindButton = document.getElementById('rewind');
    var fastForwardButton = document.getElementById('fast-forward');
    
    // Only set up event listeners if elements exist
    if (!playPauseButton || !resetButton || !rewindButton || !fastForwardButton) {
        console.warn('Some control buttons not found in the DOM. Control setup may be incomplete.');
    }
    
    let playbackSpeed = 1.0; // Normal speed
    
    // Initialize time offset if not already defined
    window.timeOffset = window.timeOffset || 0;
    
    // Set up play/pause button if it exists
    if (playPauseButton) {
        playPauseButton.addEventListener('click', () => {
            isPlaying = !isPlaying;
            playPauseButton.textContent = isPlaying ? 'Pause' : 'Play';
            if (isPlaying) {
                if (clock) clock.start();
            } else {
                if (clock) clock.stop();
            }
        });
    }
    
    // Set up reset button if it exists
    if (resetButton) {
        resetButton.addEventListener('click', () => {
            if (typeof THREE !== 'undefined') {
                clock = new THREE.Clock();
            }
            window.timeOffset = 0; // Reset the time offset
            isPlaying = true;
            if (playPauseButton) playPauseButton.textContent = 'Pause';
            playbackSpeed = 1.0; // Reset speed to normal
        });
    }
    
    // Set up rewind button if it exists
    if (rewindButton) {
        rewindButton.addEventListener('click', () => {
            // Decrease the time offset by 10 seconds (but don't go below negative animation duration)
            window.timeOffset = Math.max(window.timeOffset - 10, -ANIMATION_DURATION);
            
            // Ensure playing state
            isPlaying = true;
            if (playPauseButton) playPauseButton.textContent = 'Pause';
        });
    }
    
    // Set up fast-forward button if it exists
    if (fastForwardButton) {
        fastForwardButton.addEventListener('click', () => {
            // Increase the time offset by 10 seconds (but don't exceed animation duration)
            window.timeOffset = Math.min(window.timeOffset + 10, ANIMATION_DURATION);
            
            // Ensure playing state
            isPlaying = true;
            if (playPauseButton) playPauseButton.textContent = 'Pause';
        });
    }
}

// Initialize and start animation
init();
animate();
