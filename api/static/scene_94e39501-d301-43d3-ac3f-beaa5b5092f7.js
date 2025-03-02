// Scientific visualization: Understanding Ethene's Double Bond Structure
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
// Create materials
const carbonMaterial = new THREE.MeshPhongMaterial({ 
  color: 0x333333, 
  transparent: true, 
  opacity: 0.9, 
  shininess: 100 
});
const hydrogenMaterial = new THREE.MeshPhongMaterial({ 
  color: 0xffffff, 
  transparent: true, 
  opacity: 0.9 
});
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999 });
const sigmaMaterial = new THREE.MeshPhongMaterial({ 
  color: 0x00ffff, 
  emissive: 0x00ffff, 
  emissiveIntensity: 0.7 
});
const piMaterial = new THREE.MeshPhongMaterial({ 
  color: 0xffa500, 
  transparent: true, 
  opacity: 0.6, 
  emissive: 0xffa500, 
  emissiveIntensity: 0.5 
});

// Create atom geometries
const carbonGeo = new THREE.SphereGeometry(0.3, 32, 32);
const hydrogenGeo = new THREE.SphereGeometry(0.2, 32, 32);

// Create a group for the Ethene molecule
const etheneMolecule = new THREE.Group();
window.etheneMolecule = etheneMolecule;

// Define positions for atoms (all in the XY plane)
// Place the two carbon atoms along the X-axis.
const c1Pos = new THREE.Vector3(-1, 0, 0);
const c2Pos = new THREE.Vector3( 1, 0, 0);

// For the CH bonds, we use a 120° separation relative to the double bond axis.
// Offsets for carbon 1 hydrogens (rotating relative to the direction from c1 to c2):
const angleOffset = Math.PI * 2 / 3; // 120° in radians (~2.094)
const hydrogenDistance = 0.8; 
// For c1, the CH bonds are directed roughly backwards relative to c2
const h1Offset = new THREE.Vector3(
  Math.cos(Math.PI - angleOffset/2) * hydrogenDistance,
  Math.sin(Math.PI - angleOffset/2) * hydrogenDistance,
  0
); // Upper hydrogen for c1
const h2Offset = new THREE.Vector3(
  Math.cos(Math.PI + angleOffset/2) * hydrogenDistance,
  Math.sin(Math.PI + angleOffset/2) * hydrogenDistance,
  0
); // Lower hydrogen for c1

// For carbon 2, the CH bonds point away from c1
const h3Offset = new THREE.Vector3(
  Math.cos(0 + angleOffset/2) * hydrogenDistance,
  Math.sin(0 + angleOffset/2) * hydrogenDistance,
  0
); // Upper hydrogen for c2
const h4Offset = new THREE.Vector3(
  Math.cos(0 - angleOffset/2) * hydrogenDistance,
  Math.sin(0 - angleOffset/2) * hydrogenDistance,
  0
); // Lower hydrogen for c2

// Create carbon atoms
const c1 = new THREE.Mesh(carbonGeo, carbonMaterial);
c1.position.copy(c1Pos);
const c2 = new THREE.Mesh(carbonGeo, carbonMaterial);
c2.position.copy(c2Pos);

// Create hydrogen atoms by cloning their geometry and positioning them relative to their carbons
const h1 = new THREE.Mesh(hydrogenGeo, hydrogenMaterial);
h1.position.copy(c1Pos.clone().add(h1Offset));
const h2 = new THREE.Mesh(hydrogenGeo, hydrogenMaterial);
h2.position.copy(c1Pos.clone().add(h2Offset));
const h3 = new THREE.Mesh(hydrogenGeo, hydrogenMaterial);
h3.position.copy(c2Pos.clone().add(h3Offset));
const h4 = new THREE.Mesh(hydrogenGeo, hydrogenMaterial);
h4.position.copy(c2Pos.clone().add(h4Offset));

// Utility function to create a bond (cylinder) between two points
function createBond(start, end, radius, material) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  // CylinderGeometry is created along the Y-axis by default.
  const bondGeometry = new THREE.CylinderGeometry(radius, radius, length, 16);
  const bond = new THREE.Mesh(bondGeometry, material);
  
  // Position: move to midpoint between start and end.
  const midpoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  bond.position.copy(midpoint);
  
  // Align the bond with the direction vector.
  bond.quaternion.setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    direction.clone().normalize()
  );
  
  return bond;
}

// Create bonds between carbons and hydrogens (use a smaller radius for CH bonds)
const chBondRadius = 0.05;
const bondC1H1 = createBond(c1.position, h1.position, chBondRadius, bondMaterial);
const bondC1H2 = createBond(c1.position, h2.position, chBondRadius, bondMaterial);
const bondC2H3 = createBond(c2.position, h3.position, chBondRadius, bondMaterial);
const bondC2H4 = createBond(c2.position, h4.position, chBondRadius, bondMaterial);

// Create the double bond between carbons with its two components

// 1. Sigma bond: a cylinder along the internuclear axis (use a thicker radius)
const sigmaBondRadius = 0.1;
const sigmaBond = createBond(c1.position, c2.position, sigmaBondRadius, sigmaMaterial);

// 2. Pi bonds: two diffuse electron clouds above and below the molecular plane.
// We'll represent these as semi-transparent circular disks. The midpoint of the double bond:
const midPoint = new THREE.Vector3().addVectors(c1.position, c2.position).multiplyScalar(0.5);
// Create a circle geometry; note that the circle lies in the XY plane by default.
const circleGeo = new THREE.CircleGeometry(0.25, 32);

// Upper pi cloud: offset upward along Z
const piBondUpper = new THREE.Mesh(circleGeo, piMaterial);
piBondUpper.position.copy(midPoint);
piBondUpper.position.z += 0.2; // adjust offset above the plane

// Lower pi cloud: offset downward along Z
const piBondLower = new THREE.Mesh(circleGeo, piMaterial);
piBondLower.position.copy(midPoint);
piBondLower.position.z -= 0.2; // adjust offset below the plane

// Group all molecular parts together
etheneMolecule.add(c1, c2, h1, h2, h3, h4);
etheneMolecule.add(bondC1H1, bondC1H2, bondC2H3, bondC2H4);
etheneMolecule.add(sigmaBond, piBondUpper, piBondLower);

// Add the complete Ethene molecule to the scene
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
    
    // Get the adjusted time with offset applied
    const adjustedTime = updateUI();
    
    // We'll make this time variable available to the animation code
    window.animationTime = adjustedTime;
    
    // Animation created by the AnimationAgent
// Get current elapsed time and delta time for smooth animations
const elapsedTime = window.animationTime || clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Retrieve the ethene molecule model from the scene (ensure it exists)
const etheneMolecule = scene.getObjectByName("etheneMolecule");

// Helper function to fade an object in or out over a given duration.
// It assumes that object's material supports opacity and transparency.
function fadeObject(object, fadeIn, startTime, duration) {
  if (!object) return;
  // Calculate local progress for this fade segment
  const localTime = elapsedTime - startTime;
  const progress = Math.min(Math.max(localTime / duration, 0), 1);
  // Traverse all meshes in the object
  object.traverse(child => {
    if (child.isMesh && child.material) {
      child.material.transparent = true;
      child.material.opacity = fadeIn ? progress : (1 - progress);
      // Hide object if fully faded out
      child.visible = fadeIn ? (child.material.opacity > 0) : (child.material.opacity > 0);
    }
  });
}

// ─────────────────────────────────────────────────────────────
// [00:00 - 00:20] Overview of Ethene Molecule
if (elapsedTime >= 0 && elapsedTime < 20) {
  if (etheneMolecule) {
    // Ensure the molecule is visible and fade it in over 2 seconds at scene start
    fadeObject(etheneMolecule, true, 0, 2.0);

    // Use the provided update function if available, otherwise a simple rotation
    if (typeof window.updateEtheneMolecule === 'function') {
      // Inform the update function that we are in 'overview' mode
      window.updateEtheneMolecule(elapsedTime, { mode: "overview" });
    } else {
      // Fallback: slowly rotate to give a 3D view
      etheneMolecule.rotation.y += 0.2 * deltaTime;
    }
  }

  // Position the camera to show a full view of the molecule
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 15);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(0, 0, 0);
  }
  // (Caption displayed externally: "Ethene's planar structure: 2 carbons and 4 hydrogens arranged symmetrically.")
}

// ─────────────────────────────────────────────────────────────
// [00:20 - 00:40] Zoom in and Highlight Carbon-Carbon Double Bond
else if (elapsedTime >= 20 && elapsedTime < 40) {
  if (etheneMolecule) {
    // Optionally use an update function to highlight the double bond with a glowing effect
    if (typeof window.updateEtheneMolecule === 'function') {
      window.updateEtheneMolecule(elapsedTime, { mode: "highlightDoubleBond" });
    } else {
      // Fallback: Increase emissive intensity for a glowing effect on the molecule (if material supports it)
      etheneMolecule.traverse(child => {
        if (child.isMesh && child.material && child.material.emissive) {
          child.material.emissiveIntensity = 1.5;
        }
      });
    }
  }

  // Move camera to zoom in on the carbon-carbon bond area (assumed near the center)
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 7);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(0, 0, 0);
  }
  // (Caption displayed externally: "Carbon-carbon bond highlighted, showing the elegant double bond connection.")
}

// ─────────────────────────────────────────────────────────────
// [00:40 - 01:00] Split View of the Double Bond (Sigma and Pi Components)
else if (elapsedTime >= 40 && elapsedTime < 60) {
  if (etheneMolecule) {
    // Trigger a split view mode that separates the visual components of the double bond.
    if (typeof window.updateEtheneMolecule === 'function') {
      window.updateEtheneMolecule(elapsedTime, { mode: "splitView" });
    } else {
      // Fallback: tilt the molecule slightly to hint at separation of components
      etheneMolecule.rotation.y += 0.1 * deltaTime;
      etheneMolecule.rotation.x = Math.sin(elapsedTime * 0.5) * 0.1;
    }
  }

  // Gradually adjust the camera to a neutral central view that accommodates a split view
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 7);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(0, 0, 0);
  }
  // (Caption displayed externally: "Double bond components: Sigma and pi bonds arise from different orbital overlaps.")
}

// ─────────────────────────────────────────────────────────────
// [01:00 - 01:20] Close-Up on the Sigma Bond
else if (elapsedTime >= 60 && elapsedTime < 80) {
  if (etheneMolecule) {
    // Focus on the sigma bond component: head-on overlap represented as a cylinder between carbons.
    if (typeof window.updateSigmaBond === 'function') {
      window.updateSigmaBond(elapsedTime);
    } else if (typeof window.updateEtheneMolecule === 'function') {
      // Alternatively, pass a mode parameter if a dedicated sigma update is not available
      window.updateEtheneMolecule(elapsedTime, { mode: "sigmaBond" });
    }
  }
  
  // Adjust the camera to a close-up view on the sigma bond region
  if (camera) {
    // Assume sigma bond lies along the x-axis; move camera slightly from the side
    const targetPosition = new THREE.Vector3(2, 0, 5);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(0, 0, 0);
  }
  // (Caption displayed externally: "Sigma bond: Head-on overlap creates a concentrated electron density along the internuclear axis.")
}

// ─────────────────────────────────────────────────────────────
// [01:20 - 01:40] Focus on the Pi Bond
else if (elapsedTime >= 80 && elapsedTime < 100) {
  if (etheneMolecule) {
    // Shift the focus to the pi bond: side-to-side overlap with electron clouds above and below the molecular plane.
    if (typeof window.updatePiBond === 'function') {
      window.updatePiBond(elapsedTime);
    } else if (typeof window.updateEtheneMolecule === 'function') {
      window.updateEtheneMolecule(elapsedTime, { mode: "piBond" });
    }
  }
  
  // Adjust the camera to reveal electron clouds above and below the molecule.
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 3, 6);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  // (Caption displayed externally: "Pi bond: Side-by-side orbital overlap forms electron clouds above and below the molecular plane.")
}

// ─────────────────────────────────────────────────────────────
// [01:40 - 02:00] Final Integrated View of Ethene Molecule
else if (elapsedTime >= 100 && elapsedTime < 120) {
  if (etheneMolecule) {
    // Bring both bonding types back together into the integrated ethene structure.
    if (typeof window.updateEtheneMolecule === 'function') {
      window.updateEtheneMolecule(elapsedTime, { mode: "integratedView" });
    } else {
      // Fallback: slowly reset rotations and visual effects
      etheneMolecule.rotation.x += 0.05 * deltaTime;
      etheneMolecule.rotation.y += 0.05 * deltaTime;
    }
    
    // Ensure the glowing effect on the double bond is normalized
    etheneMolecule.traverse(child => {
      if (child.isMesh && child.material && child.material.emissive) {
        child.material.emissiveIntensity = 1.0;
      }
    });
    
    // Fade the molecule fully in if needed (in case any part was faded out)
    fadeObject(etheneMolecule, true, 100, 2.0);
  }
  
  // Set camera to a balanced overview position of the full molecule
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 10);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  // (Caption displayed externally: "Ethene's stability stems from its combined sigma and pi bonds in a flat molecular plane.")
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
