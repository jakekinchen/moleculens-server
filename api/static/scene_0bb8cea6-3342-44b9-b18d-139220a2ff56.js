// Scientific visualization: Exploring Propane's 3D Molecular Structure
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
// Create atom materials
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x333333, shininess: 50 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 100 });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x00ffff, emissive: 0x00ffff, emissiveIntensity: 0.5, shininess: 100 });

// Create a common sphere geometry for atoms
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create a group for the propane molecule
const propane = new THREE.Group();
window.propaneMolecule = propane;

// Create Carbon atoms for propane (C3 in a linear chain)
// Using positions along the x-axis
const c1 = new THREE.Mesh(atomGeometry, carbonMaterial);
c1.position.set(-1.5, 0, 0);
c1.scale.set(0.5, 0.5, 0.5);

const c2 = new THREE.Mesh(atomGeometry, carbonMaterial);
c2.position.set(0, 0, 0);
c2.scale.set(0.5, 0.5, 0.5);

const c3 = new THREE.Mesh(atomGeometry, carbonMaterial);
c3.position.set(1.5, 0, 0);
c3.scale.set(0.5, 0.5, 0.5);

// Create Hydrogen atoms
// For left methyl group (CH3 on C1)
const h1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h1.position.set(-1.5, 0.8, 0.8);
h1.scale.set(0.3, 0.3, 0.3);

const h2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h2.position.set(-1.5, -0.8, 0.8);
h2.scale.set(0.3, 0.3, 0.3);

const h3 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h3.position.set(-1.5, 0, -1.0);
h3.scale.set(0.3, 0.3, 0.3);

// For central methylene group (CH2 on C2)
const h4 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h4.position.set(0, 1, 0.5);
h4.scale.set(0.3, 0.3, 0.3);

const h5 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h5.position.set(0, -1, 0.5);
h5.scale.set(0.3, 0.3, 0.3);

// For right methyl group (CH3 on C3)
const h6 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h6.position.set(1.5, 0.8, 0.8);
h6.scale.set(0.3, 0.3, 0.3);

const h7 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h7.position.set(1.5, -0.8, 0.8);
h7.scale.set(0.3, 0.3, 0.3);

const h8 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h8.position.set(1.5, 0, -1.0);
h8.scale.set(0.3, 0.3, 0.3);

// Function to create a bond (a cylinder connecting two points)
function createBond(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);
  
  // Position the bond at the midpoint between start and end
  bond.position.copy(start).add(direction.multiplyScalar(0.5));
  
  // Orient the cylinder to align with the direction vector
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
  return bond;
}

// Create bonds between atoms
const bonds = [];

// Carbon‑Carbon bonds
bonds.push(createBond(c1.position, c2.position));
bonds.push(createBond(c2.position, c3.position));

// Bonds from C1 to its hydrogens
bonds.push(createBond(c1.position, h1.position));
bonds.push(createBond(c1.position, h2.position));
bonds.push(createBond(c1.position, h3.position));

// Bonds from C2 to its hydrogens
bonds.push(createBond(c2.position, h4.position));
bonds.push(createBond(c2.position, h5.position));

// Bonds from C3 to its hydrogens
bonds.push(createBond(c3.position, h6.position));
bonds.push(createBond(c3.position, h7.position));
bonds.push(createBond(c3.position, h8.position));

// Add all atoms and bonds to the molecule group
propane.add(c1, c2, c3, h1, h2, h3, h4, h5, h6, h7, h8);
bonds.forEach(bond => propane.add(bond));

// Add the propane molecule group to the scene
scene.add(propane);

// Note: Animations for slow rotation and subtle bond expansion/contraction can be handled by the render loop externally.
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
// Get elapsed time and delta time for smooth animations
const elapsedTime = window.animationTime || clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Retrieve the propane molecule object
const molecule = scene.getObjectByName("propaneMolecule");

// Helper function for fading objects in and out
function fadeObject(object, fadeIn, duration = 2.0) {
  if (!object) return;
  object.traverse(child => {
    if (child.isMesh && child.material) {
      child.material.transparent = true;
      // Compute progress within the duration
      let progress = Math.min(1, (elapsedTime - (child.userData.fadeStart || 0)) / duration);
      if (fadeIn) {
        child.material.opacity = progress;
        child.visible = true;
      } else {
        child.material.opacity = 1 - progress;
        if (child.material.opacity <= 0.01) {
          child.visible = false;
        }
      }
    }
  });
}

// ─────────────────────────────────────────────
// At 00:00 - Introducing Propane Structure
if (elapsedTime >= 0 && elapsedTime < 20) {
  if (molecule) {
    // Ensure molecule is visible
    molecule.visible = true;
    
    // Fade in the entire molecule over the first 2 seconds of this phase
    if (!molecule.userData.fadeStart) { 
      molecule.userData.fadeStart = elapsedTime;
    }
    fadeObject(molecule, true, 2.0);
    
    // Call the update function if provided by the geometry (for initial atomic appearance)
    if (typeof window.updatePropaneMolecule === 'function') {
      window.updatePropaneMolecule(elapsedTime, "intro");
    } else {
      // Fallback: slight slow rotation 
      molecule.rotation.y += 0.2 * deltaTime;
    }
  }
  
  // Set initial camera position for introduction
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 30);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(molecule ? molecule.position : new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional) Display caption: "Introducing Propane Structure"
}

// ─────────────────────────────────────────────
// At 00:20 - Carbon Backbone Highlight
else if (elapsedTime >= 20 && elapsedTime < 40) {
  if (molecule) {
    // Call update function for highlighting central carbon-carbon bonds
    if (typeof window.updatePropaneMolecule === 'function') {
      window.updatePropaneMolecule(elapsedTime, "highlight");
    } else {
      // Fallback: slowly zoom in and accentuate bonds (if bonds are a child object)
      molecule.rotation.y += 0.1 * deltaTime;
    }
    
    // If there is a bonds child inside the molecule, fade it in for a glowing effect
    const bonds = molecule.getObjectByName("bonds");
    if (bonds) {
      if (!bonds.userData.fadeStart) {
        bonds.userData.fadeStart = elapsedTime;
      }
      fadeObject(bonds, true, 2.0);
      bonds.visible = true;
    }
  }
  
  // Move camera to emphasize the central bonds area
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 18);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(molecule ? molecule.position : new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional) Display caption: "Carbon Backbone Highlight"
}

// ─────────────────────────────────────────────
// At 00:40 - 3D Rotation & Orientation
else if (elapsedTime >= 40 && elapsedTime < 60) {
  if (molecule) {
    // Trigger dynamic rotation and enforce color contrasts via update function if available
    if (typeof window.updatePropaneMolecule === 'function') {
      window.updatePropaneMolecule(elapsedTime, "rotate");
    } else {
      // Fallback: rotate slowly to illustrate 3D structure
      molecule.rotation.y += 0.5 * deltaTime;
    }
  }
  
  // Animate camera orbiting around the molecule for a dynamic perspective
  if (camera && molecule) {
    let angle = (elapsedTime - 40) * 0.03;
    const radius = 18;
    camera.position.x = molecule.position.x + radius * Math.sin(angle);
    camera.position.z = molecule.position.z + radius * Math.cos(angle);
    camera.position.y = molecule.position.y + 5; // slight elevation for depth
    camera.lookAt(molecule.position);
  }
  
  // (Optional) Display caption: "3D Rotation & Orientation"
}

// ─────────────────────────────────────────────
// At 01:00 - Hydrogen Bonding Angles
else if (elapsedTime >= 60 && elapsedTime < 80) {
  if (molecule) {
    // Focus on hydrogen atoms and animate bonding angles
    if (typeof window.updatePropaneMolecule === 'function') {
      window.updatePropaneMolecule(elapsedTime, "bondAngles");
    } else {
      // Fallback: apply a subtle pulsating scale to hydrogen atom groups
      // Assuming hydrogen atoms are grouped under children named "hydrogen"
      const hydrogens = molecule.getObjectByName("hydrogen");
      if (hydrogens) {
        let scaleFactor = 1 + 0.05 * Math.sin((elapsedTime - 60) * Math.PI);
        hydrogens.scale.set(scaleFactor, scaleFactor, scaleFactor);
      }
    }
  }
  
  // Keep the camera focused; a slight dolly in might emphasize bonding angles
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 16);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(molecule ? molecule.position : new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional) Display caption: "Hydrogen Bonding Angles"
}

// ─────────────────────────────────────────────
// At 01:20 - Spatial Geometry & Stability
else if (elapsedTime >= 80 && elapsedTime < 100) {
  if (molecule) {
    // Reassemble full molecule view and reveal atomic labels
    if (typeof window.updatePropaneMolecule === 'function') {
      window.updatePropaneMolecule(elapsedTime, "labels");
    } else {
      // Fallback: gently reset any local transformations applied earlier
      molecule.rotation.y += 0.2 * deltaTime;
    }
    
    // Fade in any label objects if they exist as children
    const labels = molecule.getObjectByName("labels");
    if (labels) {
      labels.visible = true;
      if (!labels.userData.fadeStart) {
        labels.userData.fadeStart = elapsedTime;
      }
      fadeObject(labels, true, 2.0);
    }
  }
  
  // Camera sweeps around the molecule to showcase symmetry
  if (camera && molecule) {
    let sweepAngle = (elapsedTime - 80) * 0.02;
    const distance = 20;
    camera.position.x = molecule.position.x + distance * Math.sin(sweepAngle);
    camera.position.z = molecule.position.z + distance * Math.cos(sweepAngle);
    camera.position.y = molecule.position.y + 5;
    camera.lookAt(molecule.position);
  }
  
  // (Optional) Display caption: "Spatial Geometry & Stability"
}

// ─────────────────────────────────────────────
// At 01:40 - Final Overview of Propane
else if (elapsedTime >= 100 && elapsedTime < 120) {
  if (molecule) {
    // Final slow rotation and glowing labels overview
    if (typeof window.updatePropaneMolecule === 'function') {
      window.updatePropaneMolecule(elapsedTime, "final");
    } else {
      molecule.rotation.y += 0.1 * deltaTime;
    }
    
    // Ensure labels remain softly glowing
    const labels = molecule.getObjectByName("labels");
    if (labels) {
      labels.visible = true;
      // Optionally add a slow pulsation effect to the label opacity
      let glow = 0.7 + 0.3 * Math.sin((elapsedTime - 100) * 0.05);
      labels.traverse(child => {
        if (child.isMesh && child.material) {
          child.material.opacity = glow;
          child.material.transparent = true;
        }
      });
    }
  }
  
  // Update camera for a final, captivating view with a gradient background focus
  if (camera && molecule) {
    let finalAngle = (elapsedTime - 100) * 0.015;
    const finalRadius = 22;
    camera.position.x = molecule.position.x + finalRadius * Math.sin(finalAngle);
    camera.position.z = molecule.position.z + finalRadius * Math.cos(finalAngle);
    camera.position.y = molecule.position.y + 8;
    camera.lookAt(molecule.position);
  }
  
  // Optional: Transition scene background to a gradient (if supported)
  if (scene.background && typeof window.updateSceneBackground === 'function') {
    window.updateSceneBackground(elapsedTime);
  }
  
  // (Optional) Display caption: "Final Overview of Propane"
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
    
    // Constrain elapsed time to be between 0 and TOTAL_DURATION
    const constrainedTime = Math.max(0, Math.min(elapsedTime, TOTAL_DURATION));
    
    // Update progress bar
    const progressBar = document.getElementById('progress-bar');
    const progress = constrainedTime / TOTAL_DURATION;
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
    const playPauseButton = document.getElementById('play-pause');
    const resetButton = document.getElementById('reset');
    const rewindButton = document.getElementById('rewind');
    const fastForwardButton = document.getElementById('fast-forward');
    
    let playbackSpeed = 1.0; // Normal speed
    
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
        window.timeOffset = 0; // Reset the time offset
        isPlaying = true;
        playPauseButton.textContent = 'Pause';
        playbackSpeed = 1.0; // Reset speed to normal
    });
    
    // Define a global variable to track the time offset
    window.timeOffset = 0;
    
    rewindButton.addEventListener('click', () => {
        // Decrease the time offset by 10 seconds (but don't go below negative total duration)
        window.timeOffset = Math.max(window.timeOffset - 10, -TOTAL_DURATION);
        
        // Ensure playing state
        isPlaying = true;
        playPauseButton.textContent = 'Pause';
    });
    
    fastForwardButton.addEventListener('click', () => {
        // Increase the time offset by 10 seconds (but don't exceed total duration)
        window.timeOffset = Math.min(window.timeOffset + 10, TOTAL_DURATION);
        
        // Ensure playing state
        isPlaying = true;
        playPauseButton.textContent = 'Pause';
    });
}

// Initialize and start animation
init();
animate();
