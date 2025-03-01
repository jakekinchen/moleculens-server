// Scientific visualization: Unveiling Methanol's 3D Structure
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
// Create a canvas to draw a smooth gradient with an overlaid grid
const canvas = document.createElement('canvas');
canvas.width = 128;
canvas.height = 128;
const context = canvas.getContext('2d');

// Create a vertical gradient for the background (light gray with hints of blue)
const gradient = context.createLinearGradient(0, 0, 0, canvas.height);
gradient.addColorStop(0, "rgba(200,210,230,0.3)"); // top: light, subtle blue-gray
gradient.addColorStop(1, "rgba(200,210,230,0.1)"); // bottom: more transparent
context.fillStyle = gradient;
context.fillRect(0, 0, canvas.width, canvas.height);

// Draw grid lines over the gradient
context.strokeStyle = "rgba(180,200,220,0.5)"; // subtle grid line color
context.lineWidth = 1;
const gridSpacing = 16;
for (let i = 0; i <= canvas.width; i += gridSpacing) {
  // vertical lines
  context.beginPath();
  context.moveTo(i, 0);
  context.lineTo(i, canvas.height);
  context.stroke();

  // horizontal lines
  context.beginPath();
  context.moveTo(0, i);
  context.lineTo(canvas.width, i);
  context.stroke();
}

// Create a texture from the canvas and configure it to repeat
const gridTexture = new THREE.CanvasTexture(canvas);
gridTexture.wrapS = THREE.RepeatWrapping;
gridTexture.wrapT = THREE.RepeatWrapping;
gridTexture.repeat.set(10, 10); // Repeat texture pattern to fill the dimensions

// Create material with the generated texture; using MeshBasicMaterial for ambient appearance
const gridMaterial = new THREE.MeshBasicMaterial({
  map: gridTexture,
  transparent: true,
  opacity: 0.5, // subdued visibility to ensure it's unobtrusive
});

// Create a plane geometry with dimensions 100 x 100 and a very thin depth (1 unit interpreted as thickness)
const gridGeometry = new THREE.PlaneGeometry(100, 100, 1, 1);

// Create the mesh, oriented to lie horizontally in the X-Y plane (facing positive Z)
const gridMesh = new THREE.Mesh(gridGeometry, gridMaterial);
gridMesh.position.z = -0.5; // Slight offset in depth to position it as a background

// Store a global reference to the grid mesh and add it to the scene
window.gridBackground = gridMesh;
scene.add(gridMesh);

// GeometryAgent LLM-generated code
// Create atom materials
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x000000 });
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, emissive: 0xff5555, emissiveIntensity: 0.5 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff });

// Create bond materials
// Metallic finish for bonds from carbon
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x888888, shininess: 100 });
// For the oxygen-hydrogen bond, use a thicker material and note that a gradient effect is desired.
const ohBondMaterial = new THREE.MeshPhongMaterial({ color: 0xaaaaaa, shininess: 100 });
ohBondMaterial.userData.gradient = true; // flag for gradient effect

// Create a shared sphere geometry for atoms (radius 1, will be scaled later)
const sphereGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create a group for the methanol molecule
const methanol = new THREE.Group();
window.methanolMolecule = methanol;

// Define positions (units chosen for an approx 10-unit diameter overall structure)
// Central carbon at origin
const carbonPosition = new THREE.Vector3(0, 0, 0);
// Oxygen positioned along positive x-axis from carbon
const oxygenPosition = new THREE.Vector3(2.2, 0, 0);
// Three hydrogens attached to carbon (CH3)
// Positioned roughly in a tetrahedral arrangement relative to carbon
const hC1Position = new THREE.Vector3(-0.8, 0.8, 0);
const hC2Position = new THREE.Vector3(-0.8, -0.8, 0);
const hC3Position = new THREE.Vector3(0, 0, -1.2);
// One hydrogen attached to oxygen (O-H)
const hOPosition = new THREE.Vector3(3.0, 0.6, 0);

// Create atom meshes
// Carbon atom (central hub)
const carbon = new THREE.Mesh(sphereGeometry, carbonMaterial);
carbon.position.copy(carbonPosition);
carbon.scale.set(0.7, 0.7, 0.7);

// Oxygen atom with glow
const oxygen = new THREE.Mesh(sphereGeometry, oxygenMaterial);
oxygen.position.copy(oxygenPosition);
oxygen.scale.set(0.65, 0.65, 0.65);

// Hydrogen atoms attached to carbon with pop-up animation flag
const hydrogenC1 = new THREE.Mesh(sphereGeometry, hydrogenMaterial);
hydrogenC1.position.copy(hC1Position);
hydrogenC1.scale.set(0.3, 0.3, 0.3);
hydrogenC1.userData.popUp = true;

const hydrogenC2 = new THREE.Mesh(sphereGeometry, hydrogenMaterial);
hydrogenC2.position.copy(hC2Position);
hydrogenC2.scale.set(0.3, 0.3, 0.3);
hydrogenC2.userData.popUp = true;

const hydrogenC3 = new THREE.Mesh(sphereGeometry, hydrogenMaterial);
hydrogenC3.position.copy(hC3Position);
hydrogenC3.scale.set(0.3, 0.3, 0.3);
hydrogenC3.userData.popUp = true;

// Hydrogen atom attached to oxygen (with pop-up animation flag)
const hydrogenO = new THREE.Mesh(sphereGeometry, hydrogenMaterial);
hydrogenO.position.copy(hOPosition);
hydrogenO.scale.set(0.3, 0.3, 0.3);
hydrogenO.userData.popUp = true;

// Function to create a bond (cylinder) from point 'start' to point 'end'
// Additional parameter bondRadius allows different thickness
function createBond(start, end, bondRadius, material) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  const bondGeometry = new THREE.CylinderGeometry(bondRadius, bondRadius, length, 16);
  const bond = new THREE.Mesh(bondGeometry, material);
  
  // Position bond at midpoint between start and end
  bond.position.copy(start).lerp(end, 0.5);
  
  // Align bond with the direction vector
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize());
  
  return bond;
}

// Create bonds
// Bond from carbon to oxygen (using cylindrical structure with metallic finish)
const bondCtoO = createBond(carbonPosition, oxygenPosition, 0.15, bondMaterial);

// Bonds from carbon to each hydrogen (CH3 bonds)
const bondCtoH1 = createBond(carbonPosition, hC1Position, 0.1, bondMaterial);
const bondCtoH2 = createBond(carbonPosition, hC2Position, 0.1, bondMaterial);
const bondCtoH3 = createBond(carbonPosition, hC3Position, 0.1, bondMaterial);

// Oxygen-hydrogen bond with enhanced thickness and gradient effect, and a gentle vibration animation flag
const bondOtoH = createBond(oxygenPosition, hOPosition, 0.2, ohBondMaterial);
bondOtoH.userData.vibrate = true; // flag to indicate dynamic vibration animation

// Add all atoms and bonds to the methanol group
methanol.add(carbon, oxygen, hydrogenC1, hydrogenC2, hydrogenC3, hydrogenO);
methanol.add(bondCtoO, bondCtoH1, bondCtoH2, bondCtoH3, bondOtoH);

// Set group-level userData for overall animations (slow rotation and periodic zoom on carbon and oxygen)
methanol.userData.rotate = true;
methanol.userData.zoomFocus = [carbon.id, oxygen.id]; // IDs for atoms to be focused on during zoom

// Optionally, you can attach a simple rotation animation to the methanol molecule
// (This function should be called in your main animation loop, passing the delta time)
window.animateMethanol = function(delta) {
  methanol.rotation.y += delta * 0.2; // slow continuous rotation
};

// Finally, add the methanol molecule group to the scene
scene.add(methanol);

// Note: The additional animation effects (pop-up for hydrogens, bond vibration, and zooming on atoms)
// should be implemented in the main render loop or via tween libraries using the userData flags.
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

// Find main objects using scene.getObjectByName()
const gridBackground = scene.getObjectByName("gridBackground");
const methanolMolecule = scene.getObjectByName("methanolMolecule");

//────────────────────────────────────────────
// At 00:00 - Introduction with grid & rotating methanol
if (elapsedTime < 20) {
  // Animate grid background (if it has any update behavior)
  if (gridBackground) {
    // Example: slowly pulse opacity to create a soft-focus effect
    gridBackground.traverse(child => {
      if (child.isMesh && child.material && child.material.opacity !== undefined) {
        child.material.opacity = Math.min(1, 0.2 + Math.abs(Math.sin(elapsedTime)) * 0.8);
        child.material.transparent = true;
      }
    });
  }
  
  // Animate the methanol molecule
  if (methanolMolecule) {
    // Use the provided geometry update function if available
    if (typeof window.animateMethanol === 'function') {
      window.animateMethanol(elapsedTime);
    } else {
      // Fallback: slowly rotate the molecule
      methanolMolecule.rotation.y += 0.2 * deltaTime;
    }
  }
  
  // Ensure camera is positioned for the overall view
  if (camera) {
    const targetPos = new THREE.Vector3(0, 0, 10);
    camera.position.lerp(targetPos, 0.05);
    camera.lookAt(0, 0, 0);
  }
  
  // (Caption: "Methanol molecule illuminated." can be rendered via overlay UI)
}

//────────────────────────────────────────────
// At 00:20 - Zoom in on central carbon & tetrahedral bonds
else if (elapsedTime >= 20 && elapsedTime < 40) {
  if (methanolMolecule) {
    // Use/update molecule via provided function if available
    if (typeof window.animateMethanol === 'function') {
      window.animateMethanol(elapsedTime);
    } else {
      methanolMolecule.rotation.y += 0.15 * deltaTime;
    }
    
    // Find the central carbon atom by name (assuming it is named "carbon")
    const carbon = methanolMolecule.getObjectByName("carbon");
    if (carbon && camera) {
      // Calculate a target position for the camera focusing on carbon (with a slight offset)
      const carbonWorldPos = new THREE.Vector3();
      carbon.getWorldPosition(carbonWorldPos);
      // Offset slightly back for a closer zoom-in view
      const targetCameraPos = carbonWorldPos.clone().add(new THREE.Vector3(0, 0, 3));
      camera.position.lerp(targetCameraPos, 0.05);
      camera.lookAt(carbonWorldPos);
    }
    
    // (Caption: "Central carbon and tetrahedral bonds." could update UI accordingly)
  }
}

//────────────────────────────────────────────
// At 00:40 - Introduce hydrogen atoms with pop-up animation
else if (elapsedTime >= 40 && elapsedTime < 60) {
  if (methanolMolecule) {
    // Continue base molecule rotation/update
    if (typeof window.animateMethanol === 'function') {
      window.animateMethanol(elapsedTime);
    } else {
      methanolMolecule.rotation.y += 0.1 * deltaTime;
    }
    
    // Array of names for the three hydrogen atoms attached directly to carbon
    const hydrogenNames = ["hydrogen1", "hydrogen2", "hydrogen3"];    
    // Each hydrogen pops up sequentially: start times 40, 45, and 50 seconds respectively
    hydrogenNames.forEach((name, index) => {
      const hydrogen = methanolMolecule.getObjectByName(name);
      if (hydrogen) {
        const popUpStart = 40 + index * 5; // 40s, 45s, 50s
        const popUpDuration = 5; // each takes 5 seconds to fully appear
        if (elapsedTime >= popUpStart) {
          // Calculate scale factor: from 0 to original scale (assume original scale 1)
          const progress = Math.min((elapsedTime - popUpStart) / popUpDuration, 1);
          hydrogen.scale.set(progress, progress, progress);
        }
      }
    });
    
    // (Caption: "Hydrogen atoms arranged spatially." can be displayed via UI)
  }
}

//────────────────────────────────────────────
// At 01:00 - Highlight oxygen focus with glowing effect and thickening bond
else if (elapsedTime >= 60 && elapsedTime < 80) {
  if (methanolMolecule) {
    // Update overall molecule (if applicable)
    if (typeof window.animateMethanol === 'function') {
      window.animateMethanol(elapsedTime);
    } else {
      methanolMolecule.rotation.y += 0.1 * deltaTime;
    }
    
    // Focus on oxygen atom (assumed to be named "oxygen")
    const oxygen = methanolMolecule.getObjectByName("oxygen");
    if (oxygen) {
      // Gradually move the oxygen slightly forward (local z-axis)
      const targetZ = 0.5; // adjust offset as needed
      oxygen.position.z = THREE.MathUtils.lerp(oxygen.position.z, targetZ, 0.05);
      // Increase glow: if material supports emissive intensity
      if (oxygen.material) {
        const glow = 0.5 + 0.5 * Math.abs(Math.sin(elapsedTime * 2));
        oxygen.material.emissive = new THREE.Color(0xff0000);
        oxygen.material.emissiveIntensity = glow;
      }
    }
    
    // Adjust the oxygen-carbon bonding cylinder (assumed to be named "oxygenBond")
    const oxygenBond = methanolMolecule.getObjectByName("oxygenBond");
    if (oxygenBond) {
      // Gradually increase the bond's thickness (e.g., by scaling its radius in y-axis)
      oxygenBond.scale.x = THREE.MathUtils.lerp(oxygenBond.scale.x, 1.5, 0.05);
      oxygenBond.scale.y = THREE.MathUtils.lerp(oxygenBond.scale.y, 1.5, 0.05);
      oxygenBond.scale.z = THREE.MathUtils.lerp(oxygenBond.scale.z, 1.5, 0.05);
    }
    
    // (Caption: "Oxygen highlighted showing polarity." can be updated in the UI)
  }
}

//────────────────────────────────────────────
// At 01:20 - Animate the O-H group with vibration and dipole arrows
else if (elapsedTime >= 80 && elapsedTime < 100) {
  if (methanolMolecule) {
    // Continue using base update function
    if (typeof window.animateMethanol === 'function') {
      window.animateMethanol(elapsedTime);
    } else {
      methanolMolecule.rotation.y += 0.1 * deltaTime;
    }
    
    // Animate the bond between oxygen and its hydrogen (assumed "OH_bond")
    const ohBond = methanolMolecule.getObjectByName("OH_bond");
    if (ohBond) {
      // Apply a gentle vibration effect using a sine wave
      ohBond.position.x += 0.005 * Math.sin(elapsedTime * 50);
      ohBond.position.y += 0.005 * Math.cos(elapsedTime * 50);
    }
    
    // For the hydrogen attached to oxygen (assumed "OH_hydrogen"), a slight scale bounce can be used
    const ohHydrogen = methanolMolecule.getObjectByName("OH_hydrogen");
    if (ohHydrogen) {
      const scaleBounce = 1 + 0.1 * Math.sin(elapsedTime * 10);
      ohHydrogen.scale.set(scaleBounce, scaleBounce, scaleBounce);
    }
    
    // Animate electron dipole arrows (assumed named "dipole")
    const dipole = methanolMolecule.getObjectByName("dipole");
    if (dipole) {
      // Oscillate opacity to simulate dynamic charge distribution
      dipole.traverse(child => {
        if (child.isMesh && child.material && child.material.opacity !== undefined) {
          child.material.opacity = 0.5 + 0.5 * Math.abs(Math.sin(elapsedTime * 5));
          child.material.transparent = true;
        }
      });
    }
    
    // (Caption: "Dynamic O-H group and charge." can be shown via overlay)
  }
}

//────────────────────────────────────────────
// At 01:40 - Camera orbits around the methanol molecule with informative labels
else if (elapsedTime >= 100 && elapsedTime < 120) {
  if (methanolMolecule) {
    // Continue the provided molecule animation routine if available
    if (typeof window.animateMethanol === 'function') {
      window.animateMethanol(elapsedTime);
    } else {
      methanolMolecule.rotation.y += 0.05 * deltaTime;
    }
  }
  
  // Orbit the camera in a circular path around the center (assumed at 0,0,0)
  if (camera) {
    const orbitRadius = 10;
    // Progress from 0 to full circle over 20 seconds (from t=100s to t=120s)
    const orbitProgress = (elapsedTime - 100) / 20;
    const angle = orbitProgress * Math.PI * 2;
    camera.position.x = orbitRadius * Math.cos(angle);
    camera.position.z = orbitRadius * Math.sin(angle);
    camera.position.y = 3; // slight elevation for a better 3D perspective
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // Optionally, update informative labels: assume they are children of methanolMolecule 
  // with names like "label_carbon", "label_oxygen", etc.
  const labelNames = ["label_carbon", "label_oxygen", "label_hydrogen1", "label_hydrogen2", "label_hydrogen3", "label_OH"];
  labelNames.forEach(labelName => {
    const label = methanolMolecule.getObjectByName(labelName);
    if (label && label.material) {
      // Fade in labels gradually
      const opacityTarget = 1;
      label.material.opacity = THREE.MathUtils.lerp(label.material.opacity || 0, opacityTarget, 0.05);
      label.material.transparent = true;
    }
  });
  
  // (Caption: "Full 3D view with informative labels." to be rendered via UI overlay)
}

//────────────────────────────────────────────
// At 02:00 and beyond - Final frame and fade out summary graphic
else if (elapsedTime >= 120) {
  // Slowly fade out the methanol molecule and grid background
  if (methanolMolecule) {
    methanolMolecule.traverse(child => {
      if (child.isMesh && child.material && child.material.opacity !== undefined) {
        child.material.opacity = THREE.MathUtils.lerp(child.material.opacity, 0, 0.02);
        child.material.transparent = true;
      }
    });
  }
  if (gridBackground) {
    gridBackground.traverse(child => {
      if (child.isMesh && child.material && child.material.opacity !== undefined) {
        child.material.opacity = THREE.MathUtils.lerp(child.material.opacity, 0, 0.02);
        child.material.transparent = true;
      }
    });
  }
  
  // Optionally, reposition camera slowly to frame the closing graphic
  if (camera) {
    const closingPos = new THREE.Vector3(0, 0, 15);
    camera.position.lerp(closingPos, 0.02);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // (Caption: "Summary: Methanol's key features." and closing graphic 
  // "Methanol: A 3D Journey" should be displayed via an external UI overlay)
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
