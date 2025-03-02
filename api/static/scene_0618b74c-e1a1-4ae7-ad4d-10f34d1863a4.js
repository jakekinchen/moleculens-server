// Scientific visualization: Visualizing Acetic Acid: Structure and Function
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
// Create a group for the lab environment
const labEnvironment = new THREE.Group();
window.labEnvironment = labEnvironment;

// Dimensions for the lab table surface
const tableLength = 100;
const tableWidth = 50;

// Create the lab table plane (the base surface)
// The plane is created in the XZ plane, so we rotate it to lie horizontally.
const planeGeometry = new THREE.PlaneGeometry(tableLength, tableWidth);
const planeMaterial = new THREE.MeshPhongMaterial({ 
  color: 0x808080,       // neutral gray
  shininess: 0,          // matte finish
  emissive: 0x050505     // subtle emissive effect for soft glow
});
const tablePlane = new THREE.Mesh(planeGeometry, planeMaterial);
tablePlane.rotation.x = -Math.PI / 2;
labEnvironment.add(tablePlane);

// Create a grid helper to overlay the clean grid pattern.
// The grid helper is square by default; so we scale it along one axis to match the table width ratio.
const gridSize = tableLength;
const gridDivisions = 10;
const gridHelper = new THREE.GridHelper(gridSize, gridDivisions, 0xffffff, 0xffffff);
// Fade the grid lines slightly
gridHelper.material.opacity = 0.5;
gridHelper.material.transparent = true;
// Scale the grid helper to cover the table dimensions (tableWidth/tableLength on the Z axis)
gridHelper.scale.set(1, 1, tableWidth / tableLength);
// Offset slightly above the table plane to avoid z-fighting
gridHelper.position.y = 0.01;
labEnvironment.add(gridHelper);

// Add the lab environment group to the scene
scene.add(labEnvironment);

// GeometryAgent LLM-generated code
(function() {
    // Materials
    const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x808080, shininess: 100 });
    const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 100 });
    const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 100 });
    const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, emissive: 0x333333, shininess: 100 });
    const electronCloudMaterial = new THREE.MeshPhongMaterial({ color: 0x00ffff, transparent: true, opacity: 0.4 });
    const overlayMaterial = new THREE.MeshPhongMaterial({ color: 0x00ffff, transparent: true, opacity: 0.5, emissive: 0x00ffff, shininess: 100 });
    const waterMaterial = new THREE.MeshPhongMaterial({ color: 0x66ccff, transparent: true, opacity: 0.6 });

    // Atom common geometry
    const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

    // Utility function to create a cylinder (bond) between two points
    function createBond(start, end, material) {
        const direction = new THREE.Vector3().subVectors(end, start);
        const length = direction.length();
        const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
        const bond = new THREE.Mesh(bondGeometry, material);
        // Position: midway between start and end
        bond.position.copy(start).addScaledVector(direction, 0.5);
        // Align cylinder with the direction vector
        bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.normalize());
        return bond;
    }

    // Function to create a single acetic acid molecule (CH3COOH)
    // Returns a THREE.Group with properties referencing key atoms for later interactions.
    function createAceticAcid() {
        const molecule = new THREE.Group();

        // --- Create atoms ---
        // Carbon atoms
        const c1 = new THREE.Mesh(atomGeometry, carbonMaterial);
        c1.position.set(0, 0, 0);
        c1.scale.set(0.5, 0.5, 0.5);

        const c2 = new THREE.Mesh(atomGeometry, carbonMaterial);
        c2.position.set(1.5, 0, 0);
        c2.scale.set(0.5, 0.5, 0.5);

        // Oxygen atoms
        // O1: carbonyl oxygen (C=O)
        const o1 = new THREE.Mesh(atomGeometry, oxygenMaterial);
        o1.position.set(2.7, 0, 0);
        o1.scale.set(0.55, 0.55, 0.55);

        // O2: hydroxyl oxygen (O–H)
        const o2 = new THREE.Mesh(atomGeometry, oxygenMaterial);
        o2.position.set(1.5, 1.2, 0);
        o2.scale.set(0.55, 0.55, 0.55);

        // Hydrogen atoms (attached to CH3 and hydroxyl group)
        const h1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
        h1.position.set(-0.8, 0.8, 0);
        h1.scale.set(0.3, 0.3, 0.3);

        const h2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
        h2.position.set(-0.8, -0.8, 0);
        h2.scale.set(0.3, 0.3, 0.3);

        const h3 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
        h3.position.set(0, 0, -0.8);
        h3.scale.set(0.3, 0.3, 0.3);

        // Hydrogen on hydroxyl group
        const hOH = new THREE.Mesh(atomGeometry, hydrogenMaterial);
        hOH.position.set(1.5, 1.8, 0);
        hOH.scale.set(0.3, 0.3, 0.3);

        // --- Create bonds ---
        const bonds = [];
        bonds.push(createBond(c1.position, c2.position, bondMaterial));      // C1-C2
        bonds.push(createBond(c2.position, o1.position, bondMaterial));      // C2=O1
        bonds.push(createBond(c2.position, o2.position, bondMaterial));      // C2-O2
        bonds.push(createBond(c1.position, h1.position, bondMaterial));      // C1-H1
        bonds.push(createBond(c1.position, h2.position, bondMaterial));      // C1-H2
        bonds.push(createBond(c1.position, h3.position, bondMaterial));      // C1-H3
        bonds.push(createBond(o2.position, hOH.position, bondMaterial));     // O2-HOH

        // --- Add electron cloud density effects ---
        // For carbonyl group (around O1)
        const electronCloudGeo = new THREE.SphereGeometry(0.8, 16, 16);
        const cloud1 = new THREE.Mesh(electronCloudGeo, electronCloudMaterial);
        cloud1.position.copy(o1.position);
        // For hydroxyl group (around O2)
        const cloud2 = new THREE.Mesh(electronCloudGeo, electronCloudMaterial);
        cloud2.position.copy(o2.position);

        // Add atoms and bonds to molecule group
        molecule.add(c1, c2, o1, o2, h1, h2, h3, hOH, cloud1, cloud2, ...bonds);

        // Store key atom references for later (e.g., hydrogen bonding in dimer formation)
        molecule.userData = {
            // For hydrogen bonding: typically the hydroxyl hydrogen and carbonyl oxygen
            hOH: hOH,
            o1: o1,
            o2: o2
        };

        return molecule;
    }

    // --- Create the simulation group that consolidates molecules and solvent effects ---
    const aceticAcidSimulation = new THREE.Group();
    window.aceticAcidSimulation = aceticAcidSimulation; // expose globally

    // Create two acetic acid molecules for dimer formation
    const molecule1 = createAceticAcid();
    // First molecule remains at origin (relative to its group)
    aceticAcidSimulation.add(molecule1);

    const molecule2 = createAceticAcid();
    // Offset the second molecule to simulate approach for dimer formation.
    molecule2.position.set(0.5, -2.5, 0);
    aceticAcidSimulation.add(molecule2);

    // --- Create dimer hydrogen bond overlay(s) ---
    // Using hydrogen-bonding between molecule1's hydroxyl hydrogen and molecule2's carbonyl oxygen, and vice versa.
    const overlayBonds = [];
    // Compute world positions for accurate bond placement:
    const pos1 = new THREE.Vector3();
    molecule1.userData.hOH.getWorldPosition(pos1);
    const pos2 = new THREE.Vector3();
    molecule2.userData.o1.getWorldPosition(pos2);
    overlayBonds.push(createBond(pos1, pos2, overlayMaterial));

    // Second hydrogen bond: from molecule2's hydroxyl hydrogen to molecule1's carbonyl oxygen.
    const pos3 = new THREE.Vector3();
    molecule2.userData.hOH.getWorldPosition(pos3);
    const pos4 = new THREE.Vector3();
    molecule1.userData.o1.getWorldPosition(pos4);
    overlayBonds.push(createBond(pos3, pos4, overlayMaterial));

    // Add hydrogen bond overlays to simulation group
    overlayBonds.forEach(bond => aceticAcidSimulation.add(bond));

    // --- Create solvent (water molecules) simulation ---
    // Generate a few water molecules around the acetic acid molecules.
    const waterGroup = new THREE.Group();
    const waterGeometry = new THREE.SphereGeometry(0.3, 16, 16);
    for (let i = 0; i < 8; i++) {
        const water = new THREE.Mesh(waterGeometry, waterMaterial);
        // Randomly position water molecules in a region around the simulation group.
        water.position.set(
            -2 + Math.random() * 4,
            -2 + Math.random() * 4,
            -2 + Math.random() * 4
        );
        waterGroup.add(water);
    }
    aceticAcidSimulation.add(waterGroup);

    // --- Add the complete acetic acid simulation group to the scene ---
    scene.add(aceticAcidSimulation);

    // (Any dynamic animations such as gentle rotation, pulsating electron clouds, dimer transitions,
    // and solvent diffusion should be implemented in the animation loop elsewhere.)
})();
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
// Calculate elapsed time (in seconds) and delta time for smooth animations
const elapsedTime = window.animationTime || clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Helper function to fade an object in/out over a specified duration (in seconds)
function fadeObject(object, fadeIn, duration = 2.0) {
  if (!object) return;
  object.traverse(child => {
    if (child.isMesh && child.material) {
      child.material.transparent = true;
      // Determine the progress based on elapsedTime modulated within duration
      // (Assumes fade starts at the beginning of the current timeline section)
      let localTime = (elapsedTime % duration);
      let progress = Math.min(1, localTime / duration);
      if (fadeIn) {
        child.material.opacity = progress;
        child.visible = true;
      } else {
        child.material.opacity = Math.max(0, 1 - progress);
        if (child.material.opacity === 0) {
          child.visible = false;
        }
      }
    }
  });
}

// ──────────────────────────────────────────────────────────────
// TIMECODE 00:00 - Introduction: Lab Environment & Title
if (elapsedTime >= 0 && elapsedTime < 20) {
  // Animate Lab Environment (digital lab table, grid, light effects, title and 2D formula)
  const labEnv = scene.getObjectByName("labEnvironment");
  if (labEnv) {
    labEnv.visible = true;
    if (typeof window.updateLabEnvironment === "function") {
      window.updateLabEnvironment(elapsedTime);
    } else {
      // Fallback: a subtle slow rotation for visual effect
      labEnv.rotation.y += 0.1 * deltaTime;
    }
    // Fade in lab environment elements for a smooth appearance
    fadeObject(labEnv, true, 2.0);
  }
  
  // Ensure acetic acid simulation remains hidden at this stage
  const acidSim = scene.getObjectByName("aceticAcidSimulation");
  if (acidSim) {
    acidSim.visible = false;
  }
  
  // Camera: smooth pan over the lab table
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 5, 10);
    camera.position.lerp(targetPosition, 0.02);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional: Overlay caption "Introduction: Acetic Acid Overview" via HTML/CSS.)
}

// ──────────────────────────────────────────────────────────────
// TIMECODE 00:20 - 3D Molecular Structure Revealed
else if (elapsedTime >= 20 && elapsedTime < 40) {
  // Bring in the Acetic Acid 3D simulation (ball-and-stick model, color-coded atoms, rotating formula)
  const acidSim = scene.getObjectByName("aceticAcidSimulation");
  if (acidSim) {
    acidSim.visible = true;
    fadeObject(acidSim, true, 2.0);
    // Use provided update function if available
    if (typeof window.updateAceticAcidSimulation === "function") {
      window.updateAceticAcidSimulation(elapsedTime, { mode: "structure" });
    } else {
      // Fallback: gentle rotation to showcase the molecule
      acidSim.rotation.y += 0.3 * deltaTime;
    }
  }
  
  // Camera: zoom in closer to focus on the 3D model
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 2, 5);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional: Display caption "3D Molecular Structure Revealed".)
}

// ──────────────────────────────────────────────────────────────
// TIMECODE 00:40 - Bonding & Polarity in Action
else if (elapsedTime >= 40 && elapsedTime < 60) {
  // Focus on molecular bonds, glowing polar bonds and electron cloud density effects
  const acidSim = scene.getObjectByName("aceticAcidSimulation");
  if (acidSim) {
    if (typeof window.updateAceticAcidSimulation === "function") {
      // Signal the simulation to highlight bonds and electron density around groups
      window.updateAceticAcidSimulation(elapsedTime, { mode: "bonds" });
    } else {
      // Fallback: slow rotation
      acidSim.rotation.y += 0.2 * deltaTime;
    }
  }
  
  // Camera: adjust to slightly tighter focus on the bond regions
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 2.5, 4);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional: Overlay caption "Bonding & Polarity in Action".)
}

// ──────────────────────────────────────────────────────────────
// TIMECODE 01:00 - Dimer Formation via H Bonds
else if (elapsedTime >= 60 && elapsedTime < 80) {
  // Animate two acetic acid molecules forming a dimer with hydrogen bonds
  const acidSim = scene.getObjectByName("aceticAcidSimulation");
  if (acidSim) {
    if (typeof window.updateAceticAcidSimulation === "function") {
      window.updateAceticAcidSimulation(elapsedTime, { mode: "dimer" });
    } else {
      acidSim.rotation.y += 0.1 * deltaTime;
    }
    // (Assumption: The simulation's internal update function handles molecule duplication,
    //  hydrogen bond overlays and transparency effects.)
  }
  
  // Camera: transition to a view that captures both molecules forming the dimer
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 3, 8);
    camera.position.lerp(targetPosition, 0.04);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional: Display caption "Dimer Formation via H Bonds".)
}

// ──────────────────────────────────────────────────────────────
// TIMECODE 01:20 - Acetic Acid in Vinegar
else if (elapsedTime >= 80 && elapsedTime < 100) {
  // Extend the simulation to include solvent interaction (water molecules, diffusion dynamics)
  const acidSim = scene.getObjectByName("aceticAcidSimulation");
  if (acidSim) {
    if (typeof window.updateAceticAcidSimulation === "function") {
      window.updateAceticAcidSimulation(elapsedTime, { mode: "solvent" });
    } else {
      acidSim.rotation.y += 0.05 * deltaTime;
    }
  }
  
  // Camera: widen the perspective to include the solvent environment
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 5, 12);
    camera.position.lerp(targetPosition, 0.03);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional: Overlay caption "Acetic Acid in Vinegar".)
}

// ──────────────────────────────────────────────────────────────
// TIMECODE 01:40 - Summary: Properties & Reactions Recap
else if (elapsedTime >= 100 && elapsedTime < 120) {
  // Transition to a recap montage summarizing key atomic features, bonding angles, etc.
  const acidSim = scene.getObjectByName("aceticAcidSimulation");
  if (acidSim) {
    if (typeof window.updateAceticAcidSimulation === "function") {
      window.updateAceticAcidSimulation(elapsedTime, { mode: "recap" });
    } else {
      acidSim.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Camera: slowly pull back to reveal the complete scene and recap overlays
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 8, 15);
    camera.position.lerp(targetPosition, 0.02);
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // (Optional: Overlay caption "Summary: Properties & Reactions".)
}

// (Optional: Actions for elapsedTime >= 120 seconds, e.g. loop or stop animation.)
    
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
