// Scientific visualization: The Linear Structure of Carbon Dioxide
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
(function(){
    // Create materials
    const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x333333, shininess: 100 });
    const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 100 });
    const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, shininess: 80 });
    // Glow effects can be achieved with emissive properties
    carbonMaterial.emissive = new THREE.Color(0x111111);
    oxygenMaterial.emissive = new THREE.Color(0x220000);
    
    // Create sphere geometry for atoms
    const sphereSegments = 32;
    const carbonGeometry = new THREE.SphereGeometry(0.3, sphereSegments, sphereSegments);
    const oxygenGeometry = new THREE.SphereGeometry(0.35, sphereSegments, sphereSegments);
    
    // Create the main group for CO2 molecule
    const co2Group = new THREE.Group();
    co2Group.name = "CO2_Molecule_Model";
    window.co2Molecule = co2Group;
    
    // Position atoms (linear: O - C - O along x-axis)
    // Central Carbon at (0,0,0)
    const carbon = new THREE.Mesh(carbonGeometry, carbonMaterial);
    carbon.position.set(0, 0, 0);
    co2Group.add(carbon);
    
    // Oxygen atoms at opposite ends (approx. 1 unit away each side for a total length ~2 units)
    const oxygenLeft = new THREE.Mesh(oxygenGeometry, oxygenMaterial);
    oxygenLeft.position.set(-1, 0, 0);
    co2Group.add(oxygenLeft);
    
    const oxygenRight = new THREE.Mesh(oxygenGeometry, oxygenMaterial);
    oxygenRight.position.set(1, 0, 0);
    co2Group.add(oxygenRight);
    
    // Function to create a cylinder (bond) between two points
    function createCylinderBond(start, end, radius) {
        const direction = new THREE.Vector3().subVectors(end, start);
        const length = direction.length();
        const cylGeometry = new THREE.CylinderGeometry(radius, radius, length, 16, 1);
        const bond = new THREE.Mesh(cylGeometry, bondMaterial);
        // Move bond so its center is at the midpoint between atoms.
        bond.position.copy(start).lerp(end, 0.5);
        // Align the bond to the direction vector
        bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
        // Save initial thickness for animation purposes.
        bond.userData = { baseRadius: radius };
        return bond;
    }
    
    // Function to create a double bond (two parallel cylinders with a slight offset)
    function createDoubleBond(start, end) {
        const offsetMagnitude = 0.08; // offset in y direction
        const bond1 = createCylinderBond(
            new THREE.Vector3(start.x, start.y + offsetMagnitude, start.z),
            new THREE.Vector3(end.x, end.y + offsetMagnitude, end.z),
            0.05
        );
        const bond2 = createCylinderBond(
            new THREE.Vector3(start.x, start.y - offsetMagnitude, start.z),
            new THREE.Vector3(end.x, end.y - offsetMagnitude, end.z),
            0.05
        );
        const group = new THREE.Group();
        group.add(bond1);
        group.add(bond2);
        // Store bonds for animation updates:
        group.userData.bonds = [bond1, bond2];
        return group;
    }
    
    // Create double bonds between carbon and oxygen atoms
    const bondLeft = createDoubleBond(carbon.position, oxygenLeft.position);
    co2Group.add(bondLeft);
    
    const bondRight = createDoubleBond(carbon.position, oxygenRight.position);
    co2Group.add(bondRight);
    
    // Utility function to create text sprites for annotations
    function makeTextSprite(message, parameters) {
        if (parameters === undefined) parameters = {};
        const fontface = parameters.hasOwnProperty("fontface") ? parameters["fontface"] : "Arial";
        const fontsize = parameters.hasOwnProperty("fontsize") ? parameters["fontsize"] : 24;
        const borderThickness = parameters.hasOwnProperty("borderThickness") ? parameters["borderThickness"] : 4;
        const borderColor = parameters.hasOwnProperty("borderColor") ? parameters["borderColor"] : { r:0, g:0, b:0, a:1.0 };
        const backgroundColor = parameters.hasOwnProperty("backgroundColor") ? parameters["backgroundColor"] : { r:255, g:255, b:255, a:1.0 };
        // create canvas element
        const canvas = document.createElement("canvas");
        const context = canvas.getContext("2d");
        context.font = fontsize + "px " + fontface;
        // Get size data (height depends only on font size)
        const metrics = context.measureText(message);
        const textWidth = metrics.width;
        // Set canvas size based on text
        canvas.width = textWidth + borderThickness * 2;
        canvas.height = fontsize * 1.4 + borderThickness * 2;
        // Need to reset font after resizing canvas
        context.font = fontsize + "px " + fontface;
        // Background color
        context.fillStyle = "rgba(" + backgroundColor.r + "," + backgroundColor.g + "," + backgroundColor.b + "," + backgroundColor.a + ")";
        // Border color
        context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + "," + borderColor.b + "," + borderColor.a + ")";
        context.lineWidth = borderThickness;
        // Draw background rectangle
        context.fillRect(borderThickness, borderThickness, textWidth, fontsize * 1.4);
        // Draw text
        context.fillStyle = "rgba(0, 0, 0, 1.0)";
        context.fillText(message, borderThickness, fontsize + borderThickness);
        // canvas contents will be used for a texture
        const texture = new THREE.CanvasTexture(canvas);
        texture.needsUpdate = true;
        const spriteMaterial = new THREE.SpriteMaterial({ map: texture });
        const sprite = new THREE.Sprite(spriteMaterial);
        // Scale sprite based on canvas size (arbitrarily scaled)
        sprite.scale.set(1.5, 0.75, 1);
        return sprite;
    }
    
    // Create text labels for the atoms and add to the molecule group
    const labelCarbon = makeTextSprite("C", { fontsize: 32, borderThickness: 2, backgroundColor: { r: 255, g: 255, b: 255, a: 0.8 } });
    labelCarbon.position.set(0, 0.6, 0);
    co2Group.add(labelCarbon);
    
    const labelOxygenLeft = makeTextSprite("O", { fontsize: 32, borderThickness: 2, backgroundColor: { r: 255, g: 255, b: 255, a: 0.8 } });
    labelOxygenLeft.position.set(-1, 0.6, 0);
    co2Group.add(labelOxygenLeft);
    
    const labelOxygenRight = makeTextSprite("O", { fontsize: 32, borderThickness: 2, backgroundColor: { r: 255, g: 255, b: 255, a: 0.8 } });
    labelOxygenRight.position.set(1, 0.6, 0);
    co2Group.add(labelOxygenRight);
    
    // Create an angular measurement overlay using a line and text sprite.
    const angleGroup = new THREE.Group();
    // Line from oxygen left to oxygen right passing through carbon
    const lineMaterial = new THREE.LineBasicMaterial({ color: 0xffffff });
    const lineGeometry = new THREE.BufferGeometry().setFromPoints([oxygenLeft.position, carbon.position, oxygenRight.position]);
    const angleLine = new THREE.Line(lineGeometry, lineMaterial);
    angleGroup.add(angleLine);
    
    // Add 180° text sprite above the carbon atom
    const angleLabel = makeTextSprite("180°", { fontsize: 24, borderThickness: 1, backgroundColor: { r: 0, g: 0, b: 0, a: 0.6 } });
    angleLabel.position.set(0, 0.9, 0);
    angleGroup.add(angleLabel);
    co2Group.add(angleGroup);
    
    // Create animated electron arrow overlays with ArrowHelper
    // Arrow pointing from carbon toward oxygenRight along the electron density area
    const arrowDirRight = new THREE.Vector3(1, 0.2, 0).normalize();
    const arrowRight = new THREE.ArrowHelper(arrowDirRight, new THREE.Vector3(0.2, 0.1, 0), 0.8, 0xffff00, 0.2, 0.1);
    co2Group.add(arrowRight);
    
    // Arrow pointing from carbon toward oxygenLeft
    const arrowDirLeft = new THREE.Vector3(-1, 0.2, 0).normalize();
    const arrowLeft = new THREE.ArrowHelper(arrowDirLeft, new THREE.Vector3(-0.2, 0.1, 0), 0.8, 0xffff00, 0.2, 0.1);
    co2Group.add(arrowLeft);
    
    // Create a summary overlay: a group of bullet point text sprites
    const bulletPoints = [
        "Linear geometry (180°)",
        "Double bonds with dynamic electron density",
        "Subtle vibrational bending"
    ];
    const overlayGroup = new THREE.Group();
    bulletPoints.forEach((point, index) => {
        const sprite = makeTextSprite("- " + point, { fontsize: 20, borderThickness: 1, backgroundColor: { r: 20, g:20, b:20, a:0.7 } });
        sprite.position.set(1.5, 0.8 - index * 0.4, 0);
        overlayGroup.add(sprite);
    });
    co2Group.add(overlayGroup);
    
    // Add the molecule group to the scene
    scene.add(co2Group);
    
    // Set up animation parameters on the co2Group
    co2Group.userData = {
        rotationSpeed: 0.005,    // slow overall rotation
        vibrationSpeed: 3.0,     // bending vibration speed
        vibrationAmplitude: 0.03 // bending vibration amplitude
    };
    
    // Save references to the double bonds for animation updates.
    window.co2DoubleBonds = [bondLeft, bondRight];
    
    // Create an update function for animations (to be called in your render loop)
    window.updateCO2Molecule = function(time) {
        // Overall slow rotation around the Y-axis
        co2Group.rotation.y += co2Group.userData.rotationSpeed;
        
        // Simulate subtle bending vibration (oscillation around the Z-axis)
        co2Group.rotation.z = Math.sin(time * 0.002 * co2Group.userData.vibrationSpeed) * co2Group.userData.vibrationAmplitude;
        
        // Animate double bonds: modulate thickness and color to simulate electron density changes
        window.co2DoubleBonds.forEach(function(bondGroup){
            bondGroup.userData.bonds.forEach(function(bond){
                // Calculate modulation factor from a sine wave based on time.
                const modulation = 0.5 + 0.5 * Math.sin(time * 0.005);
                // Update scale in the radial direction (x and z) to simulate thickening
                const baseRadius = bond.userData.baseRadius;
                bond.scale.x = 1 + 0.5 * modulation;
                bond.scale.z = 1 + 0.5 * modulation;
                
                // Update color shifting between the original and a highlight (e.g., yellowish)
                const color = new THREE.Color();
                color.lerpColors(new THREE.Color(0x999999), new THREE.Color(0xffff00), modulation);
                bond.material.color = color;
            });
        });
        
        // Optionally, animate arrow helpers (rotate them slightly for dynamic effect)
        arrowRight.setDirection(arrowDirRight.clone().applyAxisAngle(new THREE.Vector3(0,0,1), 0.01 * Math.sin(time*0.005)));
        arrowLeft.setDirection(arrowDirLeft.clone().applyAxisAngle(new THREE.Vector3(0,0,1), -0.01 * Math.sin(time*0.005)));
    };
    
    // Optional: For energy waveform overlay, one could add additional animated geometries here.
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
    
    // Animation created by the AnimationAgent
// Get elapsed time and delta time
const elapsedTime = clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Retrieve objects from the scene
const co2Molecule = scene.getObjectByName("co2Molecule");
const co2DoubleBonds = scene.getObjectByName("co2DoubleBonds");
const cameraObj = camera; // assuming camera is in scope

// ────────────────────────────────────────────────────────────
// At 00:00 - Introducing CO2 molecule structure
if (elapsedTime >= 0 && elapsedTime < 20) {
  // Fade background from dark (black) to a softer dark hue
  const introFraction = elapsedTime / 20;
  // Assuming scene.background is a THREE.Color instance
  if (scene.background && scene.background.isColor) {
    // Lerp from black to a subtle dark blue-gray (0x222244)
    const startColor = new THREE.Color(0x000000);
    const endColor = new THREE.Color(0x222244);
    scene.background.copy(startColor).lerp(endColor, introFraction);
  }

  // Animate the CO2 molecule appearance using the provided update function if available
  if (co2Molecule) {
    if (typeof window.updateCO2Molecule === 'function') {
      // Pass the introduction stage flag as a second parameter
      window.updateCO2Molecule(elapsedTime, "intro");
    } else {
      // Fallback: gradually increase the opacity of all meshes in the molecule
      co2Molecule.traverse(child => {
        if (child.isMesh && child.material && child.material.opacity !== undefined) {
          child.material.opacity = Math.min(1, introFraction);
          child.material.transparent = child.material.opacity < 1;
        }
      });
    }
  }
  
  // Optionally, set camera to a wider view at the start
  if (cameraObj) {
    // Start far out (e.g., position z = 8) and prepare for a zoom-in later
    cameraObj.position.lerp(new THREE.Vector3(0, 0, 8), 0.02);
    cameraObj.lookAt(0, 0, 0);
  }
}

// ────────────────────────────────────────────────────────────
// At 00:20 - Carbon center with oxygen endpoints (Camera zooms in)
else if (elapsedTime >= 20 && elapsedTime < 40) {
  if (co2Molecule) {
    if (typeof window.updateCO2Molecule === 'function') {
      // Pass stage flag "zoom" to possibly highlight atom colors and labels
      window.updateCO2Molecule(elapsedTime, "zoom");
    } else {
      // Fallback: a subtle continuous rotation of the molecule
      co2Molecule.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Smoothly zoom the camera in closer to the molecule (e.g., to z = 5)
  if (cameraObj) {
    const targetPos = new THREE.Vector3(0, 0, 5);
    cameraObj.position.lerp(targetPos, 0.02);
    cameraObj.lookAt(0, 0, 0);
  }
}

// ────────────────────────────────────────────────────────────
// At 00:40 - Visualizing double bonds clearly
else if (elapsedTime >= 40 && elapsedTime < 60) {
  // Focus on the bonding details by enhancing double bonds
  if (co2Molecule) {
    if (typeof window.updateCO2Molecule === 'function') {
      // Signal the bond focus stage (e.g., to thicken bonds and animate electron arrows)
      window.updateCO2Molecule(elapsedTime, "bonds");
    } else {
      // Fallback: apply a pulsating effect to the double bonds
      if (co2DoubleBonds) {
        // Increase scale to simulate thickening
        const scaleFactor = 1 + 0.2 * Math.sin((elapsedTime - 40) * Math.PI); // oscillates around 1
        co2DoubleBonds.scale.set(scaleFactor, scaleFactor, scaleFactor);
        // Optionally adjust color intensity (if material supports it)
        co2DoubleBonds.traverse(child => {
          if (child.isMesh && child.material && child.material.emissive) {
            const intensity = 0.5 + 0.5 * Math.sin((elapsedTime - 40) * Math.PI);
            child.material.emissive = new THREE.Color(0xff8888).multiplyScalar(intensity);
          }
        });
      }
    }
  }
}

// ────────────────────────────────────────────────────────────
// At 01:00 - Linear 180° molecular alignment (Molecule rotates)
else if (elapsedTime >= 60 && elapsedTime < 80) {
  if (co2Molecule) {
    if (typeof window.updateCO2Molecule === 'function') {
      // Signal the rotation stage; an angular measurement overlay could be handled internally
      window.updateCO2Molecule(elapsedTime, "rotate");
    } else {
      // Fallback: rotate the molecule slowly about the y-axis
      co2Molecule.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // (Optional) Overlay: An angular measurement tool can be updated here if available.
}

// ────────────────────────────────────────────────────────────
// At 01:20 - Bending vibrations illustrated (Molecular vibrations)
else if (elapsedTime >= 80 && elapsedTime < 100) {
  if (co2Molecule) {
    if (typeof window.updateCO2Molecule === 'function') {
      // Signal the vibration stage to simulate subtle bending motions and energy waveform graphics
      window.updateCO2Molecule(elapsedTime, "vibrate");
    } else {
      // Fallback: apply a slight oscillation rotation around the z-axis to simulate vibration
      co2Molecule.rotation.z = 0.05 * Math.sin((elapsedTime - 80) * 4);
    }
  }
}

// ────────────────────────────────────────────────────────────
// At 01:40 - CO2 structure summary revealed (Camera pulls back to show summary)
else if (elapsedTime >= 100 && elapsedTime < 120) {
  if (co2Molecule) {
    if (typeof window.updateCO2Molecule === 'function') {
      // Signal the summary stage for any internal overlay updates (like bullet-point summaries)
      window.updateCO2Molecule(elapsedTime, "summary");
    } else {
      // Fallback: continue a gentle rotation
      co2Molecule.rotation.y += 0.05 * deltaTime;
    }
  }
  
  // Gradually pull back the camera to display the complete structure along with any summary overlays
  if (cameraObj) {
    const pullBackPos = new THREE.Vector3(0, 0, 8);
    cameraObj.position.lerp(pullBackPos, 0.02);
    cameraObj.lookAt(0, 0, 0);
  }
  
  // (Optional) If there is a DOM overlay for the summary, fade it in here
  const summaryOverlay = document.getElementById("summaryOverlay");
  if (summaryOverlay) {
    // Fade in between 1:40 and 2:00 (100 - 120 seconds)
    const overlayAlpha = Math.min(1, (elapsedTime - 100) / 20);
    summaryOverlay.style.opacity = overlayAlpha;
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
