// Scientific visualization: Acetylene’s Triple Bond Explained
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
// Create a group for the acetylene molecule model
const acetylene = new THREE.Group();
window.acetylene = acetylene;

// Materials
const carbonMaterial = new THREE.MeshPhongMaterial({
  color: 0x555555,           // glowing dark grey
  emissive: 0x111111, 
  transparent: true,
  opacity: 0.9,
  shininess: 100
});

const sigmaMaterial = new THREE.MeshPhongMaterial({
  color: 0xffffff,           // bright white with subtle metallic sheen
  specular: 0x888888,
  shininess: 200,
  emissive: 0xaaaaaa
});

const piMaterial = new THREE.MeshPhongMaterial({
  color: 0xadd8e6,           // light blue with semi-transparent glow
  transparent: true,
  opacity: 0.7,
  emissive: 0x1111ff,
  side: THREE.DoubleSide
});

// Geometry for carbon atoms (spheres)
// Each carbon sphere should have a diameter of 2 units (radius = 1)
const carbonGeometry = new THREE.SphereGeometry(1, 32, 32);

// Position the two carbon atoms so that the bond length is approximately 3 units
const leftCarbon = new THREE.Mesh(carbonGeometry, carbonMaterial);
leftCarbon.position.set(-1.5, 0, 0);

const rightCarbon = new THREE.Mesh(carbonGeometry, carbonMaterial);
rightCarbon.position.set(1.5, 0, 0);

// Add carbon atoms to the group
acetylene.add(leftCarbon, rightCarbon);

// Create the sigma bond as a cylinder connecting the two carbons
// CylinderGeometry is created along the Y axis by default, so we rotate it to align along the X axis.
const sigmaLength = 3; // bond length
const sigmaRadius = 0.1;
const sigmaGeometry = new THREE.CylinderGeometry(sigmaRadius, sigmaRadius, sigmaLength, 32);
const sigmaBond = new THREE.Mesh(sigmaGeometry, sigmaMaterial);
sigmaBond.rotation.z = Math.PI / 2; // rotate so the cylinder extends along the X axis
// Position the sigma bond in between the two carbon atoms
sigmaBond.position.set(0, 0, 0);
acetylene.add(sigmaBond);

// Create the two pi bonds as flat, planar electron density clouds.
// They are modeled using BoxGeometry and oriented perpendicular to the bond (X axis).

// First pi bond: lying in the x-y plane
// Dimensions: length 3 (x), width 0.3 (y), and very thin thickness 0.05 (z)
const piGeometry1 = new THREE.BoxGeometry(sigmaLength, 0.3, 0.05);
const piBond1 = new THREE.Mesh(piGeometry1, piMaterial);
piBond1.position.set(0, 0, 0);
acetylene.add(piBond1);

// Second pi bond: lying in the x-z plane
// Dimensions: length 3 (x), width 0.3 (z), and very thin thickness 0.05 (y)
const piGeometry2 = new THREE.BoxGeometry(sigmaLength, 0.05, 0.3);
const piBond2 = new THREE.Mesh(piGeometry2, piMaterial);
piBond2.position.set(0, 0, 0);
acetylene.add(piBond2);

// Optional: Store separate references for the sigma and pi bonds for later dynamic animations
window.acetyleneSigmaBond = sigmaBond;
window.acetylenePiBonds = [piBond1, piBond2];

// The acetylene molecule is intended to continuously rotate and support various overlay animations.
// The animation code (such as slow rotation, electron arrow flows, zoom transitions, energy bars, etc.)
// should be implemented in the main render loop or in separate functions.

// Add the complete acetylene molecule to the scene
scene.add(acetylene);

// GeometryAgent LLM-generated code
(function() {
  // Utility: Create a bond (cylinder) between two points
  function createBond(start, end, radius, material) {
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();
    const bondGeometry = new THREE.CylinderGeometry(radius, radius, length, 16);
    const bond = new THREE.Mesh(bondGeometry, material);
    // position at midpoint
    const midpoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
    bond.position.copy(midpoint);
    // orient the bond
    bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
    return bond;
  }

  // Utility: Create a text sprite (using a canvas texture)
  function makeTextSprite(message, parameters) {
    if (parameters === undefined) parameters = {};
    var fontface = parameters.hasOwnProperty("fontface") ? parameters["fontface"] : "Arial";
    var fontsize = parameters.hasOwnProperty("fontsize") ? parameters["fontsize"] : 48;
    var borderThickness = parameters.hasOwnProperty("borderThickness") ? parameters["borderThickness"] : 4;
    var borderColor = parameters.hasOwnProperty("borderColor") ? parameters["borderColor"] : { r:0, g:0, b:0, a:1.0 };
    var backgroundColor = parameters.hasOwnProperty("backgroundColor") ? parameters["backgroundColor"] : { r:255, g:255, b:255, a:1.0 };

    var canvas = document.createElement('canvas');
    var context = canvas.getContext('2d');
    context.font = "Bold " + fontsize + "px " + fontface;
    // get size data (approx.)
    var metrics = context.measureText(message);
    var textWidth = metrics.width;
    canvas.width = textWidth + borderThickness * 2;
    canvas.height = fontsize * 1.4 + borderThickness * 2;
    // background color
    context.fillStyle = "rgba(" + backgroundColor.r + "," + backgroundColor.g + "," + backgroundColor.b + "," + backgroundColor.a + ")";
    context.fillRect(0, 0, canvas.width, canvas.height);
    // border
    context.strokeStyle = "rgba(" + borderColor.r + "," + borderColor.g + "," + borderColor.b + "," + borderColor.a + ")";
    context.lineWidth = borderThickness;
    context.strokeRect(0, 0, canvas.width, canvas.height);
    // text color
    context.fillStyle = "rgba(0, 0, 0, 1.0)";
    context.textAlign = "center";
    context.textBaseline = "middle";
    context.font = "Bold " + fontsize + "px " + fontface;
    context.fillText(message, canvas.width / 2, canvas.height / 2);
    var texture = new THREE.CanvasTexture(canvas);
    texture.minFilter = THREE.LinearFilter;
    var spriteMaterial = new THREE.SpriteMaterial({ map: texture });
    var sprite = new THREE.Sprite(spriteMaterial);
    // scale sprite based on font size
    sprite.scale.set(2, 1, 1);
    return sprite;
  }

  // Materials (high gloss, reflective, subtle transparency)
  const materialParams = { shininess: 100, transparent: true, opacity: 0.9 };
  const singleBondMaterial = new THREE.MeshPhongMaterial(Object.assign({ color: 0xadd8e6 }, materialParams)); // soft blue
  const doubleBondMaterial = new THREE.MeshPhongMaterial(Object.assign({ color: 0x00cc00 }, materialParams)); // vibrant green
  const tripleBondMaterial = new THREE.MeshPhongMaterial(Object.assign({ color: 0xff0000 }, materialParams)); // intense red

  // Glow material for triple bond glow effect
  const glowMaterial = new THREE.MeshBasicMaterial({
    color: 0xff0000,
    transparent: true,
    opacity: 0.5,
    blending: THREE.AdditiveBlending,
    depthWrite: false
  });

  // Common sphere geometry for atoms (representing carbon atoms)
  const atomGeometry = new THREE.SphereGeometry(0.3, 32, 32);
  
  // Model parameters
  const modelScale = 3;  // approx. 2-3 units in length
  const halfScale = modelScale / 2;
  const leftPos = new THREE.Vector3(-halfScale, 0, 0);
  const rightPos = new THREE.Vector3(halfScale, 0, 0);

  // Offsets for multiple bonds in double and triple models
  const offset = 0.15;

  // Create group for all bond comparison models
  const bondComparisonGroup = new THREE.Group();
  window.Bond_Comparison_Models = bondComparisonGroup;

  // 1. Single Bond Model
  const singleModel = new THREE.Group();
  // Atoms
  const singleAtom1 = new THREE.Mesh(atomGeometry, singleBondMaterial);
  singleAtom1.position.copy(leftPos);
  const singleAtom2 = new THREE.Mesh(atomGeometry, singleBondMaterial);
  singleAtom2.position.copy(rightPos);
  // Bond - single cylinder
  const singleBond = createBond(leftPos, rightPos, 0.1, singleBondMaterial);
  // Label
  const singleLabel = makeTextSprite("Single Bond", { fontsize: 32, backgroundColor: {r:255,g:255,b:255,a:0.8} });
  singleLabel.position.set(0, 1.2, 0);
  // Directional arrow (from left atom to right atom)
  const singleDir = new THREE.Vector3().subVectors(rightPos, leftPos).normalize();
  const singleArrow = new THREE.ArrowHelper(singleDir, leftPos, modelScale, 0x000000);
  
  singleModel.add(singleAtom1, singleAtom2, singleBond, singleLabel, singleArrow);
  // Position single bond model to the left
  singleModel.position.set(-5, 0, 0);
  bondComparisonGroup.add(singleModel);
  window.singleBondModel = singleModel;

  // 2. Double Bond Model
  const doubleModel = new THREE.Group();
  // Atoms
  const doubleAtom1 = new THREE.Mesh(atomGeometry, doubleBondMaterial);
  doubleAtom1.position.copy(leftPos);
  const doubleAtom2 = new THREE.Mesh(atomGeometry, doubleBondMaterial);
  doubleAtom2.position.copy(rightPos);
  // Bonds - two parallel cylinders with slight offset in z-axis
  const leftPosOffsetUp = leftPos.clone().setZ(offset);
  const rightPosOffsetUp = rightPos.clone().setZ(offset);
  const leftPosOffsetDown = leftPos.clone().setZ(-offset);
  const rightPosOffsetDown = rightPos.clone().setZ(-offset);
  const doubleBond1 = createBond(leftPosOffsetUp, rightPosOffsetUp, 0.08, doubleBondMaterial);
  const doubleBond2 = createBond(leftPosOffsetDown, rightPosOffsetDown, 0.08, doubleBondMaterial);
  // Label
  const doubleLabel = makeTextSprite("Double Bond", { fontsize: 32, backgroundColor: {r:255,g:255,b:255,a:0.8} });
  doubleLabel.position.set(0, 1.2, 0);
  // Directional arrow
  const doubleArrow = new THREE.ArrowHelper(singleDir, leftPos, modelScale, 0x000000);
  
  doubleModel.add(doubleAtom1, doubleAtom2, doubleBond1, doubleBond2, doubleLabel, doubleArrow);
  // Position double bond model center
  doubleModel.position.set(0, 0, 0);
  bondComparisonGroup.add(doubleModel);
  window.doubleBondModel = doubleModel;

  // 3. Triple Bond Model
  const tripleModel = new THREE.Group();
  // Atoms
  const tripleAtom1 = new THREE.Mesh(atomGeometry, tripleBondMaterial);
  tripleAtom1.position.copy(leftPos);
  const tripleAtom2 = new THREE.Mesh(atomGeometry, tripleBondMaterial);
  tripleAtom2.position.copy(rightPos);
  // Bonds - three cylinders with offsets in y-direction
  const leftPosTop = leftPos.clone().setY(offset);
  const rightPosTop = rightPos.clone().setY(offset);
  const leftPosMid = leftPos.clone();
  const rightPosMid = rightPos.clone();
  const leftPosBot = leftPos.clone().setY(-offset);
  const rightPosBot = rightPos.clone().setY(-offset);
  const tripleBondTop = createBond(leftPosTop, rightPosTop, 0.07, tripleBondMaterial);
  const tripleBondMid = createBond(leftPosMid, rightPosMid, 0.07, tripleBondMaterial);
  const tripleBondBot = createBond(leftPosBot, rightPosBot, 0.07, tripleBondMaterial);
  // Glow effect added to the middle bond: overlay a slightly larger, pulsating cylinder.
  const glowGeometry = new THREE.CylinderGeometry(0.09, 0.09, modelScale, 16);
  const tripleGlowMesh = new THREE.Mesh(glowGeometry, glowMaterial);
  // Align the glow with the bond direction (default cylinder is along Y) so we perform same rotation as bond
  tripleGlowMesh.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), new THREE.Vector3().subVectors(rightPosMid, leftPosMid).normalize());
  // Position the glow at the midpoint
  const midPoint = new THREE.Vector3().addVectors(leftPosMid, rightPosMid).multiplyScalar(0.5);
  tripleGlowMesh.position.copy(midPoint);
  window.tripleGlowMesh = tripleGlowMesh; // store for animation pulsation

  // Label
  const tripleLabel = makeTextSprite("Triple Bond", { fontsize: 32, backgroundColor: {r:255,g:255,b:255,a:0.8} });
  tripleLabel.position.set(0, 1.2, 0);
  // Directional arrow
  const tripleArrow = new THREE.ArrowHelper(singleDir, leftPos, modelScale, 0x000000);
  
  tripleModel.add(tripleAtom1, tripleAtom2, tripleBondTop, tripleBondMid, tripleBondBot, tripleGlowMesh, tripleLabel, tripleArrow);
  // Position triple bond model to the right
  tripleModel.position.set(5, 0, 0);
  bondComparisonGroup.add(tripleModel);
  window.tripleBondModel = tripleModel;

  // Add the complete bond comparison group to the global scene
  scene.add(bondComparisonGroup);

  // Continuous rotation animation (to be called externally in the render loop)
  window.BondComparisonModelsAnimate = function() {
    const time = Date.now() * 0.005;
    singleModel.rotation.y = time;
    doubleModel.rotation.y = time;
    tripleModel.rotation.y = time;
    // Pulsate the triple bond glow (scale oscillation)
    if (window.tripleGlowMesh) {
      const scaleFactor = 1 + 0.1 * Math.sin(time * 3);
      tripleGlowMesh.scale.set(scaleFactor, scaleFactor, scaleFactor);
    }
  };
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
// Get elapsed time in seconds
const elapsedTime = clock.getElapsedTime();

// Use deltaTime for smooth animations (dt is assumed to be provided as an argument)
// Note: Replace "deltaTime" with your actual deltaTime variable if named differently

// -------------------
// At 00:00 - Scene Opening & Acetylene Molecule Appearance
// -------------------
if (elapsedTime < 20) {
  // Animate smooth camera pan from an initial offset to a central view
  if (window.camera) {
    // For example, move camera.x from -5 to 0 over 20 seconds and adjust z position slightly
    const panProgress = elapsedTime / 20;
    window.camera.position.x = -5 * (1 - panProgress); // from -5 to 0
    window.camera.position.z = 10 - 5 * panProgress; // from 10 to 5
    window.camera.lookAt(0, 0, 0);
  }
  
  // Animate the appearance of the acetylene molecule (scale up and glow effect)
  if (window.Acetylene_Molecule_Model) {
    // For a smooth "pop-in" during the first 2 seconds
    if (elapsedTime < 2) {
      const scaleProgress = elapsedTime / 2;
      window.Acetylene_Molecule_Model.scale.set(scaleProgress, scaleProgress, scaleProgress);
      // If the model has a material with emissive properties, brighten it gradually
      if (window.Acetylene_Molecule_Model.material) {
        window.Acetylene_Molecule_Model.material.emissiveIntensity = scaleProgress;
      }
    } else {
      // Ensure full size after the initial appearance
      window.Acetylene_Molecule_Model.scale.set(1, 1, 1);
    }
    // (Additional caption text "Meet the acetylene molecule" can be handled by UI overlays)
  }
}

// -------------------
// At 00:20 - Zoom into Bond Region & Animate Electron Clouds
// -------------------
if (elapsedTime >= 20 && elapsedTime < 40) {
  const progress = (elapsedTime - 20) / 20; // 0 to 1 as time goes from 20s to 40s
  
  // Camera zoom into the bond connection by moving along the z-axis
  if (window.camera) {
    // Zoom from z = 10 to z = 5 gradually
    window.camera.position.z = 10 - 5 * progress;
    window.camera.lookAt(0, 0, 0);
  }
  
  // Animate electron density clouds in the acetylene molecule's triple bond:
  // Assume the model has child objects for sigma and pi bonds.
  if (window.Acetylene_Molecule_Model) {
    // Sigma bond (cylindrical electron cloud) pulsing effect
    if (window.Acetylene_Molecule_Model.tripleBondSigma) {
      const pulse = 1 + 0.2 * Math.sin(progress * Math.PI * 4);
      window.Acetylene_Molecule_Model.tripleBondSigma.scale.set(pulse, pulse, pulse);
    }
    // Pi bonds (planar electron clouds) rotating slightly to simulate electron overlap
    if (window.Acetylene_Molecule_Model.tripleBondPi1) {
      window.Acetylene_Molecule_Model.tripleBondPi1.rotation.z += 0.5 * deltaTime;
    }
    if (window.Acetylene_Molecule_Model.tripleBondPi2) {
      window.Acetylene_Molecule_Model.tripleBondPi2.rotation.z -= 0.5 * deltaTime;
    }
    // Optionally animate flowing arrows if available
    if (window.Acetylene_Molecule_Model.electronArrows) {
      window.Acetylene_Molecule_Model.electronArrows.rotation.y += 0.5 * deltaTime;
    }
  }
  // (UI overlay caption "Visualizing the triple bond" should be handled separately)
}

// -------------------
// At 00:40 - Split-Screen View: Sigma vs. Pi Bonds
// -------------------
if (elapsedTime >= 40 && elapsedTime < 60) {
  const progress = (elapsedTime - 40) / 20; // 0 to 1 over 20 seconds
  
  // Rearrange parts of the molecule to emphasize split-screen:
  if (window.Acetylene_Molecule_Model) {
    // Move sigma bond to the left
    if (window.Acetylene_Molecule_Model.tripleBondSigma) {
      window.Acetylene_Molecule_Model.tripleBondSigma.position.x = -2 * progress;
    }
    // Move one pi bond to the right and upward
    if (window.Acetylene_Molecule_Model.tripleBondPi1) {
      window.Acetylene_Molecule_Model.tripleBondPi1.position.x = 2 * progress;
      window.Acetylene_Molecule_Model.tripleBondPi1.position.y = 1 * progress;
    }
    // Move the other pi bond to the right and downward
    if (window.Acetylene_Molecule_Model.tripleBondPi2) {
      window.Acetylene_Molecule_Model.tripleBondPi2.position.x = 2 * progress;
      window.Acetylene_Molecule_Model.tripleBondPi2.position.y = -1 * progress;
    }
    // Rotate animated arrows to highlight electron motion (if present)
    if (window.Acetylene_Molecule_Model.electronArrows) {
      window.Acetylene_Molecule_Model.electronArrows.rotation.z += 0.5 * deltaTime;
    }
  }
  // (Caption "Breaking down sigma & pi bonds" is managed by the UI overlay)
}

// -------------------
// At 01:00 - Side View with Distance Markers and Energy Bars
// -------------------
if (elapsedTime >= 60 && elapsedTime < 80) {
  const progress = (elapsedTime - 60) / 20; // 0 to 1 over this interval
  
  // Transition camera to a side view
  if (window.camera) {
    // Move the camera from its current position to a more lateral (side) position
    window.camera.position.x = 0 + 8 * progress;  // from 0 to 8 on the x-axis
    window.camera.position.z = 5 + 3 * progress;  // from 5 to 8 on the z-axis
    window.camera.lookAt(0, 0, 0);
  }
  
  // Animate dynamic overlays: distance markers and energy bars
  if (window.Acetylene_Molecule_Model) {
    if (window.Acetylene_Molecule_Model.distanceMarkers && window.Acetylene_Molecule_Model.distanceMarkers.material) {
      window.Acetylene_Molecule_Model.distanceMarkers.material.opacity = progress;  // Fade in markers
    }
    if (window.Acetylene_Molecule_Model.energyBars) {
      // Increase the scale on y-axis to visualize stronger bond energy
      window.Acetylene_Molecule_Model.energyBars.scale.y = 1 + progress; 
    }
  }
  // (Caption "Bond length and energy" to be handled via overlay)
}

// -------------------
// At 01:20 - Interactive Comparison: Rotating Bond Models
// -------------------
if (elapsedTime >= 80 && elapsedTime < 100) {
  // Rotate all bond comparison models for visual analysis
  if (window.Bond_Comparison_Models) {
    window.Bond_Comparison_Models.rotation.y += 0.5 * deltaTime;
    
    // Emphasize the triple bond model with glow and pulsation.
    if (window.Bond_Comparison_Models.tripleBond) {
      // Pulsate glow with a sinusoidal function (emissive intensity modulation)
      const glow = 1 + 0.5 * Math.sin(elapsedTime * 5);
      if (window.Bond_Comparison_Models.tripleBond.material) {
        window.Bond_Comparison_Models.tripleBond.material.emissiveIntensity = glow;
      }
      // Slight pulsation of its scale
      const scaleFactor = 1 + 0.1 * Math.sin(elapsedTime * 5);
      window.Bond_Comparison_Models.tripleBond.scale.set(scaleFactor, scaleFactor, scaleFactor);
    }
  }
  // (Caption "Comparing bond strength" is added via UI)
}

// -------------------
// At 01:40 - Demonstrating Reactivity with Energy Waves & Collisions
// -------------------
if (elapsedTime >= 100 && elapsedTime < 120) {
  const progress = (elapsedTime - 100) / 20;
  
  if (window.Acetylene_Molecule_Model) {
    // Animate energy waves along the triple bond (assuming an energyWave mesh exists)
    if (window.Acetylene_Molecule_Model.energyWave) {
      // Create a pulsating opacity effect
      window.Acetylene_Molecule_Model.energyWave.material.opacity = Math.abs(Math.sin(elapsedTime * 3));
      // Expand the scale to simulate a traveling energy pulse
      const scalePulse = 1 + progress;
      window.Acetylene_Molecule_Model.energyWave.scale.set(scalePulse, scalePulse, scalePulse);
    }
    // Animate collision simulation: if there are particles or indicators, spin them
    if (window.Acetylene_Molecule_Model.collisionParticles) {
      window.Acetylene_Molecule_Model.collisionParticles.rotation.x += 0.5 * deltaTime;
      window.Acetylene_Molecule_Model.collisionParticles.rotation.y += 0.5 * deltaTime;
    }
  }
  // (Caption "Reactivity and energy dynamics" is displayed via overlay)
}

// -------------------
// At 02:00 - Conclusion: Full Molecule Rotation & Camera Pan-Out
// -------------------
if (elapsedTime >= 120) {
  // Gently rotate the acetylene molecule
  if (window.Acetylene_Molecule_Model) {
    window.Acetylene_Molecule_Model.rotation.y += 0.2 * deltaTime;
  }
  
  // Slowly pan the camera out (increase the z-distance)
  if (window.camera) {
    window.camera.position.z += 0.5 * deltaTime;
    window.camera.lookAt(0, 0, 0);
  }
  
  // Optional: fade-out effect using an overlay if available
  if (window.sceneFade && window.sceneFade.material) {
    // Increase the opacity gradually for a soft fade-out
    window.sceneFade.material.opacity = Math.min(window.sceneFade.material.opacity + 0.01 * deltaTime, 1);
  }
  // (Caption "Summary of acetylene’s bond" would be rendered as a text overlay)
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
