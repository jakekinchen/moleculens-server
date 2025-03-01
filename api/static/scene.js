// Scientific visualization: Electron Interaction with a Water Molecule
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
// Create materials
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 10 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 10 });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, shininess: 5 });

// Create geometries for atoms
const oxygenGeometry = new THREE.SphereGeometry(0.5, 32, 32);
const hydrogenGeometry = new THREE.SphereGeometry(0.25, 32, 32);

// Create a group for the water molecule
const waterMolecule = new THREE.Group();
window.waterMolecule = waterMolecule;

// Oxygen atom at the center
const oxygen = new THREE.Mesh(oxygenGeometry, oxygenMaterial);
oxygen.position.set(0, 0, 0);
waterMolecule.add(oxygen);

// Define bond length and bond angle (approximate values)
// For a water molecule: bond length ~1.0 (arbitrary unit) and bond angle ~104.45°
// Each hydrogen is offset from the oxygen by half the bond angle ~52.225° from the oxygen's horizontal axis.
const bondLength = 1.0;
const angleDeg = 52.225;
const angleRad = THREE.MathUtils.degToRad(angleDeg);

// Calculate positions for the two hydrogen atoms in the X-Y plane.
// Both hydrogen atoms will be at a distance "bondLength" from the oxygen:
const h1Pos = new THREE.Vector3(
    bondLength * Math.cos(angleRad),
    bondLength * Math.sin(angleRad),
    0
);
const h2Pos = new THREE.Vector3(
    bondLength * Math.cos(angleRad),
    -bondLength * Math.sin(angleRad),
    0
);

// Create hydrogen atoms
const hydrogen1 = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);
hydrogen1.position.copy(h1Pos);
waterMolecule.add(hydrogen1);

const hydrogen2 = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);
hydrogen2.position.copy(h2Pos);
waterMolecule.add(hydrogen2);

// Function to create a bond (cylinder) between two points
function createBond(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  const bondGeometry = new THREE.CylinderGeometry(0.08, 0.08, length, 16);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);
  
  // Position the bond at the midpoint between start and end
  bond.position.copy(start).lerp(end, 0.5);
  
  // Align the bond to point from start to end
  bond.quaternion.setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    direction.normalize()
  );
  
  return bond;
}

// Create bonds between the oxygen and the two hydrogens
const bond1 = createBond(oxygen.position, hydrogen1.position);
const bond2 = createBond(oxygen.position, hydrogen2.position);

// Add bonds to the water molecule group
waterMolecule.add(bond1, bond2);

// Add the water molecule group to the scene
scene.add(waterMolecule);

// GeometryAgent LLM-generated code
// Create a large plane geometry to serve as the background
const backgroundGeometry = new THREE.PlaneGeometry(2000, 2000);

// Create shader uniforms for the gradient colors
// "color_initial" is a subtle light blue tint and "color_final" is neutral white
const bgUniforms = {
  color1: { value: new THREE.Color(0xadd8e6) }, // subtle tint (light blue)
  color2: { value: new THREE.Color(0xffffff) }  // neutral
};

// Create a custom ShaderMaterial to generate a smooth vertical gradient
const backgroundMaterial = new THREE.ShaderMaterial({
  uniforms: bgUniforms,
  vertexShader: `
    varying vec2 vUv;
    void main() {
      vUv = uv;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }
  `,
  fragmentShader: `
    uniform vec3 color1;
    uniform vec3 color2;
    varying vec2 vUv;
    void main() {
      // Mix from color1 at the bottom (vUv.y = 0) to color2 at the top (vUv.y = 1)
      gl_FragColor = vec4(mix(color1, color2, vUv.y), 1.0);
    }
  `,
  side: THREE.DoubleSide // Ensure the background is visible from the back as well
});

// Create the mesh from the geometry and the shader material
const backgroundMesh = new THREE.Mesh(backgroundGeometry, backgroundMaterial);

// Position the background mesh back along the z-axis so it fills the entire scene
backgroundMesh.position.set(0, 0, -1000);

// Store a global reference if needed
window.backgroundMesh = backgroundMesh;

// Add the background mesh to the scene
scene.add(backgroundMesh);

// GeometryAgent LLM-generated code
// Create a group to hold the electron and its glow effect
const electronGroup = new THREE.Group();
window.electron = electronGroup; // store global reference

// Create the small electron geometry
const electronGeometry = new THREE.SphereGeometry(0.1, 32, 32);

// Create an emissive material with a cyan/light blue color for the electron
const electronMaterial = new THREE.MeshPhongMaterial({ 
  color: 0x00ffff, 
  emissive: 0x00ffff, 
  emissiveIntensity: 1.0 
});

// Create the electron mesh and add it to the group
const electronMesh = new THREE.Mesh(electronGeometry, electronMaterial);
electronGroup.add(electronMesh);

// OPTIONAL: Add a glow effect by layering a slightly larger, transparent, additive-blended sphere
const glowGeometry = new THREE.SphereGeometry(0.15, 32, 32);
const glowMaterial = new THREE.MeshBasicMaterial({
    color: 0x00ffff,
    transparent: true,
    opacity: 0.5,
    blending: THREE.AdditiveBlending,
    depthWrite: false
});
const glowMesh = new THREE.Mesh(glowGeometry, glowMaterial);
electronGroup.add(glowMesh);

// Create a curved trajectory path for the electron's motion
// The electron starts from the left and curves toward the water molecule and then deflects.
// Here we use a cubic bezier curve with chosen control points.
// Adjust control points as needed to match the precise scene choreography.
window.electronCurve = new THREE.CubicBezierCurve3(
    new THREE.Vector3(-10, 0, 0),  // start at left side of the frame
    new THREE.Vector3(-5, 2, 0),   // first control point for initial upward curvature
    new THREE.Vector3(0, -2, 0),   // second control point guiding the approach toward the water molecule
    new THREE.Vector3(2, 0, 0)     // endpoint after deflection post-interaction
);

// Optionally, one could visualize the curve for debugging:
// const curvePoints = window.electronCurve.getPoints(50);
// const curveGeometry = new THREE.BufferGeometry().setFromPoints(curvePoints);
// const curveMaterial = new THREE.LineBasicMaterial({ color: 0xff00ff });
// const curveLine = new THREE.Line(curveGeometry, curveMaterial);
// scene.add(curveLine);

// Add the electron group to the scene
scene.add(electronGroup);

// GeometryAgent LLM-generated code
// Define a curved trajectory for the electron trail using a CatmullRom curve
const trailCurve = new THREE.CatmullRomCurve3([
  new THREE.Vector3(0, 0, 0),
  new THREE.Vector3(1, 1, 0),
  new THREE.Vector3(2, 0, 1),
  new THREE.Vector3(3, -1, 0),
  new THREE.Vector3(4, 0, -1)
]);

// Create the tube geometry along the curve. 
// The parameters are: path, tubularSegments, radius, radialSegments, closed
const trailGeometry = new THREE.TubeGeometry(trailCurve, 50, 0.05, 8, false);

// Create a ShaderMaterial to simulate a glowing trail with a fading effect along its length.
// The fragment shader uses the V coordinate (vUv.y) from the tube's UV mapping 
// to reduce the opacity (fading effect) along the tube. 
const trailMaterial = new THREE.ShaderMaterial({
  uniforms: {
    glowColor: { value: new THREE.Color(0x00ffff) } // cyan/light blue color
  },
  vertexShader: `
    varying vec2 vUv;
    void main() {
      vUv = uv;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }
  `,
  fragmentShader: `
    uniform vec3 glowColor;
    varying vec2 vUv;
    void main() {
      // Fade the opacity along the tube length: full intensity at the start (vUv.y = 0)
      // and diminishing intensity at the end (vUv.y = 1)
      float alpha = 1.0 - vUv.y;
      gl_FragColor = vec4(glowColor, alpha);
    }
  `,
  transparent: true,
  blending: THREE.AdditiveBlending,
  depthWrite: false
});

// Create the mesh for the electron trail
const electronTrail = new THREE.Mesh(trailGeometry, trailMaterial);
window.electronTrail = electronTrail;

// Add the glowing electron trail to the scene
scene.add(electronTrail);

// GeometryAgent LLM-generated code
// Create a group for the electromagnetic field lines and store a global reference
const fieldLinesGroup = new THREE.Group();
window.electromagneticFieldLines = fieldLinesGroup;

// Create a thin, semi-transparent line material (light yellow)
const lineMaterial = new THREE.LineBasicMaterial({
  color: 0xffffaa,
  transparent: true,
  opacity: 0.5
});

// We'll create multiple field lines that arc between the electron and water molecule.
// For demonstration, we'll assume the electron is around (-1, 0, 0) and the water molecule is around (1, 0, 0).
// Each line will be defined by a set of vertices along a curve with a gentle perturbation applied later for animation.
const numLines = 20;
const numPoints = 20;
for (let i = 0; i < numLines; i++) {
  const positions = [];       // flat array for BufferGeometry positions
  const basePositions = [];   // store original positions to use in animation
  
  // Random parameters to vary the shape of each field line
  const theta = Math.random() * Math.PI * 2;
  const phi = Math.random() * Math.PI;
  // This radius will control the maximum deviation; vary it slightly per line.
  const radius = 0.5 + Math.random() * 0.5;
  
  // Generate numPoints along the curve
  for (let j = 0; j < numPoints; j++) {
    const t = j / (numPoints - 1);
    // Compute a base point between the electron and water molecule
    const start = new THREE.Vector3(-1, 0, 0);
    const end = new THREE.Vector3(1, 0, 0);
    const basePoint = new THREE.Vector3().lerpVectors(start, end, t);
    
    // Add a perturbation perpendicular to the line to simulate a curving electromagnetic field.
    // The offset is strongest near the middle of the line.
    const offsetStrength = radius * (1 - Math.abs(0.5 - t) * 2);
    const offset = new THREE.Vector3(
      Math.cos(theta) * Math.sin(phi),
      Math.sin(theta) * Math.sin(phi),
      Math.cos(phi)
    ).multiplyScalar(offsetStrength);
    basePoint.add(offset);
    
    // Save the vertex position
    positions.push(basePoint.x, basePoint.y, basePoint.z);
    basePositions.push(basePoint.clone());
  }
  
  // Create the BufferGeometry from the positions array
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  // Save the base positions and a random phase for animation in the geometry's userData.
  geometry.userData.basePositions = basePositions;
  geometry.userData.phase = Math.random() * Math.PI * 2;
  
  // Create the line object and add it to the group
  const line = new THREE.Line(geometry, lineMaterial);
  fieldLinesGroup.add(line);
}

// Add the group of field lines to the scene
scene.add(fieldLinesGroup);

// Define an animation function to gently pulsate the field lines,
// mimicking a fluctuating electromagnetic field.
// This function should be called from the main render loop.
window.animateFieldLines = function(time) {
  // Iterate over each field line
  fieldLinesGroup.children.forEach(function(line) {
    const geometry = line.geometry;
    const basePositions = geometry.userData.basePositions;
    const phase = geometry.userData.phase;
    const positions = geometry.attributes.position.array;
    
    // Modify each vertex to apply a small oscillation (here along the Y axis)
    for (let i = 0; i < basePositions.length; i++) {
      const basePos = basePositions[i];
      // Compute a small offset using a sine function: adjust amplitude and speed as desired.
      const offset = Math.sin(time * 0.002 + phase + i) * 0.1;
      
      positions[i * 3 + 0] = basePos.x;
      positions[i * 3 + 1] = basePos.y + offset;
      positions[i * 3 + 2] = basePos.z;
    }
    
    geometry.attributes.position.needsUpdate = true;
  });
};

// GeometryAgent LLM-generated code
// Create a group for the Electron_Cloud_Polarization effect
const electronCloudPolarization = new THREE.Group();
window.electronCloudPolarization = electronCloudPolarization;

// Create a canvas to generate a radial gradient texture (red to a warmer hue)
const canvas = document.createElement('canvas');
canvas.width = 256;
canvas.height = 256;
const ctx = canvas.getContext('2d');

const centerX = canvas.width / 2;
const centerY = canvas.height / 2;
const radius = canvas.width / 2;
const gradient = ctx.createRadialGradient(centerX, centerY, 0, centerX, centerY, radius);
gradient.addColorStop(0, 'rgba(255, 0, 0, 1)');       // solid red at the center
gradient.addColorStop(1, 'rgba(255, 153, 51, 0)');      // fading to a warmer (orange) tone at the edge

ctx.fillStyle = gradient;
ctx.fillRect(0, 0, canvas.width, canvas.height);

// Create a texture from the canvas
const gradientTexture = new THREE.CanvasTexture(canvas);
gradientTexture.needsUpdate = true;

// Create a material for the effect using the gradient texture.
// Using transparency and additive blending gives a glowing, soft effect.
const polarizationMaterial = new THREE.MeshPhongMaterial({
  map: gradientTexture,
  transparent: true,
  opacity: 0.7,
  side: THREE.DoubleSide,
  blending: THREE.AdditiveBlending,
  depthWrite: false
});

// Create a sphere geometry that will encompass the oxygen atom region.
// The radius (1.2) is chosen to slightly oversize a typical oxygen represented elsewhere.
const polarizationGeometry = new THREE.SphereGeometry(1.2, 32, 32);
const polarizationMesh = new THREE.Mesh(polarizationGeometry, polarizationMaterial);

// Position the mesh so that its center aligns with the oxygen atom's position.
// Adjust this position if your water molecule oxygen atom is located elsewhere.
polarizationMesh.position.set(0, 0, 0);

// Add the mesh to the group
electronCloudPolarization.add(polarizationMesh);

// Add the polarization group to the scene
scene.add(electronCloudPolarization);

// Store animation parameters in the group's userData for pulsation effect.
electronCloudPolarization.userData = { 
  baseScale: 1,
  amplitude: 0.1,
  speed: 2  // pulsation speed multiplier
};

// Set up an animation loop to create subtle expansion and contraction.
if (!window.electronCloudPolarizationAnimate) {
  window.electronCloudPolarizationAnimate = function () {
    requestAnimationFrame(window.electronCloudPolarizationAnimate);
    const time = performance.now() * 0.001;
    const scale = electronCloudPolarization.userData.baseScale + 
                  electronCloudPolarization.userData.amplitude * Math.sin(time * electronCloudPolarization.userData.speed);
    polarizationMesh.scale.set(scale, scale, scale);
  };
  window.electronCloudPolarizationAnimate();
}

// GeometryAgent LLM-generated code
// Create a group for the Energy Ripple Waves effect
const energyRippleGroup = new THREE.Group();
window.energyRippleWaves = energyRippleGroup;

// Create a ring geometry representing the initial ripple
const innerRadius = 0.5;
const outerRadius = 0.55;
const thetaSegments = 32;
const ringGeometry = new THREE.RingGeometry(innerRadius, outerRadius, thetaSegments);

// Create a soft, light material that fades outward
const rippleMaterial = new THREE.MeshPhongMaterial({
  color: 0xffffff,
  side: THREE.DoubleSide,
  transparent: true,
  opacity: 1,
  emissive: 0xffffff,
  emissiveIntensity: 0.5
});

// Create the ripple mesh and orient it horizontally
const rippleMesh = new THREE.Mesh(ringGeometry, rippleMaterial);
rippleMesh.rotation.x = Math.PI / 2;
energyRippleGroup.add(rippleMesh);

// Store animation parameters in the group
energyRippleGroup.userData = {
  startTime: performance.now(),
  duration: 1000, // Duration of the ripple's propagation in milliseconds
  initialScale: 1,
  finalScale: 5
};

// Add the energy ripple group to the scene
scene.add(energyRippleGroup);

// Define an update function to animate the ripple effect
window.updateEnergyRippleWaves = function() {
  const currentTime = performance.now();
  const elapsed = currentTime - energyRippleGroup.userData.startTime;
  const t = elapsed / energyRippleGroup.userData.duration;
  
  if (t < 1) {
    // Scale the ripple from its initial size to a larger size
    const scale = energyRippleGroup.userData.initialScale + t * (energyRippleGroup.userData.finalScale - energyRippleGroup.userData.initialScale);
    rippleMesh.scale.set(scale, scale, scale);
    
    // Fade out the opacity over time for a diminishing amplitude effect
    rippleMaterial.opacity = 1 - t;
  } else {
    // Once the animation is complete, remove the ripple from the scene
    scene.remove(energyRippleGroup);
    
    // Optionally reset or restart the effect here if needed
  }
};

// Note: Call window.updateEnergyRippleWaves() on each frame (e.g., within your main render/animation loop)

// GeometryAgent LLM-generated code
// Create a group for vibrational arrows and store it globally
const vibrationalArrows = new THREE.Group();
window.vibrationalArrows = vibrationalArrows;

// Define a contrasting accent color for the arrows (light blue)
const arrowColor = 0x66ccff;

// Utility function to create a vibrational arrow using THREE.ArrowHelper.
// This arrow is built along a defined hydrogen bond direction.
function createVibrationalArrow(from, to) {
  const direction = new THREE.Vector3().subVectors(to, from).normalize();
  const length = new THREE.Vector3().subVectors(to, from).length();
  
  // Create the arrow helper with a scaled head and shaft.
  // ArrowHelper internally creates two objects: a line (shaft) and a cone (arrow head)
  const arrow = new THREE.ArrowHelper(direction, from, length, arrowColor, 0.2 * length, 0.1 * length);
  
  // Save parameters for animation: original length, oscillation amplitude, and a phase offset
  arrow.userData = {
    initialLength: length,
    amplitude: 0.05 * length,
    phase: Math.random() * Math.PI * 2 // random phase for variation
  };
  
  return arrow;
}

// For the purposes of this vibrational visualization, assume the water molecule has 
// an Oxygen at the origin and two Hydrogens positioned roughly as in a water molecule.
const oxygen = new THREE.Vector3(0, 0, 0);
const hydrogen1 = new THREE.Vector3(0.95, 0.3, 0);
const hydrogen2 = new THREE.Vector3(0.95, -0.3, 0);

// Create vibrational arrows along the hydrogen bonds (from oxygen to each hydrogen)
const vibrationalArrow1 = createVibrationalArrow(oxygen, hydrogen1);
const vibrationalArrow2 = createVibrationalArrow(oxygen, hydrogen2);

// Add the arrows to the vibrational arrows group
vibrationalArrows.add(vibrationalArrow1, vibrationalArrow2);

// Add the complete vibrational arrows group to the scene
scene.add(vibrationalArrows);

// Animation update function.
// Call window.animateVibrationalArrows(time) within your render loop (time in milliseconds)
window.animateVibrationalArrows = function(time) {
  const t = time / 1000; // convert time to seconds for smoother animation
  
  vibrationalArrows.children.forEach(arrow => {
    const { initialLength, amplitude, phase } = arrow.userData;
    // Oscillatory adjustment based on a sine wave for a subtle pulsating effect
    const oscillation = Math.sin(t * 4 + phase) * amplitude; // frequency factor 4 is chosen arbitrarily
    const newLength = initialLength + oscillation;
    
    // Adjust the line (shaft) of the ArrowHelper.
    // The arrow's shaft is drawn along the local Y axis.
    if (arrow.line) {
      // Scale the Y-axis to reflect the new length
      arrow.line.scale.y = newLength / initialLength;
    }
    
    // Adjust the cone (arrow head) position along the Y axis.
    if (arrow.cone) {
      arrow.cone.position.y = newLength;
    }
  });
};
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
// Get elapsed time in seconds and deltaTime for smooth animations
const elapsedTime = clock.getElapsedTime();
const deltaTime = clock.getDelta();

// ------------------------------
// At 00:00 – Initialization of Water Molecule and Background
// ------------------------------
if (window.Water_Molecule) {
  // Ensure water molecule is centered at the beginning
  window.Water_Molecule.position.set(0, 0, 0);
  // Optionally, reset any vibration offsets if needed
}
if (window.Background && window.Background.material) {
  // Set the background to an initial subtle tint (assumed preset)
  window.Background.material.opacity = 1;
}
// (Caption “Water molecule in focus” assumed handled elsewhere)

// ------------------------------
// At 00:20 – Electron Approaching with Trail & Field Lines
// ------------------------------
if (elapsedTime >= 20 && elapsedTime < 40) {
  // Calculate normalized progress (0 to 1) for this segment
  const t = (elapsedTime - 20) / 20;
  
  // Animate Electron along a curved trajectory from the left to the water molecule (center)
  if (window.Electron) {
    // Starting far left, moving to x=0, with a curved (sinusoidal) offset in y
    window.Electron.position.x = -20 + t * 20;            // x: -20 to 0
    window.Electron.position.y = 10 * Math.sin(t * Math.PI);  // y: curve upward then downward
    window.Electron.position.z = 0;
  }
  
  // Sync the electron trail with the electron's position
  if (window.Electron_Trail) {
    window.Electron_Trail.position.copy(window.Electron.position);
    // Optionally, enhance trail opacity gradually
    if (window.Electron_Trail.material) {
      window.Electron_Trail.material.opacity = t;
    }
  }
  
  // Animate Electromagnetic Field Lines to become gradually visible
  if (window.Electromagnetic_Field_Lines && window.Electromagnetic_Field_Lines.material) {
    window.Electromagnetic_Field_Lines.material.opacity = Math.min(1, t);
  }
}
// (Caption “Electron approaching” assumed handled elsewhere)

// ------------------------------
// At 00:40 – Water Molecule Polarization Effects
// ------------------------------
if (elapsedTime >= 40 && elapsedTime < 60) {
  const t = (elapsedTime - 40) / 20;
  
  // Animate the electron cloud polarization effect on the oxygen atom
  if (window.Electron_Cloud_Polarization) {
    // Create a pulsating effect by modulating the scale between 1 and 1.2
    const scaleFactor = 1 + 0.2 * Math.sin(t * Math.PI * 4);
    window.Electron_Cloud_Polarization.scale.set(scaleFactor, scaleFactor, scaleFactor);
  }
  
  // Optional: Highlight the water molecule with a subtle color/emissive change to indicate polarization
  if (window.Water_Molecule && window.Water_Molecule.material) {
    // Increase emissive intensity gradually (assuming emissiveIntensity exists)
    window.Water_Molecule.material.emissiveIntensity = 0.5 + 0.5 * t;
  }
}
// (Caption “Molecule begins to polarize” assumed handled elsewhere)

// ------------------------------
// At 01:00 – Transient Electron Interaction and Energy Ripple Waves
// ------------------------------
if (elapsedTime >= 60 && elapsedTime < 80) {
  const t = (elapsedTime - 60) / 20;
  
  // Animate a brief interaction on the electron's trajectory by adding a slight deflection
  if (window.Electron) {
    // For example, nudge electron's x position slightly to the right during interaction
    window.Electron.position.x += 2 * deltaTime; // gentle lateral deflection
    // Optionally modify y if desired, e.g., a subtle oscillation
    window.Electron.position.y += 0.5 * Math.sin(t * Math.PI) * deltaTime;
  }
  
  // Animate Energy Ripple Waves over the water molecule to depict energy exchange
  if (window.Energy_Ripple_Waves) {
    // Expand the ripple's scale and fade its opacity over this segment
    const rippleScale = 1 + t; // scale increases from 1 to 2
    window.Energy_Ripple_Waves.scale.set(rippleScale, rippleScale, rippleScale);
    if (window.Energy_Ripple_Waves.material) {
      window.Energy_Ripple_Waves.material.opacity = 1 - t; // fade out the ripple as time progresses
    }
  }
}
// (Caption “Transient electron interaction” assumed handled elsewhere)

// ------------------------------
// At 01:20 – Energy Transfer and Molecular Vibrations
// ------------------------------
if (elapsedTime >= 80 && elapsedTime < 100) {
  const t = (elapsedTime - 80) / 20;
  
  // Animate Vibrational Arrows to simulate energy transfer along the hydrogen bonds
  if (window.Vibrational_Arrows && window.Vibrational_Arrows.material) {
    // Use a high-frequency pulsation based on elapsedTime to mimic rapid vibrations
    const vibrateIntensity = Math.abs(Math.sin(elapsedTime * 10));
    window.Vibrational_Arrows.material.opacity = vibrateIntensity;
  }
  
  // Impart small oscillatory motions to the water molecule (simulate bond vibrations)
  if (window.Water_Molecule) {
    // Apply a subtle vertical oscillation
    window.Water_Molecule.position.y = 0.1 * Math.sin(elapsedTime * 15);
  }
  
  // Optional: A slight camera shake for emphasis on the molecular reaction
  if (window.camera) {
    window.camera.position.z = 50 - 5 * Math.sin(elapsedTime * 2);
  }
}
// (Caption “Energy transfer and vibrations” assumed handled elsewhere)

// ------------------------------
// At 01:40 – Electron Exits and Equilibrium Restored
// ------------------------------
if (elapsedTime >= 100) {
  const t = (elapsedTime - 100) / 20;  // t increases beyond 0, used for smooth exit transitions
  
  // Animate the electron moving away on a gently curving trajectory
  if (window.Electron) {
    // Advance the electron further to the right over time
    window.Electron.position.x += 5 * deltaTime;
    // Create a soft curving motion in the y-direction
    window.Electron.position.y = 2 * Math.cos(t * Math.PI);
  }
  
  // Fade out the electron trail as the electron departs
  if (window.Electron_Trail) {
    window.Electron_Trail.position.copy(window.Electron.position);
    if (window.Electron_Trail.material) {
      window.Electron_Trail.material.opacity = Math.max(0, 1 - t);
    }
  }
  
  // Gradually bring the water molecule back to its equilibrium state
  if (window.Water_Molecule) {
    // Reduce any vibrational offset by damping the y-position oscillation
    window.Water_Molecule.position.y *= 0.98;
    if (Math.abs(window.Water_Molecule.position.y) < 0.001) {
      window.Water_Molecule.position.y = 0;
    }
    // Also, gently revert the emissive intensity to its normal level
    if (window.Water_Molecule.material && window.Water_Molecule.material.emissiveIntensity !== undefined) {
      window.Water_Molecule.material.emissiveIntensity *= 0.98;
    }
  }
  
  // Fade out the vibrational arrows over time
  if (window.Vibrational_Arrows && window.Vibrational_Arrows.material) {
    window.Vibrational_Arrows.material.opacity -= deltaTime;
    if (window.Vibrational_Arrows.material.opacity < 0) {
      window.Vibrational_Arrows.material.opacity = 0;
    }
  }
  
  // Restore the background to its neutral state (if any animation was applied previously)
  if (window.Background && window.Background.material) {
    // For example, gradually adjust opacity or tint to a neutral look
    window.Background.material.opacity = Math.min(1, 0.5 + t);
  }
  
  // Optionally, move the camera to a final resting position
  if (window.camera) {
    window.camera.position.set(0, 0, 50 + 5 * t);
  }
}
// (Caption “Equilibrium restored” assumed handled elsewhere)
    
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
