<<<<<<< HEAD
// Scientific visualization: Electron Interaction with a Water Molecule
=======
// Scientific visualization: Exploring Ethanol: A Molecular Journey
>>>>>>> 30fdb31ebc8832c4e9ed46ccd352c293bdb5e5c3
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
<<<<<<< HEAD
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
=======
const textureLoader = new THREE.TextureLoader();
textureLoader.load('textures/starry-sky-pattern.jpg', function (starTexture) {
  // Create a material with a deep black base and the starry sky texture.
  const starMaterial = new THREE.MeshBasicMaterial({
    color: 0x000000,           // deep black base color
    map: starTexture,          // starry sky pattern overlay
    side: THREE.BackSide       // so the texture is visible from inside the sphere
  });

  // Create a large sphere geometry to serve as the background.
  // The sphere completely surrounds the scene.
  const sphereGeometry = new THREE.SphereGeometry(1000, 64, 64);

  // Create the mesh for the dark, starry background.
  window.darkStarryBackground = new THREE.Mesh(sphereGeometry, starMaterial);

  // Add the background to the scene.
  scene.add(window.darkStarryBackground);
});

// Optionally, create subtle twinkling stars using a Points system.
const starCount = 1000;
const starPositions = new Float32Array(starCount * 3);
for (let i = 0; i < starCount; i++) {
  const radius = 980; // slightly less than the sphere's radius
  const theta = THREE.MathUtils.randFloatSpread(360);
  const phi = THREE.MathUtils.randFloatSpread(360);
  
  // Convert spherical to Cartesian coordinates.
  starPositions[i * 3] = radius * Math.sin(theta) * Math.cos(phi);
  starPositions[i * 3 + 1] = radius * Math.sin(theta) * Math.sin(phi);
  starPositions[i * 3 + 2] = radius * Math.cos(theta);
}

const starsGeometry = new THREE.BufferGeometry();
starsGeometry.setAttribute('position', new THREE.BufferAttribute(starPositions, 3));

const starsMaterial = new THREE.PointsMaterial({
  color: 0xffffff,    // subtle white twinkling stars
  size: 1,
  transparent: true,
  opacity: 0.8
});

// Create the Points object and add it to the scene.
window.starryTwinkles = new THREE.Points(starsGeometry, starsMaterial);
scene.add(window.starryTwinkles);

// GeometryAgent LLM-generated code
// Create atom materials with smooth, reflective surfaces
const carbonMaterial = new THREE.MeshPhongMaterial({ color: 0x333333, shininess: 100 });
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000, shininess: 100 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 100 });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999, shininess: 100 });

// Create a base sphere geometry for all atoms
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create a group for the ethanol molecule model and store it globally
const ethanolMolecule = new THREE.Group();
window.ethanolMolecule = ethanolMolecule;

// Carbon atoms (C1 and C2)
const c1 = new THREE.Mesh(atomGeometry, carbonMaterial);
c1.position.set(0, 0, 0);
c1.scale.set(0.5, 0.5, 0.5);

const c2 = new THREE.Mesh(atomGeometry, carbonMaterial);
c2.position.set(1.5, 0, 0);
c2.scale.set(0.5, 0.5, 0.5);

// Oxygen atom (O)
const o = new THREE.Mesh(atomGeometry, oxygenMaterial);
o.position.set(3, 0, 0);
o.scale.set(0.55, 0.55, 0.55);

// Hydrogen atoms
// Attached to C1 (CH3 group)
const h1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h1.position.set(-0.8, 0.8, 0);
h1.scale.set(0.3, 0.3, 0.3);

const h2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h2.position.set(-0.8, -0.8, 0);
h2.scale.set(0.3, 0.3, 0.3);

const h3 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h3.position.set(0, 0, -0.8);
h3.scale.set(0.3, 0.3, 0.3);

// Attached to C2 (CH2 group)
const h4 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h4.position.set(1.5, 0.8, 0);
h4.scale.set(0.3, 0.3, 0.3);

const h5 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h5.position.set(1.5, -0.8, 0);
h5.scale.set(0.3, 0.3, 0.3);

// Attached to Oxygen (OH group)
const h6 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
h6.position.set(3.8, 0, 0);
h6.scale.set(0.3, 0.3, 0.3);

// Function to create a bond (cylinder) between two points
function createBond(start, end) {
  // Compute direction and length between the two atom centers
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();

  // Create a cylinder geometry for the bond; adjust radius as needed
  const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);

  // Position bond midway between start and end
  bond.position.copy(start);
  bond.position.lerp(end, 0.5);

  // Align the bond with the vector between atoms
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());

  return bond;
}

// Create bonds between atoms
const bonds = [];
bonds.push(createBond(c1.position, c2.position)); // C1-C2 bond
bonds.push(createBond(c2.position, o.position));  // C2-O bond
bonds.push(createBond(c1.position, h1.position));   // C1-H bond
bonds.push(createBond(c1.position, h2.position));   // C1-H bond
bonds.push(createBond(c1.position, h3.position));   // C1-H bond
bonds.push(createBond(c2.position, h4.position));   // C2-H bond
bonds.push(createBond(c2.position, h5.position));   // C2-H bond
bonds.push(createBond(o.position, h6.position));    // O-H bond

// Add all atoms and bonds to the ethanol molecule group
ethanolMolecule.add(c1, c2, o, h1, h2, h3, h4, h5, h6);
bonds.forEach(bond => ethanolMolecule.add(bond));

// Scale the molecule appropriately and center it
ethanolMolecule.scale.set(1.2, 1.2, 1.2);

// Add the ethanol molecule group to the global scene
scene.add(ethanolMolecule);

// Optional: Animate the molecule with slow rotation and simulate a zoom effect over time
function animateEthanol() {
  requestAnimationFrame(animateEthanol);

  // Slow rotation about the Y-axis
  ethanolMolecule.rotation.y += 0.005;

  // Optional zoom effect: pulsate the scale slightly to simulate zoom dynamics
  const scaleFactor = 1 + 0.05 * Math.sin(Date.now() * 0.001);
  ethanolMolecule.scale.set(1.2 * scaleFactor, 1.2 * scaleFactor, 1.2 * scaleFactor);
}
animateEthanol();

// GeometryAgent LLM-generated code
// Create a material for the carbon atom with a grey matte finish
const carbonMaterial = new THREE.MeshLambertMaterial({ color: 0x333333 });

// Create a sphere geometry with a radius of 0.5 and sufficient detail
const carbonGeometry = new THREE.SphereGeometry(0.5, 32, 32);

// Create the carbon atom mesh
const carbonAtom = new THREE.Mesh(carbonGeometry, carbonMaterial);
carbonAtom.name = "Carbon_Atom";

// Create a group for the ethanol molecule if not already present
const ethanolMolecule = new THREE.Group();
ethanolMolecule.name = "Ethanol_Molecule_Model";
window.Ethanol_Molecule_Model = ethanolMolecule;

// (Optionally, if you have labels or need to reference the carbon atom later)
window.Carbon_Atom = carbonAtom;

// Position the carbon atom as needed within the molecule (here positioned at the group origin)
carbonAtom.position.set(0, 0, 0);

// Add the carbon atom mesh to the molecule group
ethanolMolecule.add(carbonAtom);

// Finally, add the molecule group to the scene
scene.add(ethanolMolecule);

// GeometryAgent LLM-generated code
// Create a group for the hydrogen atom(s)
const hydrogenGroup = new THREE.Group();
window.hydrogenAtomGroup = hydrogenGroup;

// Create hydrogen atom geometry and material (matte finish)
const hydrogenGeometry = new THREE.SphereGeometry(0.3, 32, 32);
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff, shininess: 0 });

// Create the hydrogen atom mesh
const hydrogenAtom = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);

// Optionally, set the hydrogen atom position if needed (default at origin)
hydrogenAtom.position.set(0, 0, 0);

// Add the hydrogen atom mesh to the group
hydrogenGroup.add(hydrogenAtom);

// Relationships: This hydrogen atom is part of Ethanol_Molecule_Model and Atom_Labels groups.
// (Additional logic for these relationships can be added elsewhere in the scene setup)

// Add the hydrogen group to the scene
scene.add(hydrogenGroup);

// GeometryAgent LLM-generated code
// Create a group for the oxygen atom
const oxygenAtomGroup = new THREE.Group();
window.oxygenAtomGroup = oxygenAtomGroup;

// Create a matte material for the oxygen atom using MeshLambertMaterial
const oxygenMaterial = new THREE.MeshLambertMaterial({ color: 0xff0000 });

// Create the sphere geometry for the oxygen atom (radius: 0.4 units)
const oxygenGeometry = new THREE.SphereGeometry(0.4, 32, 32);

// Create the oxygen atom mesh
const oxygenMesh = new THREE.Mesh(oxygenGeometry, oxygenMaterial);
oxygenMesh.name = "Oxygen_Atom";

// Add the oxygen atom mesh to its group
oxygenAtomGroup.add(oxygenMesh);

// Add the oxygen atom group to the scene
scene.add(oxygenAtomGroup);

// GeometryAgent LLM-generated code
const loader = new THREE.FontLoader();
loader.load('fonts/clear_sans_serif.json', function (font) {
  // Create a reusable base material for the labels with initial opacity 0.
  const baseTextMaterial = new THREE.MeshPhongMaterial({
    color: 0xffffff, // light white / bright for visibility
    transparent: true,
    opacity: 0,
  });

  // Function to create a text label near an atom.
  function createLabel(text, pos) {
    const textGeometry = new THREE.TextGeometry(text, {
      font: font,
      size: 0.3, // small font size
      height: 0.05,
      curveSegments: 12,
      bevelEnabled: false,
    });
    textGeometry.computeBoundingBox();

    // Create a mesh with a cloned material so each label animates independently.
    const textMesh = new THREE.Mesh(textGeometry, baseTextMaterial.clone());
    // Position the label slightly offset from the atom position.
    textMesh.position.copy(pos).add(new THREE.Vector3(0.7, 0.7, 0));
    return textMesh;
  }

  // Create a group to hold all atom labels.
  const atomLabels = new THREE.Group();
  window.atomLabels = atomLabels; // store globally if needed

  // Define positions corresponding to each atom from the ethanol molecule.
  const c1Pos = new THREE.Vector3(0, 0, 0);
  const c2Pos = new THREE.Vector3(1.5, 0, 0);
  const oPos  = new THREE.Vector3(3, 0, 0);
  const h1Pos = new THREE.Vector3(-0.8, 0.8, 0);
  const h2Pos = new THREE.Vector3(-0.8, -0.8, 0);
  const h3Pos = new THREE.Vector3(0, 0, -0.8);
  const h4Pos = new THREE.Vector3(1.5, 0.8, 0);
  const h5Pos = new THREE.Vector3(1.5, -0.8, 0);
  const h6Pos = new THREE.Vector3(3.8, 0, 0);

  // Create and add labels to denote the chemical symbols.
  atomLabels.add(createLabel("C", c1Pos));
  atomLabels.add(createLabel("C", c2Pos));
  atomLabels.add(createLabel("O", oPos));
  atomLabels.add(createLabel("H", h1Pos));
  atomLabels.add(createLabel("H", h2Pos));
  atomLabels.add(createLabel("H", h3Pos));
  atomLabels.add(createLabel("H", h4Pos));
  atomLabels.add(createLabel("H", h5Pos));
  atomLabels.add(createLabel("H", h6Pos));

  // Function to animate a subtle fade-in of the labels.
  function animateFadeIn() {
    atomLabels.traverse(function (child) {
      if (child.material && child.material.opacity !== undefined) {
        child.material.opacity += 0.01;
        if (child.material.opacity > 1) child.material.opacity = 1;
      }
    });
    // Continue animating as long as at least one label hasn't reached full opacity.
    let needsAnimation = false;
    atomLabels.traverse(function (child) {
      if (child.material && child.material.opacity < 1) {
        needsAnimation = true;
      }
    });
    if (needsAnimation) {
      requestAnimationFrame(animateFadeIn);
    }
  }

  // The labels are intended to appear at 00:15.
  // Wait 15 seconds, then add the labels to the scene and start the fade-in animation.
  setTimeout(function () {
    scene.add(atomLabels);
    animateFadeIn();
  }, 15000);
});

// GeometryAgent LLM-generated code
// Define the overlay group for the bond angle overlay and store it globally
const bondAngleOverlay = new THREE.Group();
window.Bond_Angle_Overlay = bondAngleOverlay;

// Create a semi‐transparent yellow material for overlay lines and markers
const overlayMaterial = new THREE.LineBasicMaterial({
  color: 0xffff00,
  linewidth: 2,           // Note: linewidth may not work on all systems
  transparent: true,
  opacity: 0              // start fully transparent for a smooth fade‐in
});

// For the overlay we need the positions of three key atoms from the ethanol model:
//   • First carbon (c1)    – reference for the bond angle
//   • Second carbon (c2)   – vertex of the angle
//   • Oxygen atom (O)      – one end of the angle
//
// (These positions should be consistent with the Ethanol_Molecule_Model in the scene.)
// For this overlay we assume:
const c1Pos = new THREE.Vector3(0, 0, 0);
const c2Pos = new THREE.Vector3(1.5, 0, 0);
// Here we adjust the oxygen position slightly so the bond angle is non‑linear.
const oPos  = new THREE.Vector3(2.8, 0.8, 0);

// 1. Draw an enhanced line along the C2–O bond to emphasize the bond
const bondLineGeometry = new THREE.BufferGeometry().setFromPoints([c2Pos, oPos]);
const bondLine = new THREE.Line(bondLineGeometry, overlayMaterial);
bondAngleOverlay.add(bondLine);

// 2. Create an arc to mark the bond angle at the second carbon.
//    Typically, a bond angle marker is drawn at the vertex (second carbon) between the bonds
//    going to the first carbon and to the oxygen atom.
// Compute the two vectors from c2:
const v1 = new THREE.Vector3().subVectors(c1Pos, c2Pos);  // from c2 to c1
const v2 = new THREE.Vector3().subVectors(oPos, c2Pos);     // from c2 to O

// Get the angle (in radians) between v1 and v2
const deltaAngle = v1.angleTo(v2);

// Get the starting angle from v1 (using atan2)
const startAngle = Math.atan2(v1.y, v1.x);
// The ending angle is:
const endAngle = startAngle + deltaAngle;

// Choose a radius for the arc (adjust as needed)
const radius = 0.5;

// Generate points along the arc
const arcPoints = [];
const segments = 32;
for (let i = 0; i <= segments; i++) {
  // Linearly interpolate the angle from startAngle to endAngle
  const theta = startAngle + (deltaAngle * i / segments);
  const x = c2Pos.x + radius * Math.cos(theta);
  const y = c2Pos.y + radius * Math.sin(theta);
  arcPoints.push(new THREE.Vector3(x, y, 0));
}
const arcGeometry = new THREE.BufferGeometry().setFromPoints(arcPoints);
const arcLine = new THREE.Line(arcGeometry, overlayMaterial);
bondAngleOverlay.add(arcLine);

// 3. Add small tick markers at the arc endpoints to clearly indicate the angle boundaries.
// Function to create a tick marker given a center point and an angle (to compute a perpendicular)
function createTickMarker(center, baseAngle, length = 0.1) {
  // Compute a tangent (perpendicular) vector: rotate the unit vector at baseAngle by 90°
  const tangent = new THREE.Vector3(-Math.sin(baseAngle), Math.cos(baseAngle), 0);
  // Calculate the two endpoints of the tick line
  const p1 = new THREE.Vector3().copy(center).addScaledVector(tangent, length);
  const p2 = new THREE.Vector3().copy(center).addScaledVector(tangent, -length);
  const tickGeom = new THREE.BufferGeometry().setFromPoints([p1, p2]);
  return new THREE.Line(tickGeom, overlayMaterial);
}

// Determine the positions of the arc endpoints
const arcStartPos = arcPoints[0];            // corresponds to direction v1 (from c2)
const arcEndPos = arcPoints[arcPoints.length - 1];  // corresponds to direction v2

// Create tick markers at the start and end of the arc.
const tickStart = createTickMarker(arcStartPos, startAngle);
const tickEnd = createTickMarker(arcEndPos, endAngle);

bondAngleOverlay.add(tickStart);
bondAngleOverlay.add(tickEnd);

// Add the overlay group to the global scene
scene.add(bondAngleOverlay);

// 4. (Optional) Animate a smooth fade‐in transition for the overlay.
// The overlay is static but appears with a smooth transition by gradually increasing opacity.
const targetOpacity = 0.8;   // final opacity for semi-transparency
const fadeDuration = 1000;   // duration in milliseconds
const startTime = performance.now();

function fadeIn() {
  const elapsed = performance.now() - startTime;
  const fraction = Math.min(elapsed / fadeDuration, 1);
  overlayMaterial.opacity = fraction * targetOpacity;
  if (fraction < 1) {
    requestAnimationFrame(fadeIn);
  }
}
fadeIn();

// GeometryAgent LLM-generated code
// Create a group for the Hydroxyl Glow Effect
const hydroxylGlowGroup = new THREE.Group();
window.hydroxylGlowEffect = hydroxylGlowGroup;

// The glow is meant to highlight the hydroxyl (OH) group.
// In our ethanol molecule example the oxygen atom is at (3, 0, 0) and its bonded hydrogen at (3.8, 0, 0).
// We choose the center of the glow to lie midway between those positions.
const glowCenter = new THREE.Vector3(3.4, 0, 0);

// Create a sphere geometry that will serve as the glowing outline.
// The radius is slightly larger than the atoms (set to 1.2) so it surrounds the OH group.
const glowGeometry = new THREE.SphereGeometry(1.2, 32, 32);

// Create a glowing material with bright, warm red color, moderate transparency,
// and additive blending to simulate a light-emitting effect.
const glowMaterial = new THREE.MeshBasicMaterial({
  color: 0xff4500, // vivid red / warm glow
  transparent: true,
  opacity: 0.5,
  blending: THREE.AdditiveBlending,
  depthWrite: false
});

// Create the glow mesh and position it at the calculated center.
const glowMesh = new THREE.Mesh(glowGeometry, glowMaterial);
glowMesh.position.copy(glowCenter);
hydroxylGlowGroup.add(glowMesh);

// Add the glow group to the global scene.
scene.add(hydroxylGlowGroup);

// Store pulsation parameters and create a simple pulsating animation effect.
// To animate, call window.animateHydroxylGlow(delta) inside your render loop,
// where delta is the elapsed time since the last frame.
hydroxylGlowGroup.userData.phase = 0;
window.animateHydroxylGlow = function(delta) {
  // Increase the phase based on elapsed time and desired speed.
  hydroxylGlowGroup.userData.phase += delta * 2.0;
  const scaleFactor = 1 + 0.15 * Math.sin(hydroxylGlowGroup.userData.phase);
  glowMesh.scale.set(scaleFactor, scaleFactor, scaleFactor);
};

// GeometryAgent LLM-generated code
// Create a gradient texture using a canvas
const canvas = document.createElement('canvas');
canvas.width = 256;
canvas.height = 16;
const context = canvas.getContext('2d');
// Create a linear gradient from left (positive hue) to right (negative hue)
const gradient = context.createLinearGradient(0, 0, canvas.width, 0);
gradient.addColorStop(0, "#ff4500"); // positive (e.g. OrangeRed)
gradient.addColorStop(1, "#1E90FF"); // negative (e.g. DodgerBlue)
context.fillStyle = gradient;
context.fillRect(0, 0, canvas.width, canvas.height);
// Create texture from canvas
const gradientTexture = new THREE.CanvasTexture(canvas);
gradientTexture.needsUpdate = true;

// Create a Phong material using the gradient texture.
// Using vertexColors would require custom geometry so we use the texture.
const arrowMaterial = new THREE.MeshPhongMaterial({
  map: gradientTexture,
  transparent: true
});

// Function to create an arrow mesh from shaft and head
function createArrow() {
  const arrowGroup = new THREE.Group();

  // --- Create Arrow Shaft ---
  // CylinderGeometry is by default centered at y=0. We want it to start at y=0, so translate it upward by half its height.
  const shaftHeight = 0.5;
  const shaftGeometry = new THREE.CylinderGeometry(0.05, 0.05, shaftHeight, 12, 1);
  shaftGeometry.translate(0, shaftHeight / 2, 0);
  const shaftMesh = new THREE.Mesh(shaftGeometry, arrowMaterial);
  arrowGroup.add(shaftMesh);

  // --- Create Arrow Head ---
  // ConeGeometry: also shift so base is at y=0.
  const headHeight = 0.25;
  const headGeometry = new THREE.ConeGeometry(0.1, headHeight, 12, 1);
  headGeometry.translate(0, headHeight / 2, 0);
  const headMesh = new THREE.Mesh(headGeometry, arrowMaterial);
  // Position the head on top of the shaft
  headMesh.position.y = shaftHeight;
  arrowGroup.add(headMesh);

  return arrowGroup;
}

// Create a group for the arrows that illustrate polarity around the OH group
const arrowsForPolarity = new THREE.Group();
window.arrowsForPolarity = arrowsForPolarity;

// For a clear visual explanation of molecular polarity, create two arrows facing in opposite directions.
// One arrow indicates the direction of the partial negative (oxygen) and one for the partial positive (hydrogen).
const arrowNegative = createArrow();
const arrowPositive = createArrow();

// Orientation adjustments:
// For the arrow associated with the negative partial charge (oxygen), keep it pointing upward.
arrowNegative.rotation.z = 0;  // Default arrow points up (along +Y)

// For the arrow associated with the positive partial charge (hydrogen), flip it so it points down.
arrowPositive.rotation.z = Math.PI; 

// Position the arrows around the polar OH group.
// Based on the ethanol molecule example, the Oxygen is at (3, 0, 0) and the hydrogen attached to it is at around (3.8, 0, 0).
// We'll position the arrows near the OH group.
// For the negative arrow (oxygen), place it slightly offset from the oxygen atom.
arrowNegative.position.set(3, 0.2, 0);

// For the positive arrow (hydrogen), place it near the hydrogen side.
arrowPositive.position.set(3.8, 0.2, 0);

// Add both arrows to the polarity group
arrowsForPolarity.add(arrowNegative);
arrowsForPolarity.add(arrowPositive);

// Add a simple animation behavior to the arrows: a gentle oscillation along their local Y-axis.
// This will simulate "continuous directional movement".
arrowsForPolarity.userData.phase = 0;
arrowsForPolarity.onBeforeRender = function(renderer, scene, camera, geometry, material, group) {
  // Increase phase based on frame delta time (assume approximately 60fps)
  this.userData.phase += 0.02;
  // Compute a small offset oscillation (-0.05 to 0.05)
  const offset = Math.sin(this.userData.phase) * 0.05;
  // Apply the offset along arrow's Y-axis for both arrows
  arrowNegative.position.y = 0.2 + offset;
  arrowPositive.position.y = 0.2 + offset;
};

// Finally, add the arrows group to the scene.
// The group is designed to be associated with the Ethanol_Molecule_Model and Hydroxyl_Glow_Effect.
scene.add(arrowsForPolarity);

// GeometryAgent LLM-generated code
// Create materials for the water molecule
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999 });

// Create a shared sphere geometry for atoms
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create a group to hold the water molecule parts
const waterMolecule = new THREE.Group();
window.waterMolecule = waterMolecule;

// Oxygen atom at the center (V-vertex)
const oxygen = new THREE.Mesh(atomGeometry, oxygenMaterial);
oxygen.position.set(0, 0, 0);
oxygen.scale.set(0.55, 0.55, 0.55);

// Define bond length (scaled for clear interaction with ethanol) and angle for H-O-H
const bondLength = 1; // you can adjust this value as needed
const halfAngle = THREE.MathUtils.degToRad(52.25); // Half of ~104.5° angle

// Hydrogen atoms positioned in a V-shape relative to oxygen
const hydrogen1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen1.position.set(bondLength * Math.sin(halfAngle), bondLength * Math.cos(halfAngle), 0);
hydrogen1.scale.set(0.3, 0.3, 0.3);

const hydrogen2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen2.position.set(-bondLength * Math.sin(halfAngle), bondLength * Math.cos(halfAngle), 0);
hydrogen2.scale.set(0.3, 0.3, 0.3);

// Function to create a cylindrical bond between two atoms
function createBond(start, end) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  
  const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);
  
  // Position bond halfway between start and end
  bond.position.copy(start);
  bond.position.lerp(end, 0.5);
  
  // Align bond with the direction between atoms
  bond.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());
>>>>>>> 30fdb31ebc8832c4e9ed46ccd352c293bdb5e5c3
  
  return bond;
}

<<<<<<< HEAD
// Create bonds between the oxygen and the two hydrogens
const bond1 = createBond(oxygen.position, hydrogen1.position);
const bond2 = createBond(oxygen.position, hydrogen2.position);

// Add bonds to the water molecule group
waterMolecule.add(bond1, bond2);
=======
// Create bonds connecting oxygen to each hydrogen
const bond1 = createBond(oxygen.position, hydrogen1.position);
const bond2 = createBond(oxygen.position, hydrogen2.position);

// Assemble the water molecule group
waterMolecule.add(oxygen, hydrogen1, hydrogen2, bond1, bond2);
>>>>>>> 30fdb31ebc8832c4e9ed46ccd352c293bdb5e5c3

// Add the water molecule group to the scene
scene.add(waterMolecule);

// GeometryAgent LLM-generated code
<<<<<<< HEAD
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
=======
// Create a group to hold all dashed hydrogen bond lines
const dashedHBGroup = new THREE.Group();
window.dashedHBGroup = dashedHBGroup;

// Define the dashed line material
const dashedLineMaterial = new THREE.LineDashedMaterial({
  color: 0xadd8e6,    // light blue
  dashSize: 0.3,      // length of a dash
  gapSize: 0.2,       // length of the gap between dashes
  linewidth: 1        // thin line (note: linewidth is mostly effective on WebGLRenderer with special settings)
});

// Example 1: Dashed line from the ethanol hydroxyl group (assumed at (3,0,0))
// to a water molecule interaction point (example position: (1.5, 2, 0))
{
  const points = [];
  points.push(new THREE.Vector3(3, 0, 0));
  points.push(new THREE.Vector3(1.5, 2, 0));
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  // computeLineDistances is required for dashed lines to work properly
  geometry.computeLineDistances();
  const dashedLine1 = new THREE.Line(geometry, dashedLineMaterial.clone());
  dashedHBGroup.add(dashedLine1);
}

// Example 2: Dashed line from the ethanol hydroxyl group (assumed at (3,0,0))
// to another water molecule interaction point (example position: (4.5, 2, 0))
{
  const points = [];
  points.push(new THREE.Vector3(3, 0, 0));
  points.push(new THREE.Vector3(4.5, 2, 0));
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  geometry.computeLineDistances();
  const dashedLine2 = new THREE.Line(geometry, dashedLineMaterial.clone());
  dashedHBGroup.add(dashedLine2);
}

// Function to animate the dashed lines for a pulsing or intermittent appearance
window.animateDashedHydrogenBondLines = function() {
  // Loop through each line in the group and update its dash offset to create an animated effect.
  dashedHBGroup.children.forEach(line => {
    if (line.material && typeof line.material.dashOffset !== 'undefined') {
      line.material.dashOffset -= 0.01; // adjust for desired animation speed
    }
  });
};

// Add the dashed hydrogen bond lines group to the scene
scene.add(dashedHBGroup);

// Optionally, if you have a main animation loop, call window.animateDashedHydrogenBondLines() each frame.
// Example:
// function animate() {
//   requestAnimationFrame(animate);
//   window.animateDashedHydrogenBondLines();
//   renderer.render(scene, camera);
// }
// animate();

// GeometryAgent LLM-generated code
// Create a gradient shader material for the background
const vertexShader = `
  varying vec2 vUv;
  void main() {
    vUv = uv;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  }
`;

const fragmentShader = `
  uniform vec3 topColor;
  uniform vec3 bottomColor;
  uniform float opacity;
  varying vec2 vUv;
  void main() {
    // Mix bottomColor (blue) at vUv.y == 0 and topColor (purple) at vUv.y == 1
    gl_FragColor = vec4(mix(bottomColor, topColor, vUv.y), opacity);
  }
`;

// Define uniforms for our shader material
// Using blue (#0000ff) at the bottom and purple (#800080) at the top
const uniforms = {
  topColor: { value: new THREE.Color(0x800080) },
  bottomColor: { value: new THREE.Color(0x0000ff) },
  // Start with 0 opacity to later allow a gradual fade in (animation/tween can modify this uniform)
  opacity: { value: 0.0 }
};

const gradientMaterial = new THREE.ShaderMaterial({
  uniforms: uniforms,
  vertexShader: vertexShader,
  fragmentShader: fragmentShader,
  transparent: true
});

// Create a large plane that covers the background.
// The dimensions are chosen to cover the view regardless of camera settings.
const planeGeometry = new THREE.PlaneGeometry(200, 200);
const gradientMesh = new THREE.Mesh(planeGeometry, gradientMaterial);

// Position the plane so that it's behind the main scene elements.
// Adjust the z position according to your scene's requirements.
gradientMesh.position.set(0, 0, -50);

// Optionally, you can rotate the plane if needed. For a straightforward background, no rotation is needed.
window.Gradient_Background = gradientMesh;

// Add the gradient background mesh to the scene
scene.add(gradientMesh);

// To achieve the "gradual fade in" effect, animate the opacity uniform somewhere in your main render loop or using a tweening library.
// For example, using GSAP (if available), you might do:
// gsap.to(gradientMaterial.uniforms.opacity, { value: 1.0, duration: 3, delay: 1.15 });

// GeometryAgent LLM-generated code
// Create a canvas for the text overlay
const canvas = document.createElement('canvas');
canvas.width = 512;
canvas.height = 256;
const ctx = canvas.getContext('2d');

// Draw a semi-transparent dark background
ctx.fillStyle = "rgba(0, 0, 0, 0.5)";
ctx.fillRect(0, 0, canvas.width, canvas.height);

// Set up clear, legible white text
ctx.fillStyle = "white";
ctx.font = "30px Arial"; // Medium size and clear font
ctx.textAlign = "center";
ctx.textBaseline = "middle";

// Write the summary text (you can customize the content as needed)
const overlayText = "Ethanol: Structure, Polarity, & Interactivity";
ctx.fillText(overlayText, canvas.width / 2, canvas.height / 2);

// Create a texture from the canvas and update it
const texture = new THREE.CanvasTexture(canvas);

// Create a plane geometry with an appropriate aspect ratio
const geometry = new THREE.PlaneGeometry(4, 2);

// Create a material using the canvas texture with transparency enabled
const material = new THREE.MeshBasicMaterial({
  map: texture,
  transparent: true,
  opacity: 0.0, // Start invisible for fade-in effect
  depthTest: false
});

// Create the mesh and store it in a global reference for later use (if needed)
const summaryOverlay = new THREE.Mesh(geometry, material);
window.summaryOverlay = summaryOverlay;

// Position the overlay in the scene (adjust position to suit your layout)
summaryOverlay.position.set(0, 0, 0); // This position can be modified to overlay the scene as desired

// Add the overlay to the global scene
scene.add(summaryOverlay);

// Optional: Smooth fade-in animation (over 2 seconds)
const fadeDuration = 2000;
const fadeStart = performance.now();
function animateFade() {
  const now = performance.now();
  const elapsed = now - fadeStart;
  const t = Math.min(elapsed / fadeDuration, 1);
  material.opacity = t;
  if (t < 1) {
    requestAnimationFrame(animateFade);
  }
}
animateFade();
>>>>>>> 30fdb31ebc8832c4e9ed46ccd352c293bdb5e5c3
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
<<<<<<< HEAD
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
=======
// Get elapsed time in seconds
const elapsedTime = clock.getElapsedTime();

// Use deltaTime for smooth animations (assumed to be defined elsewhere)
// ------------------------------
// At 00:00 – Scene Opening: Fade in dark starry background and zoom in camera
if (elapsedTime < 20) {
  // Smooth camera zoom toward the ethanol molecule
  if (window.camera && window.Ethanol_Molecule_Model) {
    // Lerp camera z-position from 100 (start) to 50 (target) over 20 seconds
    window.camera.position.z = THREE.MathUtils.lerp(100, 50, elapsedTime / 20);
  }
  // Fade in Dark Starry Background (assumes material supports opacity)
  if (window.Dark_Starry_Background && window.Dark_Starry_Background.material) {
    window.Dark_Starry_Background.material.opacity = THREE.MathUtils.lerp(0, 1, elapsedTime / 20);
  }
  // (Voiceover introduction would play here externally)
}

// ------------------------------
// At 00:20 – Rotating view with labeled atoms
if (elapsedTime >= 20 && elapsedTime < 40) {
  const t = (elapsedTime - 20) / 20; // progress from 0 to 1 over 20s window

  // Rotate the ethanol molecule for a revealing view
  if (window.Ethanol_Molecule_Model) {
    window.Ethanol_Molecule_Model.rotation.y += 0.5 * deltaTime;
  }
  // Fade in the Atom Labels (assumes material supports opacity)
  if (window.Atom_Labels && window.Atom_Labels.material) {
    window.Atom_Labels.material.opacity = THREE.MathUtils.lerp(0, 1, t);
  }
  // Optionally, the camera could slowly orbit around the molecule if needed
  if (window.camera && window.Ethanol_Molecule_Model) {
    // Rotate camera around the molecule's Y axis (simple circular orbit)
    const radius = 50;
    const angle = t * Math.PI * 2; // complete one orbit? adjust as desired
    window.camera.position.x = window.Ethanol_Molecule_Model.position.x + radius * Math.sin(angle);
    window.camera.position.z = window.Ethanol_Molecule_Model.position.z + radius * Math.cos(angle);
    window.camera.lookAt(window.Ethanol_Molecule_Model.position);
  }
}

// ------------------------------
// At 00:40 – Detailed bonds and angles
if (elapsedTime >= 40 && elapsedTime < 60) {
  const t = (elapsedTime - 40) / 20; // progress from 0 to 1 over 20s window

  // Slightly adjust camera to focus on the bond between second carbon and oxygen
  if (window.camera && window.Ethanol_Molecule_Model && window.Bond_Angle_Overlay) {
    // Offset camera towards the region of interest (example offsets)
    window.camera.position.x = THREE.MathUtils.lerp(0, 10, t);
    window.camera.position.y = THREE.MathUtils.lerp(0, 5, t);
    window.camera.lookAt(window.Ethanol_Molecule_Model.position);
  }
  // Fade in Bond Angle Overlay that illustrates bond angles
  if (window.Bond_Angle_Overlay && window.Bond_Angle_Overlay.material) {
    window.Bond_Angle_Overlay.material.opacity = THREE.MathUtils.lerp(0, 1, t);
  }
}

// ------------------------------
// At 01:00 – Highlighting the polar OH group
if (elapsedTime >= 60 && elapsedTime < 80) {
  const t = (elapsedTime - 60) / 20; // progress from 0 to 1 over 20s window

  // Fade in and intensify the glow effect around the hydroxyl group
  if (window.Hydroxyl_Glow_Effect && window.Hydroxyl_Glow_Effect.material) {
    window.Hydroxyl_Glow_Effect.material.opacity = THREE.MathUtils.lerp(0, 1, t);
    // Optional pulsing scale effect for emphasis
    const scale = 1 + 0.1 * Math.sin(elapsedTime * 5);
    window.Hydroxyl_Glow_Effect.scale.set(scale, scale, scale);
  }
  // Fade in polarity arrows that point out partial charges
  if (window.Arrows_For_Polarity && window.Arrows_For_Polarity.material) {
    window.Arrows_For_Polarity.material.opacity = THREE.MathUtils.lerp(0, 1, t);
  }
}

// ------------------------------
// At 01:20 – Ethanol in hydrogen bonding
if (elapsedTime >= 80 && elapsedTime < 100) {
  const t = (elapsedTime - 80) / 20; // progress from 0 to 1 over 20s window

  // Simulate water molecule dynamics (e.g., floating around with a gentle oscillation)
  if (window.Water_Molecule_Model) {
    window.Water_Molecule_Model.position.x = 5 * Math.sin(elapsedTime);
    window.Water_Molecule_Model.position.y = 5 * Math.cos(elapsedTime);
  }
  // Fade in dashed hydrogen bond lines to illustrate the hydrogen bonding network
  if (window.Dashed_Hydrogen_Bond_Lines && window.Dashed_Hydrogen_Bond_Lines.material) {
    window.Dashed_Hydrogen_Bond_Lines.material.opacity = THREE.MathUtils.lerp(0, 1, t);
  }
  // Continue slow rotation of the ethanol molecule to showcase interactions
  if (window.Ethanol_Molecule_Model) {
    window.Ethanol_Molecule_Model.rotation.y += 0.2 * deltaTime;
  }
}

// ------------------------------
// At 01:40 – Summary of Ethanol Structure and final zoom out
if (elapsedTime >= 100) {
  // For a smooth transition, we map t from 0 to 1 over the first 20 seconds after 01:40
  const t = Math.min((elapsedTime - 100) / 20, 1);

  // Zoom out the camera to display the broader context of chemical interactions
  if (window.camera && window.Ethanol_Molecule_Model) {
    window.camera.position.z = THREE.MathUtils.lerp(50, 80, t);
    window.camera.lookAt(window.Ethanol_Molecule_Model.position);
  }
  // Fade in gradient background replacing the starry scene
  if (window.Gradient_Background && window.Gradient_Background.material) {
    window.Gradient_Background.material.opacity = THREE.MathUtils.lerp(0, 1, t);
  }
  // Fade in a summary overlay that recaps the ethanol molecule’s properties
  if (window.Summary_Overlay && window.Summary_Overlay.material) {
    window.Summary_Overlay.material.opacity = THREE.MathUtils.lerp(0, 1, t);
  }
  // Maintain a gentle rotation of the ethanol molecule in the final view
  if (window.Ethanol_Molecule_Model) {
    window.Ethanol_Molecule_Model.rotation.y += 0.2 * deltaTime;
  }
}
>>>>>>> 30fdb31ebc8832c4e9ed46ccd352c293bdb5e5c3
    
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
