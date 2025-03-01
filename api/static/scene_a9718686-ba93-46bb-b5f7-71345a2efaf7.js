// Scientific visualization: Exploring sp2 Hybridization in Alkenes
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
// Materials
const carbonMaterial = new THREE.MeshStandardMaterial({ 
  color: 0x333333, 
  metalness: 0.8, 
  roughness: 0.2 
});
const hydrogenMaterial = new THREE.MeshStandardMaterial({ 
  color: 0xffffff, 
  metalness: 0.2, 
  roughness: 0.8 
});
const bondMaterial = new THREE.MeshStandardMaterial({ 
  color: 0x999999, 
  metalness: 0.3, 
  roughness: 0.7 
});

// Orbital materials – sp2 orbitals will have a glowing blue tint and semi-transparency,
// and the unhybridized p orbital uses a similar look.
const sp2OrbitalMaterial = new THREE.MeshPhongMaterial({ 
  color: 0x00aaff, 
  transparent: true, 
  opacity: 0.5, 
  emissive: 0x00aaff 
});
const pOrbitalMaterial = new THREE.MeshPhongMaterial({ 
  color: 0x00aaff, 
  transparent: true, 
  opacity: 0.5, 
  emissive: 0x00aaff 
});

// Material for the diffused electron cloud (soft yellow glow)
const electronCloudMaterial = new THREE.MeshPhongMaterial({ 
  color: 0xffff66, 
  transparent: true, 
  opacity: 0.3, 
  emissive: 0xffff66 
});

// Atom Geometries
const carbonGeometry = new THREE.SphereGeometry(0.4, 32, 32);
const hydrogenGeometry = new THREE.SphereGeometry(0.2, 32, 32);

// Create a group for the alkene molecule model
const alkeneMolecule = new THREE.Group();
window.alkeneMolecule = alkeneMolecule;

// Create Carbon atoms (representing the sp2-hybridized centers)
// Position them along the x-axis to represent the C=C double bond.
const carbon1 = new THREE.Mesh(carbonGeometry, carbonMaterial);
const carbon2 = new THREE.Mesh(carbonGeometry, carbonMaterial);
carbon1.position.set(-0.75, 0, 0);
carbon2.position.set(0.75, 0, 0);
alkeneMolecule.add(carbon1, carbon2);

// Create Hydrogen atoms (for an ethene-like molecule, C2H4)
// Carbon1 hydrogens
const hydrogen1 = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);
const hydrogen2 = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);
hydrogen1.position.set(-1.5, 0.9, 0);
hydrogen2.position.set(-1.5, -0.9, 0);
// Carbon2 hydrogens
const hydrogen3 = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);
const hydrogen4 = new THREE.Mesh(hydrogenGeometry, hydrogenMaterial);
hydrogen3.position.set(1.5, 0.9, 0);
hydrogen4.position.set(1.5, -0.9, 0);
alkeneMolecule.add(hydrogen1, hydrogen2, hydrogen3, hydrogen4);

// Function to create a bond between two points using a CylinderGeometry.
function createBond(start, end, radius = 0.08) {
  const direction = new THREE.Vector3().subVectors(end, start);
  const length = direction.length();
  const bondGeometry = new THREE.CylinderGeometry(radius, radius, length, 16);
  const bond = new THREE.Mesh(bondGeometry, bondMaterial);
  
  // Position the bond midway between the two atoms.
  const midPoint = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  bond.position.copy(midPoint);
  
  // Align the bond along the vector connecting the two atoms.
  bond.quaternion.setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    direction.clone().normalize()
  );
  
  return bond;
}

// Create bonds
// Sigma bond between the two carbons
const sigmaBond = createBond(carbon1.position, carbon2.position);

// Pi bond – representing the second bond of the C=C double bond.
// Slightly offset in the z-direction so it doesn't completely overlap with the sigma bond.
const offset = new THREE.Vector3(0, 0, 0.15);
const carbon1PosPi = carbon1.position.clone().add(offset);
const carbon2PosPi = carbon2.position.clone().add(offset);
const piBond = createBond(carbon1PosPi, carbon2PosPi, 0.05);

// Bonds from carbons to hydrogens
const bondC1H1 = createBond(carbon1.position, hydrogen1.position);
const bondC1H2 = createBond(carbon1.position, hydrogen2.position);
const bondC2H1 = createBond(carbon2.position, hydrogen3.position);
const bondC2H2 = createBond(carbon2.position, hydrogen4.position);

// Add bonds to the molecule group
alkeneMolecule.add(sigmaBond, piBond, bondC1H1, bondC1H2, bondC2H1, bondC2H2);

// Function to add orbital visualizations onto a given carbon atom.
// This adds three sp2 orbitals arranged in a trigonal planar layout (lying in the plane)
// and one unhybridized p orbital oriented perpendicular to that plane.
// It also adds a dashed circle to mark the orbital plane and diffused electron cloud glows.
function addOrbitals(atom, orbitalRadius = 0.6) {
  // Create three sp2 orbitals in the xy plane.
  const anglesDeg = [0, 120, 240];
  anglesDeg.forEach(angleDeg => {
    const angleRad = THREE.MathUtils.degToRad(angleDeg);
    const orbitalGeometry = new THREE.CircleGeometry(orbitalRadius, 32);
    const orbital = new THREE.Mesh(orbitalGeometry, sp2OrbitalMaterial);
    // Position the orbital very close to the atom’s center with a slight radial offset.
    orbital.position.set(
      Math.cos(angleRad) * 0.05,
      Math.sin(angleRad) * 0.05,
      0
    );
    // The CircleGeometry is created in the xy plane already.
    atom.add(orbital);
  });
  
  // Create the unhybridized p orbital.
  // This orbital is depicted as a circular disk rotated 90° so that it is perpendicular to the sp2 plane.
  const pOrbitalGeometry = new THREE.CircleGeometry(orbitalRadius * 0.8, 32);
  const pOrbital = new THREE.Mesh(pOrbitalGeometry, pOrbitalMaterial);
  pOrbital.rotation.x = Math.PI / 2; // make it vertical (perpendicular to the xy plane)
  pOrbital.position.set(0, 0, 0.05);
  atom.add(pOrbital);
  
  // Annotate the sp2 orbital plane with a dashed circle.
  const circleGeometry = new THREE.CircleGeometry(orbitalRadius, 64);
  // Remove the center vertex to draw only the outline.
  circleGeometry.vertices.shift();
  const dashMaterial = new THREE.LineDashedMaterial({ color: 0x000000, dashSize: 0.1, gapSize: 0.05 });
  const orbitalPlane = new THREE.LineLoop(circleGeometry, dashMaterial);
  orbitalPlane.computeLineDistances();
  // Rotate the dashed circle so that it lies in the xy plane.
  orbitalPlane.rotation.x = Math.PI / 2;
  atom.add(orbitalPlane);
  
  // Add diffused electron cloud glows – one above and one below the sp2 plane.
  const cloudGeometry = new THREE.CircleGeometry(orbitalRadius * 1.2, 32);
  const cloudTop = new THREE.Mesh(cloudGeometry, electronCloudMaterial);
  cloudTop.position.set(0, 0, 0.12);
  atom.add(cloudTop);
  
  const cloudBottom = new THREE.Mesh(cloudGeometry, electronCloudMaterial);
  cloudBottom.position.set(0, 0, -0.12);
  atom.add(cloudBottom);
}

// Add orbitals to each carbon atom
addOrbitals(carbon1);
addOrbitals(carbon2);

// Optional: set a userData property for rotation animation (to be used in your render loop)
alkeneMolecule.userData.rotationSpeed = 0.001; // slow, continuous rotation

// Add the complete alkene molecule group to the global scene.
scene.add(alkeneMolecule);

// GeometryAgent LLM-generated code
// Create a main group for the hybridization comparative diagram
const hybridDiagram = new THREE.Group();
window.HybridDiagram = hybridDiagram;

// Create a background panel (15 x 10 units) with a neutral gray color
const bgGeometry = new THREE.PlaneGeometry(15, 10);
const bgMaterial = new THREE.MeshBasicMaterial({ color: 0x808080 });
const bgMesh = new THREE.Mesh(bgGeometry, bgMaterial);
bgMesh.position.set(0, 0, -0.1); // slightly behind to avoid z-fighting with overlays
hybridDiagram.add(bgMesh);

// Utility function to create a text sprite for annotations
function createTextSprite(message, parameters = {}) {
    const fontface = parameters.hasOwnProperty("fontface") ? parameters["fontface"] : "Arial";
    const fontsize = parameters.hasOwnProperty("fontsize") ? parameters["fontsize"] : 48;
    const borderThickness = parameters.hasOwnProperty("borderThickness") ? parameters["borderThickness"] : 4;
    const textColor = parameters.hasOwnProperty("textColor") ? parameters["textColor"] : { r:255, g:255, b:255 };
    
    const canvas = document.createElement('canvas');
    const context = canvas.getContext('2d');
    context.font = fontsize + "px " + fontface;
    
    // get size data (approx.)
    const metrics = context.measureText(message);
    const textWidth = metrics.width;
    
    // adjust canvas size
    canvas.width = textWidth + borderThickness * 2;
    canvas.height = fontsize + borderThickness * 2;
    // re-set font after resizing canvas
    context.font = fontsize + "px " + fontface;
    
    // background (optional transparent background)
    context.fillStyle = "rgba(0, 0, 0, 0)"; // transparent
    context.fillRect(0, 0, canvas.width, canvas.height);
    
    // text color
    context.fillStyle = "rgba(" + textColor.r + "," + textColor.g + "," + textColor.b + ", 1.0)";
    context.fillText(message, borderThickness, fontsize + borderThickness - 10);
    
    // use canvas as texture
    const texture = new THREE.CanvasTexture(canvas);
    texture.needsUpdate = true;
    const spriteMaterial = new THREE.SpriteMaterial({ map: texture, transparent: true });
    const sprite = new THREE.Sprite(spriteMaterial);
    
    // scale sprite to make text clearly visible
    sprite.scale.set(3, 1.5, 1);
    return sprite;
}

// Create a group for the sp2 diagram (left panel)
// This diagram represents sp2 hybridized carbon with a double bond (vivid blue accents)
const sp2Group = new THREE.Group();
sp2Group.position.set(-3.75, 0, 0); // left of center

// Create an energy profile line for sp2 using a series of points
const sp2Points = [];
sp2Points.push(new THREE.Vector3(-3, -2, 0));
sp2Points.push(new THREE.Vector3(-2, 1, 0));
sp2Points.push(new THREE.Vector3(-1, 2, 0));
sp2Points.push(new THREE.Vector3(0, 1.5, 0));
sp2Points.push(new THREE.Vector3(1, 2.5, 0));
sp2Points.push(new THREE.Vector3(2, 0, 0));
const sp2Geometry = new THREE.BufferGeometry().setFromPoints(sp2Points);
const sp2Material = new THREE.LineBasicMaterial({ color: 0x00aaff, linewidth: 3 });
const sp2Line = new THREE.Line(sp2Geometry, sp2Material);
sp2Group.add(sp2Line);

// Add annotated bullet points for sp2 diagram
const sp2Annotation1 = createTextSprite("• Reactivity: Higher", { fontsize: 32, textColor: { r:255, g:255, b:255 } });
sp2Annotation1.position.set(-1.5, 2.8, 0);
sp2Annotation1.material.opacity = 0; // start invisible for fade-in animation

const sp2Annotation2 = createTextSprite("• Bond Energy: Lower", { fontsize: 32, textColor: { r:255, g:255, b:255 } });
sp2Annotation2.position.set(-1.5, 2.0, 0);
sp2Annotation2.material.opacity = 0;

const sp2Annotation3 = createTextSprite("• Orbital: Planar", { fontsize: 32, textColor: { r:255, g:255, b:255 } });
sp2Annotation3.position.set(-1.5, 1.2, 0);
sp2Annotation3.material.opacity = 0;

sp2Group.add(sp2Annotation1, sp2Annotation2, sp2Annotation3);

// Create a group for the sp3 diagram (right panel)
// This diagram represents sp3 hybridized carbon with single bonds (warm orange accents)
const sp3Group = new THREE.Group();
sp3Group.position.set(3.75, 0, 0); // right of center

// Create an energy profile line for sp3 using a series of points
const sp3Points = [];
sp3Points.push(new THREE.Vector3(-3, -1, 0));
sp3Points.push(new THREE.Vector3(-2, 0.5, 0));
sp3Points.push(new THREE.Vector3(-1, 1, 0));
sp3Points.push(new THREE.Vector3(0, 0.5, 0));
sp3Points.push(new THREE.Vector3(1, 1.2, 0));
sp3Points.push(new THREE.Vector3(2, -0.5, 0));
const sp3Geometry = new THREE.BufferGeometry().setFromPoints(sp3Points);
const sp3Material = new THREE.LineBasicMaterial({ color: 0xff6600, linewidth: 3 });
const sp3Line = new THREE.Line(sp3Geometry, sp3Material);
sp3Group.add(sp3Line);

// Add annotated bullet points for sp3 diagram
const sp3Annotation1 = createTextSprite("• Reactivity: Lower", { fontsize: 32, textColor: { r:255, g:255, b:255 } });
sp3Annotation1.position.set(-1.5, 1.8, 0);
sp3Annotation1.material.opacity = 0;

const sp3Annotation2 = createTextSprite("• Bond Energy: Higher", { fontsize: 32, textColor: { r:255, g:255, b:255 } });
sp3Annotation2.position.set(-1.5, 1.0, 0);
sp3Annotation2.material.opacity = 0;

const sp3Annotation3 = createTextSprite("• Orbital: Tetrahedral", { fontsize: 32, textColor: { r:255, g:255, b:255 } });
sp3Annotation3.position.set(-1.5, 0.2, 0);
sp3Annotation3.material.opacity = 0;

sp3Group.add(sp3Annotation1, sp3Annotation2, sp3Annotation3);

// Add the two panels (sp2 and sp3 groups) to the main hybrid diagram group
hybridDiagram.add(sp2Group);
hybridDiagram.add(sp3Group);

// Set initial positions for sliding-in animation (offscreen horizontally)
// These positions can later be animated using external tweening libraries (e.g., GSAP) or custom animation loops.
sp2Group.position.x -= 10; // start offscreen left
sp3Group.position.x += 10; // start offscreen right

// Similarly, annotations are initially faded-out (opacity already set to 0) and can be animated to fade in

// Finally, add the complete hybrid diagram group to the scene
scene.add(hybridDiagram);

// Optional: Store global references to groups if needed for animation later
window.sp2Diagram = sp2Group;
window.sp3Diagram = sp3Group;
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
const elapsedTime = clock.getElapsedTime();
const deltaTime = clock.getDelta();

// ────────────────────────────────────────────────────────────
// At 00:00 - Introduction to sp2 Hybridization
// Slowly rotate the alkene molecule and fade in the introductory label
if (elapsedTime < 20) {
  if (window.Alkene_Molecule_Model) {
    window.Alkene_Molecule_Model.rotation.y += 0.1 * deltaTime;
  }
  // Assume an introductory label exists as window.introLabel
  if (window.introLabel) {
    // Fade in over the first 5 seconds
    window.introLabel.material.opacity = Math.min(1, elapsedTime / 5);
  }
}

// ────────────────────────────────────────────────────────────
// At 00:20 - Visualizing Carbon's sp2 Orbitals
// Zoom in the camera to focus on one carbon atom and highlight its orbitals
if (elapsedTime >= 20 && elapsedTime < 40) {
  // Compute normalized time progress (0 to 1) between 20s and 40s
  const t = (elapsedTime - 20) / 20;
  if (window.camera) {
    // Lerp the camera position toward a close-up view (example positions)
    // (Assumes initial position ~ (0, 0, 5) and target (0, 1, 2))
    window.camera.position.lerp(new THREE.Vector3(0, 1, 2), 0.05);
    window.camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  
  // Highlight the three sp2 orbitals and the unhybridized p orbital.
  // (Assumes each orbital is a child mesh named accordingly.)
  ['sp2_Orbital1', 'sp2_Orbital2', 'sp2_Orbital3', 'p_Orbital'].forEach(name => {
    const orbital = window.Alkene_Molecule_Model.getObjectByName(name);
    if (orbital) {
      // Gradually increase the scale to emphasize the orbital (e.g., up to 1.5×)
      orbital.scale.setScalar(1 + 0.5 * t);
    }
  });
  
  // Optionally, animate directional arrows if they exist (e.g., window.orbitalArrows)
  if (window.orbitalArrows) {
    window.orbitalArrows.children.forEach(arrow => {
      arrow.material.opacity = Math.min(1, t);
    });
  }
}

// ────────────────────────────────────────────────────────────
// At 00:40 - Sigma & Pi Bond Formation
// Animate bonds forming: extend sigma bonds and form the sideways-overlapping pi bond.
if (elapsedTime >= 40 && elapsedTime < 60) {
  const t = (elapsedTime - 40) / 20; // progresses from 0 to 1 over 20s
  if (window.Alkene_Molecule_Model) {
    // Animate sigma bonds extension (assuming objects named 'sigmaBond1', 'sigmaBond2', etc.)
    ['sigmaBond1', 'sigmaBond2'].forEach(name => {
      const sigmaBond = window.Alkene_Molecule_Model.getObjectByName(name);
      if (sigmaBond) {
        // Extend along the local z-axis from scale 1 to 2
        sigmaBond.scale.z = 1 + t;
      }
    });
    // Animate the pi bond formation (assuming an object named 'piBond')
    const piBond = window.Alkene_Molecule_Model.getObjectByName('piBond');
    if (piBond) {
      // Increase the lateral scale and slightly adjust opacity to emphasize overlap
      const scaleFactor = 1 + 0.3 * t;
      piBond.scale.set(scaleFactor, scaleFactor, 1);
      if(piBond.material) {
        piBond.material.opacity = Math.min(1, t);
      }
    }
  }
}

// ────────────────────────────────────────────────────────────
// At 01:00 - Planar Geometry & Electron Cloud
// Rearrange the camera to an overhead view and reveal the planar dotted lines and diffused electron cloud glow.
if (elapsedTime >= 60 && elapsedTime < 80) {
  const t = (elapsedTime - 60) / 20;
  if (window.camera) {
    // Lerp camera to an overhead view (example: from its current position to (0, 5, 0.1))
    window.camera.position.lerp(new THREE.Vector3(0, 5, 0.1), 0.05);
    window.camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  // Animate dotted plane to appear (assumes an object named 'dottedPlane' is a child of the model)
  const dottedPlane = window.Alkene_Molecule_Model.getObjectByName('dottedPlane');
  if (dottedPlane && dottedPlane.material) {
    dottedPlane.material.opacity = Math.min(1, t);
  }
  // Animate diffused glow for the pi electron cloud (assumes an object named 'piElectronCloud')
  const piElectronCloud = window.Alkene_Molecule_Model.getObjectByName('piElectronCloud');
  if (piElectronCloud && piElectronCloud.material) {
    piElectronCloud.material.opacity = Math.min(1, t);
    // Add a subtle pulsating effect
    const pulse = 1 + 0.2 * Math.sin(elapsedTime * 3);
    piElectronCloud.scale.set(pulse, pulse, pulse);
  }
}

// ────────────────────────────────────────────────────────────
// At 01:20 - Reactivity of the Pi Bond
// Focus on the electron density with a pulsating, glowing pi bond region to indicate its reactivity.
if (elapsedTime >= 80 && elapsedTime < 100) {
  const t = (elapsedTime - 80) / 20;
  const piElectronCloud = window.Alkene_Molecule_Model.getObjectByName('piElectronCloud');
  if (piElectronCloud) {
    // Increase pulsation and brightness to emphasize reactivity
    const pulse = 1 + 0.3 * Math.sin((elapsedTime - 80) * 5);
    piElectronCloud.scale.set(pulse, pulse, pulse);
    if (piElectronCloud.material && piElectronCloud.material.emissive) {
      piElectronCloud.material.emissiveIntensity = 0.5 + 0.5 * t;
    }
  }
  // Optionally, fade in a reactivity annotation label (assume window.piBondLabel)
  if (window.piBondLabel) {
    window.piBondLabel.material.opacity = Math.min(1, (elapsedTime - 80) / 5);
  }
}

// ────────────────────────────────────────────────────────────
// At 01:40 - Comparing sp2 & sp3 Structures
// Slide in two energy diagrams side-by-side to compare the hybridization states.
// (Assumes existence of objects: window.energyDiagram_sp2 and window.energyDiagram_sp3)
if (elapsedTime >= 100 && elapsedTime < 120) {
  const t = (elapsedTime - 100) / 20;
  if (window.energyDiagram_sp2) {
    // Slide from left: starting at x = -5, moving to x = 0
    window.energyDiagram_sp2.position.x = -5 + 5 * t;
    if (window.energyDiagram_sp2.material) {
      window.energyDiagram_sp2.material.opacity = t;
    }
  }
  if (window.energyDiagram_sp3) {
    // Slide from right: starting at x = 5, moving to x = 0
    window.energyDiagram_sp3.position.x = 5 - 5 * t;
    if (window.energyDiagram_sp3.material) {
      window.energyDiagram_sp3.material.opacity = t;
    }
  }
  // Optionally, fade out details of the alkene model to emphasize the diagrams
  if (window.Alkene_Molecule_Model) {
    window.Alkene_Molecule_Model.traverse(child => {
      if (child.material && child.name !== 'energyDiagram_sp2' && child.name !== 'energyDiagram_sp3') {
        child.material.opacity = Math.max(0, 1 - t);
      }
    });
  }
}

// ────────────────────────────────────────────────────────────
// At 02:00 - Summary of sp2 Hybridization
// Pull the camera back for an overall summary and display bullet points.
if (elapsedTime >= 120 && elapsedTime < 140) {
  const t = (elapsedTime - 120) / 20;
  if (window.camera) {
    // Lerp the camera to a wide view (example position: (0, 0, 10))
    window.camera.position.lerp(new THREE.Vector3(0, 0, 10), 0.05);
    window.camera.lookAt(new THREE.Vector3(0, 0, 0));
  }
  // Fade in summary bullet points (assumes an object window.summaryText)
  if (window.summaryText && window.summaryText.material) {
    window.summaryText.material.opacity = t;
  }
}

// ────────────────────────────────────────────────────────────
// At 03:00 - Display Hybridization Comparative Diagram
// Bring in the comparative diagram that contrasts sp2 and sp3 hybridization.
if (elapsedTime >= 180) {
  // Animate diagram sliding in from below (y from -5 to 0)
  const t = Math.min(1, (elapsedTime - 180) / 5);
  if (window.Hybridization_Comparative_Diagram) {
    window.Hybridization_Comparative_Diagram.position.y = -5 + 5 * t;
    if (window.Hybridization_Comparative_Diagram.material) {
      window.Hybridization_Comparative_Diagram.material.opacity = t;
    }
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
