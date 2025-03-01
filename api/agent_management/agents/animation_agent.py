"""
Animation Agent - Uses LLM to generate Three.js animation code based on user prompts.
Specialized for scientific visualizations including molecular structures.
"""

from llm_service import LLMService

class AnimationAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service

    def get_animation_snippet(self, user_prompt: str) -> str:
        """
        Generate Three.js animation code using LLM based on user prompt.
        This code will animate the objects created by the GeometryAgent,
        with special focus on scientific visualizations like molecular structures.
        """
        prompt_for_llm = f"You are an expert Three.js animator specializing in scientific visualizations, particularly molecular structures and chemical reactions. The user wants to animate a scene with the following prompt: '{user_prompt}'\n\n" + r"""
Requirements:
- Generate JavaScript animation code for objects created by the GeometryAgent.
- Objects are available via global window references (e.g., window.molecule, window.atom1).
- Your code will be injected into an existing Three.js scene with a standard animation loop.
- Create scientifically accurate yet visually appealing animations that represent the prompt.
- Implement a updateScene(deltaTime) function to handle animation updates.
- Use deltaTime parameter for smooth, frame-rate independent animations.
- For chemical/molecular visualizations, consider appropriate molecular behaviors.

# Guide to Animating Molecular Structures

When animating molecular structures in Three.js:
1. Molecules are typically organized as Groups containing atoms (spheres) and bonds (cylinders).
2. Identify the type of molecular motion appropriate to the visualization:
   - Rotation (molecules spinning)
   - Vibration (atoms oscillating around equilibrium positions)
   - Translation (molecules moving through space)
   - Reaction transitions (bond breaking/forming, atom rearrangement)
3. For atom vibrations, modulate position with sine/cosine waves to create natural motion.
4. For reactions, consider gradual transformations between states.
5. Use camera movement strategically to highlight important aspects.

## Example: Ethanol Molecule Animation

For a prompt like "animate an ethanol molecule showing molecular vibrations":

```javascript
// AnimationAgent code for ethanol molecule visualization
let elapsedTime = 0;
const rotationSpeed = 0.15;

// Vibration parameters
const vibrationFrequency = 3.0; // Higher = faster vibration
const carbonVibrationAmount = 0.01; // Less movement for carbon (heavier)
const oxygenVibrationAmount = 0.015; // Medium movement for oxygen
const hydrogenVibrationAmount = 0.025; // More movement for hydrogen (lighter)

// Stores original positions of atoms
const originalPositions = new Map();

// Initialize original positions on first frame
let initialized = false;

function initializePositions() {
    if (typeof window.molecule === 'undefined') return;
    
    window.molecule.children.forEach((child) => {
        // Store original positions to vibrate around
        originalPositions.set(child.uuid, child.position.clone());
        
        // Determine atom type based on scale (from GeometryAgent)
        if (child.scale.x === 0.5) {
            child.userData.atomType = 'carbon';
        } else if (child.scale.x === 0.55) {
            child.userData.atomType = 'oxygen';
        } else if (child.scale.x === 0.3) {
            child.userData.atomType = 'hydrogen';
        } else {
            child.userData.atomType = 'bond';
        }
    });
    
    initialized = true;
}

function updateScene(deltaTime) {
    elapsedTime += deltaTime;
    
    // Handle molecule animation if it exists
    if (typeof window.molecule !== 'undefined') {
        // Initialize on first frame
        if (!initialized) {
            initializePositions();
        }
        
        // Gentle rotation of entire molecule
        window.molecule.rotation.y += rotationSpeed * deltaTime;
        window.molecule.rotation.x = Math.sin(elapsedTime * 0.2) * 0.1;
        
        // Animate each atom with appropriate vibration
        window.molecule.children.forEach((child) => {
            // Skip if it's a bond or we don't have original position
            if (child.userData.atomType === 'bond' || !originalPositions.has(child.uuid)) {
                return;
            }
            
            const originalPos = originalPositions.get(child.uuid);
            
            // Determine vibration amount based on atom type
            let vibrationAmount = 0.01; // default
            if (child.userData.atomType === 'carbon') {
                vibrationAmount = carbonVibrationAmount;
            } else if (child.userData.atomType === 'oxygen') {
                vibrationAmount = oxygenVibrationAmount;
            } else if (child.userData.atomType === 'hydrogen') {
                vibrationAmount = hydrogenVibrationAmount;
            }
            
            // Create unique oscillation pattern for each atom
            const offset = child.uuid.charCodeAt(0) % 10; // Use part of UUID for offset
            const time = elapsedTime * vibrationFrequency;
            
            // Apply vibration around original position with 3D variation
            child.position.x = originalPos.x + Math.sin(time + offset) * vibrationAmount;
            child.position.y = originalPos.y + Math.cos(time + offset * 2) * vibrationAmount;
            child.position.z = originalPos.z + Math.sin(time * 1.3 + offset) * vibrationAmount;
        });
        
        // Bonds automatically follow atoms since they're child objects
    }
}

// Hook into animate loop
const originalAnimateAgent = animate;
animate = function() {
    const now = performance.now() * 0.001;
    const deltaTime = now - (window.lastTime || 0);
    window.lastTime = now;

    updateScene(deltaTime);
    originalAnimateAgent();
};
```

## Example: Chemical Reaction Animation

For a prompt like "show a chemical reaction between two molecules":

```javascript
// AnimationAgent code for chemical reaction visualization
let elapsedTime = 0;
const reactionDuration = 10.0; // Seconds the reaction takes
const reactionStartTime = 2.0; // When the reaction begins

// Define reaction states and parameters
const initialDistance = 5.0; // Starting distance between molecule centers
const reactionDistance = 1.0; // Distance during reaction
const finalDistance = 3.0; // End distance after reaction

function updateScene(deltaTime) {
    elapsedTime += deltaTime;
    
    // Handle molecules if they exist
    if (typeof window.molecule1 !== 'undefined' && typeof window.molecule2 !== 'undefined') {
        // Phase 1: Initial approach
        if (elapsedTime < reactionStartTime) {
            const approachProgress = elapsedTime / reactionStartTime;
            const currentDistance = initialDistance * (1 - approachProgress) + reactionDistance * approachProgress;
            
            // Position molecules on opposite sides
            window.molecule1.position.x = -currentDistance / 2;
            window.molecule2.position.x = currentDistance / 2;
            
            // Add some rotation for visual interest
            window.molecule1.rotation.y += 0.2 * deltaTime;
            window.molecule2.rotation.y -= 0.2 * deltaTime;
        }
        // Phase 2: Reaction occurring
        else if (elapsedTime < reactionStartTime + reactionDuration) {
            const reactionProgress = (elapsedTime - reactionStartTime) / reactionDuration;
            
            // Keep molecules at reaction distance
            window.molecule1.position.x = -reactionDistance / 2;
            window.molecule2.position.x = reactionDistance / 2;
            
            // Increase vibration/rotation during reaction
            const reactionIntensity = Math.sin(reactionProgress * Math.PI);
            
            window.molecule1.rotation.y += (0.2 + reactionIntensity) * deltaTime;
            window.molecule2.rotation.y -= (0.2 + reactionIntensity) * deltaTime;
            
            // Optional: change colors during reaction for visual feedback
            if (typeof window.bondReacting !== 'undefined') {
                // Change the color of the reacting bond
                const hue = 240 + reactionProgress * 120; // Shift from blue to green
                window.bondReacting.material.color.setHSL(hue/360, 1, 0.5);
            }
        }
        // Phase 3: Post-reaction separation
        else {
            const separationProgress = Math.min(1.0, (elapsedTime - reactionStartTime - reactionDuration) / 2.0);
            const currentDistance = reactionDistance * (1 - separationProgress) + finalDistance * separationProgress;
            
            // Move molecules to final positions
            window.molecule1.position.x = -currentDistance / 2;
            window.molecule2.position.x = currentDistance / 2;
            
            // Settle into gentle rotation
            window.molecule1.rotation.y += 0.1 * deltaTime;
            window.molecule2.rotation.y -= 0.1 * deltaTime;
        }
    }
}

// Hook into animate loop
const originalAnimateAgent = animate;
animate = function() {
    const now = performance.now() * 0.001;
    const deltaTime = now - (window.lastTime || 0);
    window.lastTime = now;

    updateScene(deltaTime);
    originalAnimateAgent();
};
```

Return your code as a single JavaScript snippet tailored to the user's prompt, with no JSON wrapper or additional explanations.
"""

        llm_response = self.llm_service.generate(prompt_for_llm)
        animation_code = llm_response.strip()
        
        return f"""
// AnimationAgent LLM-generated code
{animation_code}
"""