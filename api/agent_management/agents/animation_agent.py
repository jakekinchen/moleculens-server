"""
Animation Agent - Generates animation code for Three.js scenes based on script and objects.
Specialized for scientific visualizations including molecular structures.
"""

from typing import Dict, List, Optional, Any
import re
from agent_management.models import SceneScript, OrchestrationPlan, AnimationCode, AnimationKeyframe
from agent_management.llm_service import LLMService, LLMRequest, StructuredLLMRequest

class AnimationAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service
        
    def generate_animation_code(self, 
                               script: SceneScript, 
                               object_geometries: Dict[str, Dict],
                               orchestration_plan: OrchestrationPlan) -> AnimationCode:
        """
        Generate Three.js animation code based on the scene script and generated geometries.
        
        Args:
            script: The scene script with timecodes and descriptions
            object_geometries: Dictionary of object names to their geometry code and metadata
            orchestration_plan: The orchestration plan with object details
            
        Returns:
            AnimationCode: Animation code with keyframes
        """
        # Format script timeline
        script_timeline = "\n\n".join([
            f"TIMECODE {point.timecode}:\nDESCRIPTION: {point.description}\nCAPTION: {point.caption}"
            for point in script.content
        ])
        
        # Prepare a list of available objects with their appearance times
        object_list = "\n".join([
            f"- {obj.name} (appears at {obj.appears_at}): {obj.description[:100]}..." 
            if isinstance(obj.description, str) and len(obj.description) > 100 else
            f"- {obj.name} (appears at {obj.appears_at}): {obj.description}"
            for obj in orchestration_plan.objects
        ])
        
        # Count available objects with geometries
        geometry_count = sum(1 for obj_data in object_geometries.values() 
                           if obj_data.get("status") == "success")
        
        # Format the object availability info
        objects_with_geometries = "\n".join([
            f"- {obj_name} (available and ready to animate)"
            for obj_name, obj_data in object_geometries.items()
            if obj_data.get("status") == "success" and obj_name != "_summary"
        ])
        
        # Prepare the basic prompt components
        title_section = f"## SCENE TITLE\n{script.title}"
        timeline_section = f"## SCRIPT TIMELINE\n{script_timeline}"
        objects_section = f"## AVAILABLE OBJECTS ({geometry_count} objects with generated geometries)\n{objects_with_geometries}"
        complete_objects_section = f"## COMPLETE OBJECT LIST\n{object_list}"
        
        # Use raw string for the example code to avoid f-string issues
        example_section = r"""
IMPORTANT: Here's ONE high-quality example of self-contained animation code that follows a script:

```javascript
// Get elapsed time and delta time
const elapsedTime = clock.getElapsedTime();
const deltaTime = clock.getDelta();

// ────────────────────────────────────────────────────────────
// At 00:00 - Introduction to Protein Folding
if (elapsedTime < 15) {
  // Find the protein in the scene
  const protein = scene.getObjectByName("protein_model");
  if (protein) {
    // Slow rotation to introduce the protein
    protein.rotation.y += 0.2 * deltaTime;
    
    // Fade in the protein if it has an opacity property
    protein.traverse(child => {
      if (child.isMesh && child.material && child.material.opacity !== undefined) {
        // Gradually increase opacity over 5 seconds
        child.material.opacity = Math.min(1, elapsedTime / 5);
        child.material.transparent = child.material.opacity < 1;
      }
    });
  }
  
  // Position camera for introduction
  if (camera) {
    // Start with a wider view
    const targetPosition = new THREE.Vector3(0, 0, 10);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(0, 0, 0);
  }
}

// ────────────────────────────────────────────────────────────
// At 00:15 - Primary Structure Highlight
else if (elapsedTime >= 15 && elapsedTime < 30) {
  // Calculate progress within this section
  const sectionProgress = (elapsedTime - 15) / 15;
  
  // Find relevant objects
  const protein = scene.getObjectByName("protein_model");
  const aminoAcids = scene.getObjectByName("amino_acids");
  
  if (protein && aminoAcids) {
    // Highlight the amino acid sequence
    aminoAcids.traverse(child => {
      if (child.isMesh && child.material) {
        // Pulse the emissive property
        const pulse = Math.sin(elapsedTime * 2) * 0.5 + 0.5;
        child.material.emissive = new THREE.Color(0x003366);
        child.material.emissiveIntensity = pulse;
      }
    });
    
    // Move camera to focus on the amino acid chain
    if (camera) {
      const targetPosition = new THREE.Vector3(-2, 0, 5);
      camera.position.lerp(targetPosition, 0.05);
      camera.lookAt(-2, 0, 0);
    }
  }
}

// ────────────────────────────────────────────────────────────
// At 00:30 - Secondary Structure Formation
else if (elapsedTime >= 30 && elapsedTime < 60) {
  // Calculate section progress
  const sectionProgress = (elapsedTime - 30) / 30;
  
  // Find the secondary structures
  const alphaHelix = scene.getObjectByName("alpha_helix");
  const betaSheet = scene.getObjectByName("beta_sheet");
  
  if (alphaHelix && betaSheet) {
    // Alpha helix forms first (0-15s of this section)
    if (sectionProgress < 0.5) {
      const helixProgress = sectionProgress * 2; // 0-1 during first half
      
      // Gradually reveal the alpha helix
      alphaHelix.traverse(child => {
        if (child.isMesh && child.material) {
          child.material.opacity = helixProgress;
          child.material.transparent = child.material.opacity < 1;
        }
      });
      
      // Rotate to give 3D perspective of the helix
      alphaHelix.rotation.y += 0.3 * deltaTime;
    } 
    // Beta sheet forms second (15-30s of this section)
    else {
      // Keep alpha helix visible
      alphaHelix.traverse(child => {
        if (child.isMesh && child.material) {
          child.material.opacity = 1;
          child.material.transparent = false;
        }
      });
      
      // Gradually form the beta sheet
      const sheetProgress = (sectionProgress - 0.5) * 2; // 0-1 during second half
      betaSheet.traverse(child => {
        if (child.isMesh && child.material) {
          child.material.opacity = sheetProgress;
          child.material.transparent = child.material.opacity < 1;
          
          // Expand the sheet to its full size
          child.scale.set(
            1,
            Math.min(1, sheetProgress * 1.5), // Slightly exaggerated for effect
            1
          );
        }
      });
    }
    
    // Camera movement to showcase the structures
    if (camera) {
      // Pan around to view different angles
      const angle = sectionProgress * Math.PI * 2; // Full 360° rotation
      const radius = 7 - sectionProgress * 2; // Move closer over time
      const cameraX = Math.cos(angle) * radius;
      const cameraZ = Math.sin(angle) * radius;
      
      camera.position.lerp(new THREE.Vector3(cameraX, 1, cameraZ), 0.02);
      camera.lookAt(0, 0, 0);
    }
  }
}

// Continue with remaining script sections...
```

Note how this example:
1. Divides animation into clear time-based sections matching script timecodes
2. Finds objects using scene.getObjectByName() instead of window references
3. Checks if objects exist before using them
4. Creates smooth transitions between states
5. Handles camera movements to highlight important elements
6. Uses clear section dividers and comments
"""
        
        # Combine everything into the final prompt
        prompt = f"""You are a Three.js animation expert. Your task is to create simple, self-contained, and foolproof animation code for a scientific visualization.

{title_section}

{timeline_section}

{objects_section}

{complete_objects_section}

IMPORTANT STRUCTURE INFORMATION:
The script contains 5-8 key points, each with:
- timecode: Time marker in MM:SS format (starting at 00:00, ending around 02:00)
- description: Detailed explanation of what should be happening in the 3D scene at this point
- caption: A concise text caption that would appear on screen at this time (30-50 characters)

Your task is to write the animation code that would go INSIDE the animate() function in Three.js.
DO NOT include the function declaration itself - only provide the code that would go inside the function.

CRITICAL REQUIREMENTS:
1. DO NOT use window.* references to access objects
2. ALWAYS use scene.getObjectByName() to find objects by name
3. ALWAYS check if objects exist before operating on them
4. Use clear time ranges (if elapsedTime >= start && elapsedTime < end) for each script section
5. Rely only on scene, camera, clock, and objects passed directly to the animate function
6. Each timecode section should be clearly separated with comments
7. Include section dividers for readability (like ────────────────────────────────────)

The animation should:
1. Respect the timecodes in the script
2. Add proper animations and transitions based on the script
3. Include camera movements where appropriate
4. Stage the entrance and exit of objects according to their timecodes
5. Create a cohesive visual experience that follows the script's narrative
{example_section}

Return ONLY the animation code that would go inside the animate function, not the function declaration itself.
Include comments with timecodes to make the code clear and maintainable.

REMEMBER:
- Use clear timecode comments (// At MM:SS - Description) throughout the code
- Don't include the animate function declaration 
- Always check if objects exist before using them
- Use deltaTime for smooth animations (multiply by deltaTime)
- DO NOT use any HTML references or window properties
"""

        # Create the request
        request = LLMRequest(
            user_prompt=prompt,
            system_prompt="You are an expert Three.js animator specializing in scientific visualizations. Create detailed, timeline-based animation code following the script exactly. Your code should be simple, self-contained, and foolproof, without relying on HTML references or window properties.",
            llm_config=self.llm_service.config
        )
        
        # Generate the animation code
        response = self.llm_service.generate(request)
        animation_code = response.content.strip()
        
        # Extract keyframes from the code by parsing comments with timecodes
        keyframes = self._extract_keyframes(animation_code)
        
        return AnimationCode(
            code=animation_code,
            keyframes=keyframes
        )
    
    def _extract_keyframes(self, code: str) -> List[AnimationKeyframe]:
        """
        Extract keyframes from animation code by parsing comments with timecodes.
        """
        keyframes = []
        # Look for comments with timecodes in MM:SS format
        pattern = r'//.*?(\d{2}:\d{2}).*?\n(.*?)(?=//|\n\n|$)'
        matches = re.findall(pattern, code, re.DOTALL)
        
        for timecode, actions_block in matches:
            # Clean up the actions text and split into separate actions
            actions = [line.strip() for line in actions_block.split('\n') 
                      if line.strip() and not line.strip().startswith('//')]
            
            if actions:
                keyframes.append(AnimationKeyframe(
                    timecode=timecode,
                    actions=actions
                ))
        
        # Sort keyframes by timecode
        keyframes.sort(key=lambda k: k.timecode)
        
        return keyframes
        
    def get_animation_snippet(self, user_prompt: str) -> str:
        """
        Legacy method for generating Three.js animation code using LLM based on user prompt.
        This code will animate the objects created by the GeometryAgent,
        with special focus on scientific visualizations like molecular structures.
        """
        prompt_for_llm = f"You are an expert Three.js animator specializing in scientific visualizations, particularly molecular structures and chemical reactions. The user wants to animate a scene with the following prompt: '{user_prompt}'\n\n" + r"""
Requirements:
- Generate JavaScript animation code for objects created by the GeometryAgent.
- DO NOT use global window references.
- Your code will be injected into an existing Three.js scene with a standard animation loop.
- Create scientifically accurate yet visually appealing animations that represent the prompt.
- Implement clear, self-contained animation logic.
- Use deltaTime parameter for smooth, frame-rate independent animations.
- For chemical/molecular visualizations, consider appropriate molecular behaviors.

# Guide to Creating Self-Contained Animations

When creating animations in Three.js:
1. Use scene.getObjectByName() to access objects rather than global window references
2. Always check if objects exist before manipulating them
3. Keep animation logic simple and focused on the scene graph
4. Avoid HTML DOM interactions or window properties
5. Use clock.getElapsedTime() and clock.getDelta() for timing

## Example: Self-Contained Molecular Animation

```javascript
// Get elapsed time and delta time
const elapsedTime = clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Find the molecule in the scene
const molecule = scene.getObjectByName("molecule");
if (!molecule) return; // Exit early if molecule not found

// Define animation phases based on elapsed time
if (elapsedTime < 5) {
  // Phase 1: Introduction (0-5 seconds)
  // Slowly rotate and scale up the molecule
  const introProgress = elapsedTime / 5;
  molecule.rotation.y += 0.5 * deltaTime;
  molecule.scale.setScalar(0.2 + introProgress * 0.8); // Scale from 0.2 to 1.0
  
  // Position camera for a good view
  if (camera) {
    camera.position.lerp(new THREE.Vector3(0, 2, 10), 0.05);
    camera.lookAt(molecule.position);
  }
} 
else if (elapsedTime < 15) {
  // Phase 2: Highlight atoms (5-15 seconds)
  molecule.rotation.y += 0.2 * deltaTime;
  
  // Find atoms based on their names
  molecule.traverse(child => {
    if (!child.isMesh) return;
    
    // Determine atom type by name
    const isCarbon = child.name.includes('carbon');
    const isOxygen = child.name.includes('oxygen');
    const isHydrogen = child.name.includes('hydrogen');
    
    // Skip bonds
    if (!isCarbon && !isOxygen && !isHydrogen) return;
    
    // Calculate pulse effect
    const pulseTime = elapsedTime - 5; // Adjusted time
    const pulseSpeed = isHydrogen ? 3 : (isOxygen ? 2 : 1); // Different speeds
    const pulseAmount = 0.1 * Math.sin(pulseTime * pulseSpeed);
    
    // Store original scale if not already saved
    if (!child.userData.originalScale) {
      child.userData.originalScale = child.scale.x;
    }
    
    // Apply pulsing effect
    const baseScale = child.userData.originalScale;
    child.scale.setScalar(baseScale + pulseAmount);
    
    // Highlight specific atom types at different times
    if ((isHydrogen && pulseTime < 3) || 
        (isOxygen && pulseTime >= 3 && pulseTime < 7) ||
        (isCarbon && pulseTime >= 7)) {
      // Highlight active atoms with emissive color
      if (child.material) {
        child.material.emissive = new THREE.Color(
          isHydrogen ? 0x00ffff : (isOxygen ? 0xff0000 : 0x00ff00)
        );
        child.material.emissiveIntensity = 0.3 + 0.7 * Math.sin(pulseTime * 2);
      }
    } else if (child.material) {
      // Reset non-active atoms
      child.material.emissiveIntensity = 0;
    }
  });
  
  // Move camera for a closer look at specific parts
  if (camera) {
    const targetPos = new THREE.Vector3(
      Math.sin(elapsedTime * 0.2) * 5,
      2 + Math.cos(elapsedTime * 0.3),
      8 - (elapsedTime - 5) * 0.3 // Gradually move closer
    );
    camera.position.lerp(targetPos, 0.02);
    camera.lookAt(molecule.position);
  }
}
else {
  // Phase 3: Show molecular motion (after 15 seconds)
  molecule.rotation.y += 0.3 * deltaTime;
  
  // Simulate molecular vibration
  molecule.traverse(child => {
    if (!child.isMesh) return;
    
    // Skip large parent objects
    if (child.children.length > 0) return;
    
    // Store original position if not already saved
    if (!child.userData.originalPosition) {
      child.userData.originalPosition = child.position.clone();
    }
    
    const originalPos = child.userData.originalPosition;
    const vibrationAmount = 0.03;
    const vibrationSpeed = 3;
    
    // Create unique oscillation pattern for each atom
    const offset = child.id % 10; // Use object ID for offset
    
    // Apply vibration around original position
    child.position.x = originalPos.x + Math.sin(elapsedTime * vibrationSpeed + offset) * vibrationAmount;
    child.position.y = originalPos.y + Math.cos(elapsedTime * vibrationSpeed + offset * 2) * vibrationAmount;
    child.position.z = originalPos.z + Math.sin(elapsedTime * vibrationSpeed * 1.3 + offset) * vibrationAmount;
    
    // Reset any highlights
    if (child.material && child.material.emissiveIntensity) {
      child.material.emissiveIntensity = 0;
    }
  });
  
  // Move camera to show the whole molecule in motion
  if (camera) {
    const orbitRadius = 12;
    const orbitSpeed = 0.1;
    const orbitY = 5;
    const targetPos = new THREE.Vector3(
      Math.cos(elapsedTime * orbitSpeed) * orbitRadius,
      orbitY,
      Math.sin(elapsedTime * orbitSpeed) * orbitRadius
    );
    camera.position.lerp(targetPos, 0.01);
    camera.lookAt(molecule.position);
  }
}
```

Return your code as a single JavaScript snippet tailored to the user's prompt, with no JSON wrapper or additional explanations.
"""

        llm_response = self.llm_service.generate(prompt_for_llm)
        animation_code = llm_response.content.strip()
        
        return f"""
// AnimationAgent LLM-generated code
{animation_code}
"""