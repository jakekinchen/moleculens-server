"""
Animation Agent - Generates animation code for Three.js scenes based on script and objects.
Specialized for scientific visualizations including molecular structures.
"""

from typing import Dict, List, Optional, Any
import re
from agent_management.models import SceneScript, OrchestrationPlan, AnimationCode, AnimationKeyframe
from agent_management.llm_service import LLMService, LLMRequest, StructuredLLMRequest
from agent_management.utils.code_extraction import extract_code_block

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
        
        # Extract update function names and object references from geometry code
        geometry_functions = []
        object_references = []
        
        for obj_name, obj_data in object_geometries.items():
            if obj_data.get("status") == "success" and obj_name != "_summary":
                code = obj_data.get("code", "")
                
                # Extract function assignments to window
                function_matches = re.findall(r'window\.(\w+)\s*=\s*function', code)
                for func_name in function_matches:
                    geometry_functions.append(func_name)
                
                # Extract object assignments to window
                obj_matches = re.findall(r'window\.(\w+)\s*=', code)
                for ref_name in obj_matches:
                    if ref_name not in geometry_functions:  # Avoid duplicates with functions
                        object_references.append(ref_name)
        
        # Format the object availability info, including extracted functions and references
        objects_with_geometries = "\n".join([
            f"- {obj_name} (available and ready to animate)"
            for obj_name, obj_data in object_geometries.items()
            if obj_data.get("status") == "success" and obj_name != "_summary"
        ])
        
        geometry_functions_str = "\n".join([
            f"- {func_name} (animation function provided by geometry)"
            for func_name in geometry_functions
        ])
        
        object_references_str = "\n".join([
            f"- {ref_name} (object reference provided by geometry)"
            for ref_name in object_references
        ])
        
        # Prepare the basic prompt components
        title_section = f"## SCENE TITLE\n{script.title}"
        timeline_section = f"## SCRIPT TIMELINE\n{script_timeline}"
        objects_section = f"## AVAILABLE OBJECTS ({geometry_count} objects with generated geometries)\n{objects_with_geometries}"
        complete_objects_section = f"## COMPLETE OBJECT LIST\n{object_list}"
        geometry_functions_section = f"## GEOMETRY FUNCTIONS (use these if provided)\n{geometry_functions_str}" if geometry_functions else ""
        object_references_section = f"## OBJECT REFERENCES (use these with scene.getObjectByName())\n{object_references_str}" if object_references else ""
        
        # Use raw string for the example code to avoid f-string issues
        example_section = r"""
IMPORTANT: Here's ONE high-quality example of self-contained animation code that also shows how to work with any custom functions or update methods provided by the geometry:

```javascript
// Use window.animationTime for elapsed time (includes user time offset from controls)
// This ensures the animation works with the fast-forward and rewind buttons
const elapsedTime = window.animationTime || clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Find main objects in the scene - always use this pattern
const molecule = scene.getObjectByName("molecule");
const energyDiagram = scene.getObjectByName("energy_diagram");

// Helper function for fading objects in and out
function fadeObject(object, fadeIn, duration = 1.0) {
  if (!object) return;
  
  // Set all materials to transparent
  object.traverse(child => {
    if (child.isMesh && child.material) {
      child.material.transparent = true;
      
      // Set target opacity based on fade direction
      if (fadeIn) {
        // Fade in: 0 to 1 over duration
        const progress = Math.min(1, elapsedTime / duration);
        child.material.opacity = progress;
        // Make object visible as soon as we start fading in
        child.visible = true;
      } else {
        // Fade out: 1 to 0 over duration
        const timeRemaining = Math.max(0, 1 - elapsedTime / duration);
        child.material.opacity = timeRemaining;
        // Hide object completely when fully transparent
        if (timeRemaining <= 0) {
          child.visible = false;
        }
      }
    }
  });
}

// ────────────────────────────────────────────────────────────
// At 00:00 - Introduction to Molecule
if (elapsedTime < 15) {
  if (molecule) {
    // IMPORTANT: If the geometry provides an update function, use it
    // This ensures compatibility with how the geometry was designed to animate
    if (typeof window.updateMolecule === 'function') {
      // Call the provided function with appropriate parameters
      window.updateMolecule(elapsedTime);
    } else {
      // Fallback to direct manipulation if no update function exists
      molecule.rotation.y += 0.2 * deltaTime;
    }
    
    // Fade in the molecule using our helper function
    fadeObject(molecule, true, 2.0); // true = fade in, 2.0 = duration in seconds
  }
  
  // Ensure energy diagram is hidden during this phase
  if (energyDiagram) {
    energyDiagram.visible = false;
  }
  
  // Position camera for introduction
  if (camera) {
    const targetPosition = new THREE.Vector3(0, 0, 10);
    camera.position.lerp(targetPosition, 0.05);
    camera.lookAt(0, 0, 0);
  }
}

// ────────────────────────────────────────────────────────────
// At 00:15 - Show Atomic Structure
else if (elapsedTime >= 15 && elapsedTime < 30) {
  if (molecule) {
    // Continue using any update function provided by geometry
    if (typeof window.updateMolecule === 'function') {
      // Pass parameters to control specific behaviors
      const showAtoms = true;
      window.updateMolecule(elapsedTime, showAtoms);
    } else {
      // Direct manipulation fallback
      molecule.rotation.y += 0.1 * deltaTime;
    }
  }
  
  // Introduce the energy diagram with a fade-in effect
  if (energyDiagram) {
    // Calculate time relative to this section
    const sectionTime = elapsedTime - 15;
    
    // Only start showing diagram after 2 seconds into this section
    if (sectionTime >= 2) {
      // Fade in energy diagram
      fadeObject(energyDiagram, true, 1.5);
    }
  }
}

// ────────────────────────────────────────────────────────────
// At 00:30 - Focus on Energy Diagram, Fade Out Molecule
else if (elapsedTime >= 30 && elapsedTime < 45) {
  // Fade out the molecule
  if (molecule) {
    // The false parameter means "fade out"
    fadeObject(molecule, false, 2.0);
  }
  
  // Keep energy diagram fully visible
  if (energyDiagram) {
    // Move energy diagram to center focus
    energyDiagram.position.lerp(new THREE.Vector3(0, 0, 0), 0.05);
    energyDiagram.scale.lerp(new THREE.Vector3(1.5, 1.5, 1.5), 0.05);
  }
}
  }
}

// ────────────────────────────────────────────────────────────
// At 00:30 - Energy Diagram
else if (elapsedTime >= 30 && elapsedTime < 60) {
  // Calculate section progress
  const sectionProgress = (elapsedTime - 30) / 30;
  
  if (molecule) {
    // Keep using update function if available
    if (typeof window.updateMolecule === 'function') {
      window.updateMolecule(elapsedTime);
    } else {
      molecule.rotation.y += 0.05 * deltaTime;
    }
  }
  
  if (energyDiagram) {
    // For the energy diagram, fade it in
    energyDiagram.traverse(child => {
      if (child.isMesh && child.material) {
        child.material.opacity = Math.min(1, sectionProgress * 2);
        child.material.transparent = child.material.opacity < 1;
      }
    });
    
    // Position it properly
    energyDiagram.position.lerp(new THREE.Vector3(5, 0, 0), 0.05);
  }
}
```

Note how this example:
1. Checks for objects using scene.getObjectByName()
2. Properly handles update functions provided by the geometry
3. Provides fallbacks for direct manipulation when no update functions exist
4. Creates smooth time-based animations with clear section dividers
"""
        
        # Combine everything into the final prompt
        prompt = f"""You are a Three.js animation expert. Your task is to create simple, self-contained, and foolproof animation code for a scientific visualization.

{title_section}

{timeline_section}

{objects_section}

{complete_objects_section}

{geometry_functions_section}

{object_references_section}

IMPORTANT STRUCTURE INFORMATION:
The script contains 5-8 key points, each with:
- timecode: Time marker in MM:SS format (starting at 00:00, ending around 02:00)
- description: Detailed explanation of what should be happening in the 3D scene at this point
- caption: A concise text caption that would appear on screen at this time (30-50 characters)

Your task is to write the animation code that would go INSIDE the animate() function in Three.js.
DO NOT include the function declaration itself - only provide the code that would go inside the function.

CRITICAL REQUIREMENTS:
1. If the geometry provides update functions (listed above), USE THEM INSTEAD of directly manipulating objects
2. Use scene.getObjectByName() to find objects by name
3. ALWAYS check if objects exist before operating on them
4. Use clear time ranges (if elapsedTime >= start && elapsedTime < end) for each script section
5. Each timecode section should be clearly separated with comments
6. Include section dividers for readability (like ────────────────────────────────────)

VERY IMPORTANT:
- The geometry code might have its own animation/update functions for complex behaviors
- If window.updateSomething functions are listed above, call them instead of trying to manipulate objects directly
- It's ok to use window.functionName() to call these functions, but DON'T use window.objectName to access objects

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
- IMPORTANT: Objects that aren't actively used in a scene should not be visible until needed
- IMPORTANT: Set visible=false for objects that aren't currently relevant (except main objects)
- IMPORTANT: Fade objects in (opacity 0→1) when introducing them and fade out (opacity 1→0) when no longer needed
- Use opacity transitions of 1-2 seconds for smooth appearance/disappearance
"""

        # Create the request
        request = LLMRequest(
            user_prompt=prompt,
            system_prompt="You are an expert Three.js animator specializing in scientific visualizations. Create detailed, timeline-based animation code following the script exactly. Your code should be simple, self-contained, and foolproof, with proper handling of geometry-provided update functions.",
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
        Generate Three.js animation code using LLM.
        """
        prompt_for_llm = f"You are an expert Three.js animator specializing in scientific visualizations, particularly molecular structures and chemical reactions. The user wants to animate a scene with the following prompt: '{user_prompt}'\n\n" + r"""
Requirements:
- Generate JavaScript animation code for objects created by the GeometryAgent.
- Your code will be injected into an existing Three.js scene with a standard animation loop.
- Create scientifically accurate yet visually appealing animations that represent the prompt.
- Implement clear, self-contained animation logic.
- Use deltaTime parameter for smooth, frame-rate independent animations.
- For chemical/molecular visualizations, consider appropriate molecular behaviors.
- IMPORTANT: Objects that aren't actively used should not be visible until needed.
- Set visible=false for objects that aren't currently relevant.
- Fade objects in (opacity 0→1) when introducing them and fade out (opacity 1→0) when no longer needed.
- Use opacity transitions of 1-2 seconds for smooth appearance/disappearance.

# Guide to Creating Self-Contained Animations

When creating animations in Three.js:
1. Use scene.getObjectByName() to access objects rather than global window references
2. Always check if objects exist before manipulating them
3. If the geometry has provided update functions (like window.updateMolecule), USE THEM
4. Avoid HTML DOM interactions
5. Use clock.getElapsedTime() and clock.getDelta() for timing

## Example: Self-Contained Molecular Animation with Update Function

```javascript
// Use window.animationTime for elapsed time (includes user time offset from controls)
// This ensures the animation works with the fast-forward and rewind buttons
const elapsedTime = window.animationTime || clock.getElapsedTime();
const deltaTime = clock.getDelta();

// Find the molecule in the scene
const molecule = scene.getObjectByName("molecule");
if (!molecule) return; // Exit early if molecule not found

// Check if there's a dedicated update function provided by the geometry
const hasUpdateFunction = typeof window.updateMolecule === 'function';

// Define animation phases based on elapsed time
if (elapsedTime < 5) {
  // Phase 1: Introduction (0-5 seconds)
  if (hasUpdateFunction) {
    // Use the provided update function with parameters that make sense for this phase
    window.updateMolecule(elapsedTime, { phase: "intro" });
  } else {
    // Fallback: manually animate if no update function exists
    // Slowly rotate and scale up the molecule
    const introProgress = elapsedTime / 5;
    molecule.rotation.y += 0.5 * deltaTime;
    molecule.scale.setScalar(0.2 + introProgress * 0.8); // Scale from 0.2 to 1.0
    
    // Fade in the molecule
    molecule.traverse(child => {
      if (child.isMesh && child.material) {
        // Set material to transparent for fade effect
        child.material.transparent = true;
        // Gradually increase opacity from 0 to 1 over 2 seconds
        child.material.opacity = Math.min(1, elapsedTime / 2);
      }
    });
  }
  
  // Position camera for a good view
  if (camera) {
    camera.position.lerp(new THREE.Vector3(0, 2, 10), 0.05);
    camera.lookAt(molecule.position);
  }
} 
else if (elapsedTime < 15) {
  // Phase 2: Highlight atoms (5-15 seconds)
  if (hasUpdateFunction) {
    // Pass the appropriate parameters to the update function
    window.updateMolecule(elapsedTime, { phase: "highlight", highlight: "atoms" });
  } else {
    // Manual fallback animation
    molecule.rotation.y += 0.2 * deltaTime;
    
    // Highlight atoms
    molecule.traverse(child => {
      if (!child.isMesh) return;
      
      // Highlight logic...
      // (simplified for brevity)
    });
  }
  
  // Move camera
  if (camera) {
    // Camera movement...
  }
}
else {
  // Phase 3: Show molecular motion (after 15 seconds)
  if (hasUpdateFunction) {
    window.updateMolecule(elapsedTime, { phase: "motion" });
  } else {
    // Manual fallback animation
    molecule.rotation.y += 0.3 * deltaTime;
    
    // Vibration logic...
    // (simplified for brevity)
  }
  
  // Move camera
  if (camera) {
    // Camera movement...
  }
}
```

Return your code as a single JavaScript snippet tailored to the user's prompt, with no JSON wrapper or additional explanations.
"""

        llm_response = self.llm_service.generate(prompt_for_llm)
        animation_code = extract_code_block(llm_response.content, "javascript")
        
        return f"""
// AnimationAgent LLM-generated code
{animation_code}
"""