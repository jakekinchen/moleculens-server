"""
Orchestration Agent - Breaks down scene scripts into discrete objects for visualization.
"""

from agent_management.models import SceneScript, OrchestrationPlan
from agent_management.llm_service import (
    LLMService,
    StructuredLLMRequest,
    LLMModelConfig,
    ProviderType
)

class OrchestrationAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service
    
    def generate_orchestration_plan(self, script: SceneScript) -> OrchestrationPlan:
        """
        Generate an orchestration plan that identifies all discrete objects needed for the scene.
        
        Args:
            script: The scene script to analyze
            
        Returns:
            OrchestrationPlan: A structured plan with all objects needed for the scene
        """
        # Convert the script to a structured format for the LLM prompt
        script_content = "\n\n".join([
            f"TIME {point.timecode}: {point.description}\nCAPTION: {point.caption}"
            for point in script.content
        ])
        
        # Prepare prompt parts
        prompt_intro = f"Analyze this scene script and identify all discrete 3D objects needed for the visualization:"
        prompt_title = f"SCENE TITLE: {script.title}"
        prompt_content = f"SCRIPT CONTENT:\n{script_content}"
        
        prompt_instructions = """
Based on this script, identify all discrete 3D objects that need to be generated for this visualization.
For each object, provide:
1. A unique name
2. A detailed description of what it represents
3. Key properties (size, color, etc.)
4. When it first appears in the scene (timecode in MM:SS format)
5. Relationships to other objects (if any)

Your response should be a complete list of all objects required throughout the entire scene.
Think about what objects are explicitly mentioned as well as what would be needed to create a
cohesive and complete visualization of the described scenarios.

IMPORTANT: For each object, ensure you format properties as JSON objects/dictionaries, not strings.
All properties must be in key-value pairs, and all relationships must be arrays/lists.

For example:
```json
{
  "scene_title": "DNA Structure",
  "objects": [
    {
      "name": "DNA_Helix",
      "description": "The main DNA double helix structure",
      "properties": {
        "color": "blue and red",
        "size": "large",
        "opacity": "85%"
      },
      "appears_at": "00:05",
      "relationships": ["contains DNA_Backbone", "contains BasePairs"]
    }
  ],
  "scene_overview": "Overview of the DNA structure scene"
}
```

Return a JSON object with:
1. scene_title - The title of the scene
2. objects - An array of all discrete objects needed
3. scene_overview - A high-level overview of the overall scene
"""
        
        # Combine the prompt parts
        full_prompt = f"{prompt_intro}\n\n{prompt_title}\n\n{prompt_content}\n\n{prompt_instructions}"
        
        # Create a request for orchestration planning
        request = StructuredLLMRequest(
            user_prompt=full_prompt,
            
            system_prompt="""You are an expert in 3D visualization planning and asset management. Your task is to analyze scene scripts and identify all discrete objects needed for implementation.

When analyzing scripts:
- Identify all explicit and implicit objects needed
- Create clear, unique names for each object
- Provide sufficient description for a 3D artist to model each object
- Specify key properties (dimensions, colors, textures, etc.)
- Note when each object first appears in the timeline
- Identify relationships between objects (e.g., "part of", "interacts with")
- Think holistically about the complete scene requirements

Always return properly structured JSON objects matching the requested schema exactly.
""",
            llm_config=self.llm_service.config,
            response_model=OrchestrationPlan,
        )
        
        # Get the structured response
        return self.llm_service.generate_structured(request)