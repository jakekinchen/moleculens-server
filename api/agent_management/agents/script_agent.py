"""
Script Agent - Generates structured scene scripts with timecodes, descriptions, and captions.
This is separate from the animation agent that will handle actual animation generation.
"""

from agent_management.models import SceneScript
from agent_management.llm_service import (
    LLMService,
    StructuredLLMRequest,
    LLMModelConfig,
    ProviderType
)

class ScriptAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service
    
    def generate_script(self, topic: str) -> SceneScript:
        """
        Generate a structured scene script for a scientific topic.
        
        Args:
            topic: The scientific topic to explain in the scene
            
        Returns:
            SceneScript: A structured script with title and timed content points
        """
        # Create a request for script generation
        request = StructuredLLMRequest(
            user_prompt=f"""Create a detailed script for a 3D scientific visualization about: {topic}
            
The script should cover a single coherent 3D scene that effectively explains the topic visually.
Focus on creating a cohesive visual narrative with clear transition points.

Return a JSON object with:
1. A concise, descriptive title for the scene
2. A content array containing 5-8 key points in the scene, each with:
   - timecode: Time marker in MM:SS format (starting at 00:00, ending around 02:00)
   - description: Detailed explanation of what should be happening in the 3D scene at this point
   - caption: A concise text caption that would appear on screen at this time (30-50 characters)

Make sure each time point builds logically on the previous ones to tell a complete story about the topic.
Ensure descriptions are specific enough to guide the creation of 3D visuals (objects, movements, transitions).
Keep the explanations scientific but accessible to a general audience.""",
            
            system_prompt="""You are an expert in scientific communication and 3D visualization. Your task is to create structured scripts for educational scientific scenes.

When creating scene scripts:
- Focus on clear, accurate scientific explanations
- Create logical progression of visual elements
- Design effective transitions between key concepts
- Balance visual appeal with educational value
- Use appropriate pacing (approx. 2 minutes total)
- Ensure captions complement but don't duplicate the visuals
- Structure time points to flow naturally from introduction to conclusion

Always return complete, properly formatted JSON objects matching the requested schema exactly.
""",
            llm_config=self.llm_service.config,
            response_model=SceneScript,
        )
        
        # Get the structured response
        return self.llm_service.generate_structured(request)