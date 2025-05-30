"""
Caption Agent - Uses an LLM to generate text captions, then returns JS snippet injecting those captions.
"""
import json
from typing import List, Dict, Any
from agent_management.llm_service import LLMService
from agent_management.providers.deepseek_utils import is_deepseek_model, extract_structured_output_from_deepseek

class CaptionAgent:
    def __init__(self, llm_service: LLMService):
        self.llm_service = llm_service

    def get_caption_snippet(self, user_prompt: str) -> str:
        """
        Call the LLM to generate an array of { time, text } captions,
        then return a code snippet that updates a #caption <div> accordingly.
        """
        # Prompt the LLM for some short creative captions
        llm_prompt = f"""
        Generate a JSON array of simple 'timestamped' captions (time in seconds, text string)
        describing an animation for: '{user_prompt}'.
        Example:
        [
          {{ "time": 0, "text": "Behold the starting moment..." }},
          {{ "time": 5, "text": "An incredible transformation is underway..." }}
        ]
        """
        response = self.llm_service.generate(llm_prompt)
        
        # Check if we're using a DeepSeek model and apply special extraction if needed
        content = response.content
        if is_deepseek_model(response.model):
            content = extract_structured_output_from_deepseek(content, "json")
        
        try:
            # Ensure response is a string before passing to json.loads
            captions_data: List[Dict[str, Any]] = json.loads(str(content))
        except:
            # fallback if invalid
            captions_data = [
                {"time": 0, "text": "Default start caption"},
                {"time": 5, "text": "Default second caption"}
            ]
        
        # We'll convert the array into JS code.
        # We'll reference an HTML <div id="caption"></div> for display.
        return f"""
// CaptionAgent code
const captionElement = document.getElementById('caption');
const captionStartTime = performance.now() / 1000;
let lastCaption = '';

const captions = {json.dumps(captions_data, indent=2)};

function updateCaptions() {{
    const elapsed = (performance.now() / 1000) - captionStartTime;
    let current = '';
    for (const c of captions) {{
        if (elapsed >= c.time) {{
            current = c.text;
        }}
    }}
    if (current !== lastCaption) {{
        captionElement.textContent = current;
        lastCaption = current;
    }}
}}

// Insert caption updates into animate loop
const originalAnimateCaptions = animate;
animate = function() {{
    originalAnimateCaptions();
    updateCaptions();
}};
"""