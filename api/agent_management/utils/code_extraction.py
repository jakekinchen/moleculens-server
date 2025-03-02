"""
Code extraction utilities for handling LLM responses.
"""
import re
from typing import Optional

def extract_code_block(content: str, language: Optional[str] = None) -> str:
    """
    Extract code from a response, removing code block markers and thinking sections.
    
    Args:
        content: The response content to extract code from
        language: Optional language to look for (e.g. "javascript", "python")
        
    Returns:
        The extracted code with markers removed
    """
    
    # Try language-specific code block
    if language:
        pattern = f"```{language}\s*(.*?)\s*```"
        match = re.search(pattern, content, re.DOTALL)
        if match:
            return match.group(1).strip()
    
    
    # If no code block found, return cleaned content
    return content.strip() 