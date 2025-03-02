"""
Utility functions for handling DeepSeek model responses.
"""

import re
from typing import Optional

def is_deepseek_model(model_name: str) -> bool:
    """
    Check if a model is a DeepSeek model based on its name.
    
    Args:
        model_name: The name of the model to check
        
    Returns:
        True if the model is a DeepSeek model, False otherwise
    """
    return "deepseek" in model_name.lower()

def extract_javascript_from_deepseek_response(response: str) -> str:
    """
    Extract JavaScript code from a DeepSeek model response.
    
    DeepSeek models typically return responses with a <thinking> XML tag block first,
    followed by the actual JavaScript code in a ```javascript code ``` block.
    
    Args:
        response: The raw response from the DeepSeek model
        
    Returns:
        The extracted JavaScript code, or the original response if no code block is found
    """
    # First, try to extract code from a JavaScript code block
    js_match = re.search(r"```javascript\s*([\s\S]*?)\s*```", response)
    if js_match:
        return js_match.group(1).strip()
    
    # Also try to extract code from an HTML block (which may contain JavaScript)
    html_match = re.search(r"```html\s*([\s\S]*?)\s*```", response)
    if html_match:
        return html_match.group(1).strip()
    
    
    # If no code block is found, return the original response
    return response

def extract_structured_output_from_deepseek(response: str, output_type: str = "javascript") -> str:
    """
    Extract structured output from a DeepSeek model response based on the output type.
    
    Args:
        response: The raw response from the DeepSeek model
        output_type: The type of output to extract (e.g., "javascript", "json")
        
    Returns:
        The extracted structured output, or the original response if no matching block is found
    """

    
    # Also check for <think> tags which are used in some DeepSeek models
    think_match = re.search(r"<think>([\s\S]*?)</think>", response)
    if think_match:
        # Remove the think section from the response
        response = response.replace(think_match.group(0), "").strip()
    
    # Extract based on output type
    if output_type.lower() == "javascript":
        return extract_javascript_from_deepseek_response(response)
    elif output_type.lower() == "json":
        # Extract JSON from code blocks
        json_match = re.search(r"```json\s*([\s\S]*?)\s*```", response)
        if json_match:
            return json_match.group(1).strip()
        
        # If no JSON block is found, try to find any code block
        code_match = re.search(r"```\s*([\s\S]*?)\s*```", response)
        if code_match:
            return code_match.group(1).strip()
    
    # If no matching block is found, return the original response
    return response 