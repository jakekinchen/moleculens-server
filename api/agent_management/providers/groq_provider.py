"""
Groq provider implementation for LLM service.
"""

import json
import os
import re
from typing import Dict, Any, TypeVar, Type, List, Optional, Union, cast
from ..llm_service import LLMProvider, LLMRequest, LLMResponse, StructuredLLMRequest, T, MessageRole
from .deepseek_utils import is_deepseek_model, extract_structured_output_from_deepseek

class GroqProvider(LLMProvider):
    """Groq-specific implementation"""
    
    def __init__(self, api_key: Optional[str] = None):
        try:
            import groq
        except ImportError:
            raise ImportError(
                "The 'groq' package is not installed. "
                "Please install it using: pip install groq"
            )
            
        self.client = groq.Groq(
            api_key=api_key or os.getenv("GROQ_API_KEY")
        )
        
        if not api_key and not os.getenv("GROQ_API_KEY"):
            raise ValueError(
                "Groq API key is not provided and GROQ_API_KEY environment variable is not set"
            )

    def _convert_messages(self, request: LLMRequest) -> List[Dict[str, str]]:
        """Convert our message format to Groq's format"""
        messages = []
        if request.system_prompt:
            messages.append({
                "role": "system",
                "content": request.system_prompt
            })
        messages.append({
            "role": "user",
            "content": request.user_prompt
        })
        return messages

    def generate(self, request: LLMRequest) -> LLMResponse:
        """Generate a response using Groq's API"""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)
            
            # Build parameters dictionary
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "stream": request.stream
            }
            
            # Add temperature if provided
            if request.temperature is not None:
                params["temperature"] = request.temperature
            
            # Add max_tokens if provided
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
                
            # Add top_p if provided
            if request.top_p is not None:
                params["top_p"] = request.top_p
                
            # Add any additional parameters
            params.update(request.additional_params)
                
            response = self.client.chat.completions.create(**params)
            
            if not response.choices or not response.choices[0].message or not response.choices[0].message.content:
                raise ValueError("No response content received from Groq")
            
            return LLMResponse(
                content=response.choices[0].message.content,
                model=request.llm_config.model_name,
                usage={
                    "prompt_tokens": response.usage.prompt_tokens if hasattr(response.usage, "prompt_tokens") else 0,
                    "completion_tokens": response.usage.completion_tokens if hasattr(response.usage, "completion_tokens") else 0,
                    "total_tokens": response.usage.total_tokens if hasattr(response.usage, "total_tokens") else 0
                }
            )
        except Exception as e:
            raise Exception(f"Groq API error: {str(e)}")

    def _extract_json_from_text(self, text: str) -> str:
        """
        Extract JSON from text, handling various formats the model might return.
        
        Args:
            text: The text containing JSON
            
        Returns:
            The extracted JSON string
        """
        # First, try to find JSON within code blocks
        if "```json" in text:
            # Extract content between ```json and ```
            match = re.search(r"```json\s*([\s\S]*?)\s*```", text)
            if match:
                return match.group(1).strip()
        elif "```" in text:
            # Extract content between ``` and ```
            match = re.search(r"```\s*([\s\S]*?)\s*```", text)
            if match:
                return match.group(1).strip()
                
        # If no code blocks, look for JSON objects directly
        # This regex looks for content between { and } including nested braces
        match = re.search(r"(\{[\s\S]*\})", text)
        if match:
            return match.group(1).strip()
            
        # If we still haven't found JSON, check if the model prefixed with text like
        # "Here is a JSON object that conforms to the schema:"
        if "Here is" in text and "{" in text:
            # Extract everything from the first { to the end of the text
            match = re.search(r"(\{[\s\S]*)", text)
            if match:
                # Now find the balanced JSON object
                json_text = match.group(1)
                # Count opening and closing braces to find the complete JSON object
                open_count = 0
                close_count = 0
                for i, char in enumerate(json_text):
                    if char == '{':
                        open_count += 1
                    elif char == '}':
                        close_count += 1
                        if open_count == close_count:
                            return json_text[:i+1].strip()
                            
        # If all else fails, return the original text
        return text

    def _fix_json_content(self, json_content: str) -> str:
        """
        Fix common issues with JSON content returned by the model.
        
        Args:
            json_content: The JSON content to fix
            
        Returns:
            Fixed JSON content
        """
        # Replace JavaScript Math expressions with their values
        json_content = re.sub(r'Math\.PI\s*/\s*2', '1.5707963267948966', json_content)
        json_content = re.sub(r'Math\.PI', '3.141592653589793', json_content)
        
        # Replace any other JavaScript expressions
        json_content = re.sub(r'Math\.[a-zA-Z]+\([^)]*\)', '0', json_content)
        
        # Fix common syntax errors
        # Replace single quotes with double quotes for property names
        json_content = re.sub(r"'([^']+)':", r'"\1":', json_content)
        
        # Replace single quotes with double quotes for property values
        # First, handle the case where the entire JSON is using single quotes
        if json_content.strip().startswith("'") and json_content.strip().endswith("'"):
            # This is a JSON string wrapped in single quotes
            json_content = json_content.strip()[1:-1]
        
        # Handle single-quoted strings within the JSON
        # This is more complex as we need to avoid replacing quotes within already double-quoted strings
        # We'll use a simple approach that works for most cases
        try:
            # Try to parse as is - if it works, no need for fixes
            json.loads(json_content)
        except json.JSONDecodeError:
            # If it fails, try to replace single quotes with double quotes
            # This regex finds single-quoted strings not inside double quotes
            # It's a simplified approach and may not work for all edge cases
            json_content = re.sub(r"(?<!\\)'([^']*?)(?<!\\)'", r'"\1"', json_content)
        
        # Remove trailing commas in arrays and objects
        json_content = re.sub(r',\s*}', '}', json_content)
        json_content = re.sub(r',\s*]', ']', json_content)
        
        return json_content

    def generate_structured(self, request: StructuredLLMRequest[T]) -> T:
        """Generate a structured response using Groq's API"""
        try:
            if not request.llm_config:
                raise ValueError("LLM configuration is required")

            messages = self._convert_messages(request)
            
            # Append instruction to format response as JSON
            messages[-1]["content"] += "\n\nYou must respond with a valid JSON object that conforms to the following schema:\n"
            schema = request.response_model.model_json_schema()
            messages[-1]["content"] += json.dumps(schema, indent=2)
            messages[-1]["content"] += "\n\nIMPORTANT INSTRUCTIONS:"
            messages[-1]["content"] += "\n1. Respond ONLY with a valid JSON instance, not the schema itself."
            messages[-1]["content"] += "\n2. Do not include $schema, $defs, or any schema-related fields in your response."
            messages[-1]["content"] += "\n3. Do not use JavaScript expressions like Math.PI - use the actual numeric values instead."
            messages[-1]["content"] += "\n4. Do not include any explanatory text before or after the JSON."
            
            # Build parameters dictionary
            params: Dict[str, Any] = {
                "model": request.llm_config.model_name,
                "messages": messages,
                "stream": request.stream
            }
            
            # Add temperature if provided
            if request.temperature is not None:
                params["temperature"] = request.temperature
            
            # Add max_tokens if provided
            if request.max_tokens is not None:
                params["max_tokens"] = request.max_tokens
                
            # Add top_p if provided
            if request.top_p is not None:
                params["top_p"] = request.top_p
                
            # Add any additional parameters
            params.update(request.additional_params)
                
            response = self.client.chat.completions.create(**params)
            
            if not response.choices or not response.choices[0].message or not response.choices[0].message.content:
                raise ValueError("No response content received from Groq")
            
            content = response.choices[0].message.content
            
            # Check if we're using a DeepSeek model and apply special extraction if needed
            if is_deepseek_model(request.llm_config.model_name):
                content = extract_structured_output_from_deepseek(content, "json")
            
            # Try to extract JSON from the response
            try:
                # Extract JSON from the response
                json_content = self._extract_json_from_text(content)
                
                # Fix common issues with JSON
                json_content = self._fix_json_content(json_content)
                
                # Check if the response is the schema itself rather than an instance
                if "$schema" in json_content or "$defs" in json_content:
                    # The model returned the schema instead of an instance
                    # Try to extract a valid instance from the "properties" field if it exists
                    try:
                        parsed = json.loads(json_content)
                        if "properties" in parsed and isinstance(parsed["properties"], dict):
                            # Look for a nested object that might be the actual instance
                            for key, value in parsed.items():
                                if key not in ["$schema", "$defs", "properties", "required", "title", "type"]:
                                    if isinstance(value, dict) and "properties" not in value:
                                        # This might be our instance
                                        json_content = json.dumps(value)
                                        break
                    except Exception:
                        # If we can't extract an instance, continue with the original content
                        pass
                
                # Parse the JSON
                json_response = json.loads(json_content)
                
                # Validate against the model
                result = request.response_model.model_validate(json_response)
                return result
            except json.JSONDecodeError as e:
                # If we still have a JSON decode error, try a more aggressive approach
                try:
                    # Try to manually fix the JSON by replacing Math.PI expressions
                    # This is a more targeted approach for the specific error we're seeing
                    if "Math.PI" in json_content:
                        # Find the specific line with Math.PI
                        lines = json_content.split('\n')
                        for i, line in enumerate(lines):
                            if "Math.PI" in line:
                                # Replace Math.PI with its value
                                lines[i] = line.replace("Math.PI", "3.141592653589793")
                        
                        # Join the lines back together
                        json_content = '\n'.join(lines)
                        
                        # Try parsing again
                        json_response = json.loads(json_content)
                        result = request.response_model.model_validate(json_response)
                        return result
                    else:
                        # If there's no Math.PI, try a more general approach
                        # Sometimes the model returns invalid JSON with JavaScript-like syntax
                        # Try to fix it by evaluating it as a Python dictionary
                        try:
                            # This is a last resort and potentially unsafe
                            # Only use for testing purposes
                            import ast
                            # Try to parse as a Python literal
                            fixed_json = ast.literal_eval(json_content)
                            # Convert to JSON
                            json_response = json.loads(json.dumps(fixed_json))
                            result = request.response_model.model_validate(json_response)
                            return result
                        except Exception:
                            # If that still fails, raise the original error
                            raise ValueError(f"Failed to parse JSON from Groq response: {content}\n\nError: {str(e)}")
                except Exception as ex:
                    # If that still fails, raise the original error
                    raise ValueError(f"Failed to parse JSON from Groq response: {content}\n\nError: {str(e)}\nSecondary error: {str(ex)}")
            except Exception as e:
                # Print more detailed error information for debugging
                print(f"JSON content that failed validation: {json_content}")
                print(f"Error validating against model: {str(e)}")
                print(f"Model schema: {request.response_model.model_json_schema()}")
                raise ValueError(f"Failed to validate response against model: {str(e)}\n\nResponse: {content}")
                
        except Exception as e:
            raise Exception(f"Groq structured output error: {str(e)}")
        
        # This line should never be reached due to the returns or exceptions above
        # Adding it to satisfy the linter
        raise ValueError("Unexpected error in generate_structured")