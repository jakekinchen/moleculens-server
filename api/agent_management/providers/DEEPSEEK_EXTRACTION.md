# DeepSeek Model Response Extraction

## Overview

DeepSeek models have a specific response format that requires special handling. When generating structured output (like JavaScript code or JSON), DeepSeek models typically return responses with the following structure:

1. A thinking section enclosed in either `<thinking>...</thinking>` or `<think>...</think>` XML tags containing the model's reasoning process
2. The actual structured output (e.g., JavaScript code or JSON) in a code block (```javascript ... ```, ```html ... ```, or ```json ... ```)

This document explains the implementation of special extraction functions for DeepSeek model responses.

## Implementation Details

### Utility Functions

The following utility functions have been implemented in `api/agent_management/providers/deepseek_utils.py`:

1. `is_deepseek_model(model_name: str) -> bool`
   - Checks if a model is a DeepSeek model based on its name
   - Returns `True` if the model name contains "deepseek" (case-insensitive)

2. `extract_javascript_from_deepseek_response(response: str) -> str`
   - Extracts JavaScript code from a DeepSeek model response
   - First tries to find a JavaScript code block (```javascript ... ```)
   - Also tries to find an HTML code block (```html ... ```) which may contain JavaScript
   - If neither is found, tries to find any code block (``` ... ```)
   - If no code block is found, returns the original response

3. `extract_structured_output_from_deepseek(response: str, output_type: str = "javascript") -> str`
   - Extracts structured output from a DeepSeek model response based on the output type
   - First removes any thinking section enclosed in either `<thinking>...</thinking>` or `<think>...</think>` tags
   - Then extracts the structured output based on the specified output type (e.g., "javascript" or "json")
   - If no matching block is found, returns the original response

### Integration with Agents

The following agents have been updated to use the DeepSeek extraction functions:

1. **GeometryAgent**
   - Checks if the model is a DeepSeek model using `is_deepseek_model`
   - If it is, applies `extract_structured_output_from_deepseek` to extract JavaScript code

2. **CaptionAgent**
   - Checks if the model is a DeepSeek model using `is_deepseek_model`
   - If it is, applies `extract_structured_output_from_deepseek` to extract JSON data

### Integration with Providers

The **GroqProvider** has been updated to use the DeepSeek extraction functions for structured output:

- In the `generate_structured` method, it checks if the model is a DeepSeek model using `is_deepseek_model`
- If it is, applies `extract_structured_output_from_deepseek` to extract JSON data before further processing

## Testing

A test script has been created at `api/test_deepseek_extraction.py` to verify that the DeepSeek extraction functions work correctly. The script includes:

1. Tests for the `is_deepseek_model` function
2. Tests for the `extract_javascript_from_deepseek_response` function, including HTML code blocks
3. Tests for the `extract_structured_output_from_deepseek` function, including both `<thinking>` and `<think>` tags
4. An optional test with an actual DeepSeek model (requires the `GROQ_API_KEY` environment variable to be set)

To run the tests:

```bash
python api/test_deepseek_extraction.py
```

To include the live model test:

```bash
python api/test_deepseek_extraction.py --with-model
```

## Example Response Format

Here are examples of typical DeepSeek model responses:

### Example 1: Using `<thinking>` tags

```
<thinking>
I need to create a simple sphere in Three.js.
</thinking>

Here's the code to create a sphere in Three.js:

```javascript
// Create a sphere
const geometry = new THREE.SphereGeometry(1, 32, 32);
const material = new THREE.MeshPhongMaterial({ color: 0x0000ff });
const sphere = new THREE.Mesh(geometry, material);
sphere.position.set(0, 0, 0);
scene.add(sphere);
```
```

### Example 2: Using `<think>` tags

```
<think>
I'll create a simple animation.
</think>

```html
<!DOCTYPE html>
<html>
<head>
    <title>Animation</title>
</head>
<body>
    <script>
        const cube = new THREE.Mesh(geometry, material);
        scene.add(cube);
        
        function animate() {
            requestAnimationFrame(animate);
            cube.rotation.y += 0.01;
            renderer.render(scene, camera);
        }
        
        animate();
    </script>
</body>
</html>
```
```

After extraction, the result would be the code without the thinking section and code block markers.

## Future Improvements

1. Add support for more output types (e.g., XML, CSV)
2. Implement model-specific extraction functions for other models that may have unique response formats
3. Add more robust error handling for edge cases
4. Consider adding a configuration option to enable/disable extraction for specific models 