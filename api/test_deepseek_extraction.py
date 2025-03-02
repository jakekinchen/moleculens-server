#!/usr/bin/env python3
"""
Test script to verify that the DeepSeek extraction functions work correctly.
"""

import os
import sys
from typing import Dict, Any, List

# Add the parent directory to the path so we can import the modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.agent_management.providers.deepseek_utils import (
    is_deepseek_model,
    extract_javascript_from_deepseek_response,
    extract_structured_output_from_deepseek
)
from api.agent_management.model_config import get_llm_service
from api.agent_management.llm_service import LLMRequest

def test_is_deepseek_model():
    """Test the is_deepseek_model function."""
    print("\n=== Testing is_deepseek_model ===")
    
    # Test cases
    test_cases = [
        ("deepseek-r1-distill-llama-70b", True),
        ("DeepSeek-Coder", True),
        ("llama-3.2-90b-vision-preview", False),
        ("gpt-4o", False),
        ("claude-3-7-sonnet-latest", False)
    ]
    
    for model_name, expected in test_cases:
        result = is_deepseek_model(model_name)
        print(f"Model: {model_name}, Expected: {expected}, Result: {result}")
        assert result == expected, f"Expected {expected} for {model_name}, got {result}"

def test_extract_javascript():
    """Test the extract_javascript_from_deepseek_response function."""
    print("\n=== Testing extract_javascript_from_deepseek_response ===")
    
    # Test case 1: JavaScript code block
    test_case_1 = """<thinking>
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
```"""
    
    expected_1 = """// Create a sphere
const geometry = new THREE.SphereGeometry(1, 32, 32);
const material = new THREE.MeshPhongMaterial({ color: 0x0000ff });
const sphere = new THREE.Mesh(geometry, material);
sphere.position.set(0, 0, 0);
scene.add(sphere);"""
    
    result_1 = extract_javascript_from_deepseek_response(test_case_1)
    print(f"Test case 1 result matches expected: {result_1 == expected_1}")
    assert result_1 == expected_1, f"Expected:\n{expected_1}\n\nGot:\n{result_1}"
    
    # Test case 2: Generic code block
    test_case_2 = """<thinking>
I'll create a simple cube.
</thinking>

```
// Create a cube
const geometry = new THREE.BoxGeometry(1, 1, 1);
const material = new THREE.MeshPhongMaterial({ color: 0xff0000 });
const cube = new THREE.Mesh(geometry, material);
cube.position.set(0, 0, 0);
scene.add(cube);
```"""
    
    expected_2 = """// Create a cube
const geometry = new THREE.BoxGeometry(1, 1, 1);
const material = new THREE.MeshPhongMaterial({ color: 0xff0000 });
const cube = new THREE.Mesh(geometry, material);
cube.position.set(0, 0, 0);
scene.add(cube);"""
    
    result_2 = extract_javascript_from_deepseek_response(test_case_2)
    print(f"Test case 2 result matches expected: {result_2 == expected_2}")
    assert result_2 == expected_2, f"Expected:\n{expected_2}\n\nGot:\n{result_2}"
    
    # Test case 3: No code block
    test_case_3 = """// Create a cylinder
const geometry = new THREE.CylinderGeometry(0.5, 0.5, 2, 32);
const material = new THREE.MeshPhongMaterial({ color: 0x00ff00 });
const cylinder = new THREE.Mesh(geometry, material);
cylinder.position.set(0, 0, 0);
scene.add(cylinder);"""
    
    expected_3 = test_case_3
    
    result_3 = extract_javascript_from_deepseek_response(test_case_3)
    print(f"Test case 3 result matches expected: {result_3 == expected_3}")
    assert result_3 == expected_3, f"Expected:\n{expected_3}\n\nGot:\n{result_3}"
    
    # Test case 4: HTML block with JavaScript
    test_case_4 = """<think>
I'll create an HTML page with JavaScript.
</think>

```html
<!DOCTYPE html>
<html>
<head>
    <title>Test</title>
</head>
<body>
    <script>
        const cube = new THREE.Mesh(geometry, material);
        scene.add(cube);
    </script>
</body>
</html>
```"""
    
    expected_4 = """<!DOCTYPE html>
<html>
<head>
    <title>Test</title>
</head>
<body>
    <script>
        const cube = new THREE.Mesh(geometry, material);
        scene.add(cube);
    </script>
</body>
</html>"""
    
    result_4 = extract_javascript_from_deepseek_response(test_case_4)
    print(f"Test case 4 result matches expected: {result_4 == expected_4}")
    assert result_4 == expected_4, f"Expected:\n{expected_4}\n\nGot:\n{result_4}"

def test_extract_structured_output():
    """Test the extract_structured_output_from_deepseek function."""
    print("\n=== Testing extract_structured_output_from_deepseek ===")
    
    # Test case 1: JSON with thinking section
    test_case_1 = """<thinking>
I need to create a JSON array of captions.
</thinking>

```json
[
  { "time": 0, "text": "The beginning of our journey" },
  { "time": 5, "text": "Exploring the molecular structure" },
  { "time": 10, "text": "Observing the chemical reaction" }
]
```"""
    
    expected_1 = """[
  { "time": 0, "text": "The beginning of our journey" },
  { "time": 5, "text": "Exploring the molecular structure" },
  { "time": 10, "text": "Observing the chemical reaction" }
]"""
    
    result_1 = extract_structured_output_from_deepseek(test_case_1, "json")
    print(f"Test case 1 result matches expected: {result_1 == expected_1}")
    assert result_1 == expected_1, f"Expected:\n{expected_1}\n\nGot:\n{result_1}"
    
    # Test case 2: JavaScript with thinking section
    test_case_2 = """<thinking>
I'll create a simple animation.
</thinking>

```javascript
// Create an animation
const box = new THREE.BoxGeometry(1, 1, 1);
const material = new THREE.MeshBasicMaterial({ color: 0xffffff });
const cube = new THREE.Mesh(box, material);
scene.add(cube);

function animate() {
  cube.rotation.x += 0.01;
  cube.rotation.y += 0.01;
}
```"""
    
    expected_2 = """// Create an animation
const box = new THREE.BoxGeometry(1, 1, 1);
const material = new THREE.MeshBasicMaterial({ color: 0xffffff });
const cube = new THREE.Mesh(box, material);
scene.add(cube);

function animate() {
  cube.rotation.x += 0.01;
  cube.rotation.y += 0.01;
}"""
    
    result_2 = extract_structured_output_from_deepseek(test_case_2, "javascript")
    print(f"Test case 2 result matches expected: {result_2 == expected_2}")
    assert result_2 == expected_2, f"Expected:\n{expected_2}\n\nGot:\n{result_2}"
    
    # Test case 3: JavaScript with think section (alternative tag)
    test_case_3 = """<think>
I'll create a simple animation with a different tag.
</think>

```javascript
// Create an animation with think tag
const sphere = new THREE.SphereGeometry(1, 32, 32);
const material = new THREE.MeshBasicMaterial({ color: 0x00ff00 });
const ball = new THREE.Mesh(sphere, material);
scene.add(ball);
```"""
    
    expected_3 = """// Create an animation with think tag
const sphere = new THREE.SphereGeometry(1, 32, 32);
const material = new THREE.MeshBasicMaterial({ color: 0x00ff00 });
const ball = new THREE.Mesh(sphere, material);
scene.add(ball);"""
    
    result_3 = extract_structured_output_from_deepseek(test_case_3, "javascript")
    print(f"Test case 3 result matches expected: {result_3 == expected_3}")
    assert result_3 == expected_3, f"Expected:\n{expected_3}\n\nGot:\n{result_3}"

def test_with_deepseek_model():
    """Test with an actual DeepSeek model if available."""
    print("\n=== Testing with DeepSeek model (if available) ===")
    
    # Check if the DeepSeek model is available
    try:
        # Check if GROQ_API_KEY is set
        if not os.getenv("GROQ_API_KEY"):
            print("GROQ_API_KEY environment variable not set, skipping live model test")
            return
            
        # Try to get the DeepSeek model service
        llm_service = get_llm_service("deepseek-r1-distill-llama-70b")
        
        # Create a simple prompt
        prompt = "Create a simple Three.js scene with a rotating red cube"
        
        # Create a request
        request = LLMRequest(user_prompt=prompt)
        
        # Generate a response
        print("Generating response from DeepSeek model...")
        response = llm_service.generate(request)
        
        # Print the raw response
        print("\nRaw response from DeepSeek model:")
        print("=" * 80)
        print(response.content)
        print("=" * 80)
        
        # Extract the JavaScript code
        extracted_js = extract_structured_output_from_deepseek(response.content, "javascript")
        
        # Print the extracted JavaScript code
        print("\nExtracted JavaScript code:")
        print("=" * 80)
        print(extracted_js)
        print("=" * 80)
        
        # Check if the extraction worked
        if ("<thinking>" in response.content and "<thinking>" not in extracted_js) or \
           ("<think>" in response.content and "<think>" not in extracted_js):
            print("Successfully extracted JavaScript code from DeepSeek response!")
        else:
            print("Note: No thinking tags found in the response or extraction failed")
            
    except Exception as e:
        print(f"Error testing with DeepSeek model: {str(e)}")

if __name__ == "__main__":
    # Run the tests
    test_is_deepseek_model()
    test_extract_javascript()
    test_extract_structured_output()
    
    # Only run the live model test if explicitly requested
    if len(sys.argv) > 1 and sys.argv[1] == "--with-model":
        test_with_deepseek_model()
    else:
        print("\nSkipping live model test. Run with --with-model to include it.")
    
    print("\nAll tests passed!") 