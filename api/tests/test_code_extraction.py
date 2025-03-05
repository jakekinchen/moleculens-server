"""
Test code extraction utilities.
"""
import pytest
from agent_management.utils.code_extraction import extract_code_block

def test_code_extraction():
    # Test case 1: JavaScript code block
    input_1 = '''
<thinking>
Let's create a sphere.
</thinking>

```javascript
const geometry = new THREE.SphereGeometry(1, 32, 32);
const material = new THREE.MeshPhongMaterial({ color: 0x0000ff });
const sphere = new THREE.Mesh(geometry, material);
scene.add(sphere);
```
'''
    expected_1 = '''const geometry = new THREE.SphereGeometry(1, 32, 32);
const material = new THREE.MeshPhongMaterial({ color: 0x0000ff });
const sphere = new THREE.Mesh(geometry, material);
scene.add(sphere);'''
    
    assert extract_code_block(input_1, "javascript") == expected_1
    
    # Test case 2: Generic code block
    input_2 = '''
```
const cube = new THREE.BoxGeometry(1, 1, 1);
scene.add(cube);
```
'''
    expected_2 = '''const cube = new THREE.BoxGeometry(1, 1, 1);
scene.add(cube);'''
    
    assert extract_code_block(input_2) == expected_2
    
    # Test case 3: No code block
    input_3 = '''const sphere = new THREE.SphereGeometry(1);
scene.add(sphere);'''
    
    assert extract_code_block(input_3) == input_3
    
    # Test case 4: Multiple thinking sections
    input_4 = '''
<thinking>
First, let's plan the geometry.
</thinking>

```javascript
const geometry = new THREE.BoxGeometry(1, 1, 1);
```

<think>
Now let's add material.
</think>

```javascript
const material = new THREE.MeshBasicMaterial({ color: 0xff0000 });
```
'''
    expected_4 = '''const geometry = new THREE.BoxGeometry(1, 1, 1);'''
    
    assert extract_code_block(input_4, "javascript") == expected_4 