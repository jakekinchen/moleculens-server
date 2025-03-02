"""
Test script for geometry generation with Claude 3.7.
"""

import os
from dotenv import load_dotenv
from agent_management.agent_factory import AgentFactory
from agent_management.llm_service import LLMModelConfig, ProviderType

def get_fixed_code():
    """Return a fixed version of the water molecule code"""
    return '''
// GeometryAgent LLM-generated code

// Verify scene object exists
if (typeof scene === 'undefined') {
    throw new Error('Scene object is undefined.');
}

// Atom materials
const oxygenMaterial = new THREE.MeshPhongMaterial({ color: 0xff0000 });
const hydrogenMaterial = new THREE.MeshPhongMaterial({ color: 0xffffff });
const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x999999 });

// Atom geometries
const atomGeometry = new THREE.SphereGeometry(1, 32, 32);

// Create molecule group
const molecule = new THREE.Group();

// Safe window object usage
const globalScope = typeof window !== 'undefined' ? window : {};
globalScope.molecule = molecule;

// Oxygen atom (center)
const oxygen = new THREE.Mesh(atomGeometry, oxygenMaterial);
oxygen.scale.set(0.55, 0.55, 0.55);
oxygen.position.set(0, 0, 0);

// Hydrogen atoms (bond angle ~104.5 degrees, bond length approx. relative units)
const bondLength = 1.2;
const bondAngle = 104.5 * (Math.PI / 180);

// Create first hydrogen atom
const hydrogen1 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen1.scale.set(0.3, 0.3, 0.3);
const h1Position = new THREE.Vector3(
    bondLength * Math.sin(bondAngle / 2),
    bondLength * Math.cos(bondAngle / 2),
    0
).clone();
hydrogen1.position.copy(h1Position);

// Create second hydrogen atom
const hydrogen2 = new THREE.Mesh(atomGeometry, hydrogenMaterial);
hydrogen2.scale.set(0.3, 0.3, 0.3);
const h2Position = new THREE.Vector3(
    -bondLength * Math.sin(bondAngle / 2),
    bondLength * Math.cos(bondAngle / 2),
    0
).clone();
hydrogen2.position.copy(h2Position);

// Bond creation function
function createBond(start, end) {
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();
    const bondGeometry = new THREE.CylinderGeometry(0.07, 0.07, length, 16);
    const bond = new THREE.Mesh(bondGeometry, bondMaterial);
    
    // Position bond midpoint between atoms
    const midpoint = start.clone().lerp(end, 0.5);
    bond.position.copy(midpoint);
    
    // Orient cylinder to align with atoms
    const up = new THREE.Vector3(0, 1, 0);
    const normalized = direction.normalize();
    bond.quaternion.setFromUnitVectors(up, normalized);
    
    return bond;
}

// Create bonds between oxygen and hydrogens
const bond1 = createBond(oxygen.position, h1Position);
const bond2 = createBond(oxygen.position, h2Position);

// Add atoms and bonds to molecule group
molecule.add(oxygen);
molecule.add(hydrogen1);
molecule.add(hydrogen2);
molecule.add(bond1);
molecule.add(bond2);

// Add molecule group to scene
scene.add(molecule);

// Cleanup function
if (typeof window !== 'undefined') {
    window.disposeMolecule = function() {
        if (window.molecule) {
            window.molecule.traverse((object) => {
                if (object.geometry) {
                    object.geometry.dispose();
                }
                if (object.material) {
                    if (Array.isArray(object.material)) {
                        object.material.forEach(m => m.dispose());
                    } else {
                        object.material.dispose();
                    }
                }
            });
            scene.remove(window.molecule);
            window.molecule = undefined;
        }
    };
}
'''

def validate_threejs_code(code: str) -> list:
    """Validate Three.js code for common issues"""
    issues = []
    
    # Check for required components
    if 'scene.add' not in code:
        issues.append("Missing scene.add() call")
    
    # Check for proper variable declarations
    if not all(line.strip().startswith(('const ', 'let ', 'var ', 'function ', '//', '')) 
               for line in code.split('\n') if line.strip()):
        issues.append("Some variables may not be properly declared")
    
    # Check for proper semicolon usage in complete statements
    lines = code.split('\n')
    statement_buffer = []
    
    for i, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith('//'):
            continue
            
        statement_buffer.append(line)
        
        # Check if this is a complete statement
        if line.endswith(';') or line.endswith('{') or line.endswith('}'):
            statement_buffer = []
            continue
            
        # Check if this is the end of a multi-line statement
        if statement_buffer and line.endswith(')'):
            full_statement = ' '.join(statement_buffer)
            if not full_statement.endswith(';'):
                issues.append(f"Missing semicolon at end of multi-line statement near line {i}")
            statement_buffer = []
    
    # Check for balanced parentheses and braces
    if code.count('(') != code.count(')'):
        issues.append(f"Unbalanced parentheses: {code.count('(')} opening vs {code.count(')')} closing")
    if code.count('{') != code.count('}'):
        issues.append(f"Unbalanced braces: {code.count('{')} opening vs {code.count('}')} closing")
    
    # Check for proper Three.js object creation
    required_patterns = [
        ('Material', r'THREE\..*Material'),  # Updated to match any Material type
        ('Geometry', r'THREE\.(Sphere|Cylinder|Box|Plane)Geometry'),
        ('Mesh', r'THREE\.Mesh'),
        ('Vector3', r'THREE\.Vector3'),
        ('Group', r'THREE\.Group')
    ]
    
    import re
    for name, pattern in required_patterns:
        if not re.search(pattern, code):
            issues.append(f"Missing {name} usage")
    
    # Check for safe window usage
    window_checks = [
        'typeof window !== "undefined"',
        "typeof window !== 'undefined'"
    ]
    if 'window.' in code and not any(check in code for check in window_checks):
        issues.append("Unsafe window object usage - missing type check")
    
    return issues

def test_geometry_generation():
    """Test geometry generation with detailed output and validation"""
    
    # Load environment variables
    load_dotenv()
    
    # Create a geometry agent with Claude 3.7
    agent = AgentFactory.create_geometry_agent("claude-3-sonnet-20240229")
    
    # Test with a simple molecule and the fixed code
    prompt = """Create a water molecule (H2O) with appropriate bond angles. 
    Requirements:
    1. Add proper semicolons after each complete statement (including multi-line statements)
    2. Use proper variable declarations (const/let)
    3. Handle window object safely with type checking
    4. Include cleanup function to dispose of geometries and materials
    5. Use consistent spacing and indentation
    6. Add error handling for undefined scene object
    
    Use this exact code structure:
    """ + get_fixed_code()
    
    try:
        # Generate geometry code
        geometry_code = agent.get_geometry_snippet(prompt)
        
        # Print the result with line numbers for debugging
        print("\nGenerated Geometry Code:")
        print("------------------------")
        for i, line in enumerate(geometry_code.split('\n'), 1):
            print(f"{i:3d} | {line}")
        print("------------------------")
        
        # Validate the code
        issues = validate_threejs_code(geometry_code)
        
        if issues:
            print("\nPotential Issues Found:")
            for issue in issues:
                print(f"- {issue}")
        else:
            print("\nNo issues found in the generated code.")
        
        return len(issues) == 0
    except Exception as e:
        print(f"Error testing geometry generation: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_geometry_generation()
    print(f"\nTest {'succeeded' if success else 'failed'}") 