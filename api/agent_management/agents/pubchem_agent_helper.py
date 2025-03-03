"""
Helper functions for the PubChem Agent to handle molecule data validation and conversion.
"""

from typing import Dict, Any, List

def validate_and_convert_script(script: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validates the script structure and ensures all atoms are strings.
    
    Args:
        script (Dict[str, Any]): The script data to validate
        
    Returns:
        Dict[str, Any]: The validated script with atoms converted to strings if needed
    """
    # Basic structure validation
    if not isinstance(script, dict):
        raise ValueError("Script must be a dictionary")
    if 'title' not in script:
        raise ValueError("Script must have a 'title' field")
    if 'content' not in script:
        raise ValueError("Script must have a 'content' field")
    if not isinstance(script['content'], list):
        raise ValueError("Script content must be a list")
    
    # Process each time point
    for i, time_point in enumerate(script['content']):
        if not isinstance(time_point, dict):
            raise ValueError(f"Time point at index {i} must be a dictionary")
        if 'timecode' not in time_point:
            raise ValueError(f"Time point at index {i} must have a 'timecode' field")
        if 'atoms' not in time_point:
            raise ValueError(f"Time point at index {i} must have an 'atoms' field")
        if 'caption' not in time_point:
            raise ValueError(f"Time point at index {i} must have a 'caption' field")
        
        # Convert atoms to strings if they are integers
        atoms = time_point['atoms']
        if not isinstance(atoms, list):
            atoms = [str(atoms)]
            time_point['atoms'] = atoms
        
        # Convert each atom to string if needed
        string_atoms = []
        for j, atom in enumerate(atoms):
            if not isinstance(atom, str):
                string_atoms.append(str(atom))
            else:
                string_atoms.append(atom)
        time_point['atoms'] = string_atoms
    
    return script