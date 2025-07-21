"""
Helper functions for the PubChem Agent to handle molecule data validation and conversion.
"""

from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem


def generate_atom_label_mapping(molecule_data: Dict[str, Any]) -> Dict[int, str]:
    """
    Generate a mapping from atom indices to element-based labels (e.g., C1, C2, O1).
    Uses PDB data as primary source, falls back to SMILES or element arrays.

    Args:
        molecule_data: Dictionary containing molecule information

    Returns:
        Dict[int, str]: Mapping from atom indices to element-based labels
    """
    # Default empty mapping
    mapping = {}

    # APPROACH 1: Use PDB data from SDF (most accurate)
    if molecule_data.get("sdf"):
        try:
            pdb_lines = molecule_data["sdf"].split("\n")
            atom_idx = 0

            for line in pdb_lines:
                # Match HETATM or ATOM lines that contain atom definitions in PDB format
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    # Extract the atom label (like "C1", "H2") from position 12-16
                    if len(line) >= 16:
                        atom_label = line[12:16].strip()
                        mapping[atom_idx] = atom_label
                        atom_idx += 1

            # If we found atom labels, return this mapping
            if mapping:
                return mapping
        except Exception as e:
            print(f"Error parsing PDB data: {str(e)}")
            # Fall through to alternative methods

    # APPROACH 2: Use elements array directly
    if molecule_data.get("elements") and isinstance(molecule_data["elements"], list):
        try:
            elements = molecule_data["elements"]
            # Count elements to create proper numbering
            element_counts = {}

            for idx, element in enumerate(elements):
                if element not in element_counts:
                    element_counts[element] = 1
                else:
                    element_counts[element] += 1

                # Create label like "C1", "O2", etc.
                label = f"{element}{element_counts[element]}"
                mapping[idx] = label

            # If we found elements, return this mapping
            if mapping:
                return mapping
        except Exception as e:
            print(f"Error using elements array: {str(e)}")
            # Fall through to SMILES method

    # APPROACH 3: Use SMILES string (fallback method)
    if molecule_data.get("smiles"):
        try:
            # Parse the SMILES string to get atom information
            mol = Chem.MolFromSmiles(molecule_data["smiles"])
            if mol:
                # Count atoms of each element type
                element_counts = {}

                # Create mapping from atom index to element-based label
                for atom_idx in range(mol.GetNumAtoms()):
                    atom = mol.GetAtomWithIdx(atom_idx)
                    element = atom.GetSymbol()

                    # Increment count for this element
                    if element not in element_counts:
                        element_counts[element] = 1
                    else:
                        element_counts[element] += 1

                    # Create label like "C1", "O2", etc.
                    label = f"{element}{element_counts[element]}"
                    mapping[atom_idx] = label
        except Exception as e:
            print(f"Error using SMILES method: {str(e)}")

    return mapping


def reverse_atom_label_mapping(atom_label_mapping: Dict[int, str]) -> Dict[str, int]:
    """
    Create a reverse mapping from element-based labels back to numeric indices.

    Args:
        atom_label_mapping: Dictionary mapping atom indices to element-based labels

    Returns:
        Dict[str, int]: Mapping from element-based labels to atom indices
    """
    return {label: idx for idx, label in atom_label_mapping.items()}


def validate_and_convert_script(
    script: Dict[str, Any],
    molecule_data: Optional[Dict[str, Any]] = None,
    use_element_labels: bool = False,
    convert_back_to_indices: bool = False,
) -> Dict[str, Any]:
    """
    Validates the script structure and ensures all atoms are strings.
    Can convert between numeric indices and element-based labels.

    Args:
        script (Dict[str, Any]): The script data to validate
        molecule_data (Optional[Dict[str, Any]]): Molecule data for element-based label conversion
        use_element_labels (bool): Whether to use element-based labels (C1, O1) instead of indices
        convert_back_to_indices (bool): After script agent returns, convert element-labels back to numeric indices

    Returns:
        Dict[str, Any]: The validated script with atoms converted to strings if needed
    """
    # Basic structure validation
    if not isinstance(script, dict):
        raise ValueError("Script must be a dictionary")
    if "title" not in script:
        raise ValueError("Script must have a 'title' field")
    if "content" not in script:
        raise ValueError("Script must have a 'content' field")
    if not isinstance(script["content"], list):
        raise ValueError("Script content must be a list")

    # Deep copy the script to avoid modifying the original during processing
    import copy

    script = copy.deepcopy(script)

    # Generate element-based atom labels if requested
    atom_label_mapping = {}
    reverse_mapping = {}
    if molecule_data:
        atom_label_mapping = generate_atom_label_mapping(molecule_data)
        reverse_mapping = reverse_atom_label_mapping(atom_label_mapping)

    # Process each time point
    for i, time_point in enumerate(script["content"]):
        if not isinstance(time_point, dict):
            raise ValueError(f"Time point at index {i} must be a dictionary")
        if "timecode" not in time_point:
            raise ValueError(f"Time point at index {i} must have a 'timecode' field")
        if "atoms" not in time_point:
            raise ValueError(f"Time point at index {i} must have an 'atoms' field")
        if "caption" not in time_point:
            raise ValueError(f"Time point at index {i} must have a 'caption' field")

        # Handle atoms field
        atoms = time_point["atoms"]

        # If atoms is not a list, convert it to a single-item list
        if not isinstance(atoms, list):
            atoms = [atoms]

        processed_atoms = []

        # PHASE 1: Convert everything to valid string format first
        string_atoms = []
        for atom in atoms:
            if atom is None:
                string_atoms.append("")  # Handle None values
            elif isinstance(atom, (int, float)):
                string_atoms.append(
                    str(int(atom))
                )  # Convert numbers to string integers
            elif isinstance(atom, str):
                string_atoms.append(atom)  # Keep strings as is
            else:
                string_atoms.append(str(atom))  # Convert anything else to string

        # PHASE 2: Apply the appropriate conversion based on settings
        for atom_str in string_atoms:
            # Case 1: Need element-labels for script agent
            if use_element_labels and not convert_back_to_indices:
                # Convert numeric string to element-label if possible
                if atom_str.isdigit() and int(atom_str) in atom_label_mapping:
                    processed_atoms.append(atom_label_mapping[int(atom_str)])
                # Otherwise keep as is (might already be element-label)
                else:
                    processed_atoms.append(atom_str)

            # Case 2: Need to convert element-labels back to indices
            elif convert_back_to_indices:
                # If it's an element-label, convert to numeric index
                if atom_str in reverse_mapping:
                    processed_atoms.append(str(reverse_mapping[atom_str]))
                # If it looks like an element-label but not in mapping, try to keep numeric part
                elif any(c.isalpha() for c in atom_str) and any(
                    c.isdigit() for c in atom_str
                ):
                    # Try to extract numeric part as a fallback
                    try:
                        # Extract digits from something like "C1" → "1"
                        digits = "".join(c for c in atom_str if c.isdigit())
                        if digits:
                            processed_atoms.append(digits)
                        else:
                            processed_atoms.append(atom_str)
                    except:
                        processed_atoms.append(atom_str)
                # Otherwise keep as is (likely already a numeric index)
                else:
                    processed_atoms.append(atom_str)

            # Case 3: Just pass through with standard string validation
            else:
                processed_atoms.append(atom_str)

        # Replace the atoms field with the processed values
        time_point["atoms"] = processed_atoms

    return script
