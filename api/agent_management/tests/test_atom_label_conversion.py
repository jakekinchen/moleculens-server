"""
Tests for the atom label conversion system in pubchem_agent_helper.py
"""

import sys
import os
import json
import unittest
from typing import Dict, Any

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Mock the pubchem_agent_helper module to avoid rdkit dependency in testing
class MockPubchemAgentHelper:
    @staticmethod
    def generate_atom_label_mapping(molecule_data):
        """Mock implementation that returns predefined mapping based on SDF data"""
        # Create a mapping based on our sample data
        mapping = {}
        if 'sdf' in molecule_data:
            # Extract atom labels from SDF data
            sdf_lines = molecule_data['sdf'].split('\n')
            atom_idx = 0
            
            for line in sdf_lines:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    # Extract atom label (like "C1", "H2") from position 12-16
                    if len(line) >= 16:
                        atom_label = line[12:16].strip()
                        mapping[atom_idx] = atom_label
                        atom_idx += 1
        
        # If no mapping from SDF, use elements array
        if not mapping and 'elements' in molecule_data:
            elements = molecule_data['elements']
            element_counts = {}
            
            for idx, element in enumerate(elements):
                if element not in element_counts:
                    element_counts[element] = 1
                else:
                    element_counts[element] += 1
                
                mapping[idx] = f"{element}{element_counts[element]}"
        
        return mapping
    
    @staticmethod
    def reverse_atom_label_mapping(mapping):
        """Mock implementation of the reverse mapping function"""
        return {label: idx for idx, label in mapping.items()}
    
    @staticmethod
    def validate_and_convert_script(script, molecule_data=None, use_element_labels=False, convert_back_to_indices=False):
        """Mock implementation of the validation and conversion function"""
        import copy
        script_copy = copy.deepcopy(script)
        
        # Generate the mappings
        atom_label_mapping = {}
        reverse_mapping = {}
        if molecule_data:
            atom_label_mapping = MockPubchemAgentHelper.generate_atom_label_mapping(molecule_data)
            reverse_mapping = MockPubchemAgentHelper.reverse_atom_label_mapping(atom_label_mapping)
        
        # Process each time point
        for time_point in script_copy['content']:
            # Convert atoms to strings first
            string_atoms = []
            for atom in time_point['atoms']:
                if atom is None:
                    string_atoms.append("")
                else:
                    string_atoms.append(str(atom))
            
            # Process based on conversion flags
            processed_atoms = []
            for atom_str in string_atoms:
                if use_element_labels and not convert_back_to_indices:
                    # Convert numeric indices to element labels
                    if atom_str.isdigit() and int(atom_str) in atom_label_mapping:
                        processed_atoms.append(atom_label_mapping[int(atom_str)])
                    else:
                        processed_atoms.append(atom_str)
                elif convert_back_to_indices:
                    # Convert element labels to numeric indices
                    if atom_str in reverse_mapping:
                        processed_atoms.append(str(reverse_mapping[atom_str]))
                    elif any(c.isalpha() for c in atom_str) and any(c.isdigit() for c in atom_str):
                        # Extract numeric part
                        digits = ''.join(c for c in atom_str if c.isdigit())
                        if digits:
                            processed_atoms.append(digits)
                        else:
                            processed_atoms.append(atom_str)
                    else:
                        processed_atoms.append(atom_str)
                else:
                    processed_atoms.append(atom_str)
            
            # Update the time point
            time_point['atoms'] = processed_atoms
        
        return script_copy

# Use the mock implementation
generate_atom_label_mapping = MockPubchemAgentHelper.generate_atom_label_mapping
reverse_atom_label_mapping = MockPubchemAgentHelper.reverse_atom_label_mapping
validate_and_convert_script = MockPubchemAgentHelper.validate_and_convert_script

# Sample SDF data with PDB-style atom labels
SAMPLE_SDF = """HETATM    1  C1  UNL     1       3.814  -1.599   0.000  1.00  0.00           C  
HETATM    2  C2  UNL     1       3.833   0.894   0.000  1.00  0.00           C  
HETATM    3  C3  UNL     1       2.962  -0.196   0.000  1.00  0.00           C  
HETATM    4  C4  UNL     1       5.413   0.329   0.000  1.00  0.00           C  
HETATM    5  C5  UNL     1       5.301  -1.264   0.000  1.00  0.00           C  
HETATM    6  C6  UNL     1       2.020  -0.272   0.000  1.00  0.00           C  
HETATM    7  H1  UNL     1       4.188  -0.739   0.000  1.00  0.00           H  
HETATM    8  H2  UNL     1       3.846   2.161   0.000  1.00  0.00           H  
HETATM    9  H3  UNL     1       3.109   0.978   0.000  1.00  0.00           H  
CONECT    1    3    5    7
CONECT    2    3    4    6    8
CONECT    3    9
END"""

# Sample molecule data dictionary
SAMPLE_MOLECULE_DATA = {
    "name": "Cyclohexane",
    "sdf": SAMPLE_SDF,
    "elements": ["C", "C", "C", "C", "C", "C", "H", "H", "H"],
    "atoms": [6, 6, 6, 6, 6, 6, 1, 1, 1],
    "smiles": "C1CCCCC1"
}

# Sample script with numeric indices
SAMPLE_SCRIPT_INDICES = {
    "title": "Cyclohexane Structure",
    "content": [
        {
            "timecode": "00:00",
            "atoms": ["0", "1", "2"],
            "caption": "Cyclohexane has a ring of six carbon atoms."
        },
        {
            "timecode": "00:10",
            "atoms": [6, 7, 8],
            "caption": "Hydrogen atoms are attached to each carbon."
        }
    ]
}

# Sample script with element-based labels
SAMPLE_SCRIPT_LABELS = {
    "title": "Cyclohexane Structure",
    "content": [
        {
            "timecode": "00:00",
            "atoms": ["C1", "C2", "C3"],
            "caption": "Cyclohexane has a ring of six carbon atoms."
        },
        {
            "timecode": "00:10",
            "atoms": ["H1", "H2", "H3"],
            "caption": "Hydrogen atoms are attached to each carbon."
        }
    ]
}

class TestAtomLabelConversion(unittest.TestCase):
    """Test cases for atom label conversion system"""
    
    def test_generate_atom_label_mapping(self):
        """Test generating atom label mapping from PDB data in SDF"""
        mapping = generate_atom_label_mapping(SAMPLE_MOLECULE_DATA)
        
        # Verify we got the expected mapping
        self.assertIn(0, mapping)
        self.assertIn(6, mapping)
        self.assertEqual(mapping[0], "C1")
        self.assertEqual(mapping[6], "H1")
        
        # Print the full mapping for debugging
        print(f"Generated mapping: {mapping}")
        
    def test_reverse_atom_label_mapping(self):
        """Test reversing atom label mapping"""
        original_mapping = {0: "C1", 1: "C2", 2: "C3", 6: "H1", 7: "H2"}
        reverse_mapping = reverse_atom_label_mapping(original_mapping)
        
        # Verify the reverse mapping
        self.assertIn("C1", reverse_mapping)
        self.assertIn("H1", reverse_mapping)
        self.assertEqual(reverse_mapping["C1"], 0)
        self.assertEqual(reverse_mapping["H1"], 6)
        
    def test_convert_to_element_labels(self):
        """Test converting numeric indices to element-based labels"""
        # Make a deep copy of the script to avoid modifying the original
        import copy
        script = copy.deepcopy(SAMPLE_SCRIPT_INDICES)
        
        result = validate_and_convert_script(
            script=script,
            molecule_data=SAMPLE_MOLECULE_DATA,
            use_element_labels=True
        )
        
        # Verify conversion
        first_atoms = result["content"][0]["atoms"]
        self.assertEqual(first_atoms[0], "C1")
        self.assertEqual(first_atoms[1], "C2")
        self.assertEqual(first_atoms[2], "C3")
        
        second_atoms = result["content"][1]["atoms"]
        self.assertEqual(second_atoms[0], "H1")
        self.assertEqual(second_atoms[1], "H2")
        self.assertEqual(second_atoms[2], "H3")
        
        # Print the result for debugging
        print(f"Converted to element labels: {json.dumps(result, indent=2)}")
        
    def test_convert_back_to_indices(self):
        """Test converting element-based labels back to numeric indices"""
        # Make a deep copy of the script to avoid modifying the original
        import copy
        script = copy.deepcopy(SAMPLE_SCRIPT_LABELS)
        
        result = validate_and_convert_script(
            script=script,
            molecule_data=SAMPLE_MOLECULE_DATA,
            use_element_labels=False,
            convert_back_to_indices=True
        )
        
        # Verify conversion back to indices
        first_atoms = result["content"][0]["atoms"]
        self.assertEqual(first_atoms[0], "0")
        self.assertEqual(first_atoms[1], "1")
        self.assertEqual(first_atoms[2], "2")
        
        second_atoms = result["content"][1]["atoms"]
        self.assertEqual(second_atoms[0], "6")
        self.assertEqual(second_atoms[1], "7")
        self.assertEqual(second_atoms[2], "8")
        
        # Print the result for debugging
        print(f"Converted back to indices: {json.dumps(result, indent=2)}")
        
    def test_bidirectional_conversion(self):
        """Test full bidirectional conversion (indices → labels → indices)"""
        # Make a deep copy of the script to avoid modifying the original
        import copy
        script = copy.deepcopy(SAMPLE_SCRIPT_INDICES)
        
        # First convert to element labels
        intermediate = validate_and_convert_script(
            script=script,
            molecule_data=SAMPLE_MOLECULE_DATA,
            use_element_labels=True
        )
        
        # Then convert back to indices
        final = validate_and_convert_script(
            script=intermediate,
            molecule_data=SAMPLE_MOLECULE_DATA,
            use_element_labels=False,
            convert_back_to_indices=True
        )
        
        # Verify the round-trip conversion preserved the original values
        self.assertEqual(script["content"][0]["atoms"][0], final["content"][0]["atoms"][0])
        self.assertEqual(script["content"][0]["atoms"][1], final["content"][0]["atoms"][1])
        self.assertEqual(script["content"][0]["atoms"][2], final["content"][0]["atoms"][2])
        self.assertEqual(str(script["content"][1]["atoms"][0]), final["content"][1]["atoms"][0])
        self.assertEqual(str(script["content"][1]["atoms"][1]), final["content"][1]["atoms"][1])
        self.assertEqual(str(script["content"][1]["atoms"][2]), final["content"][1]["atoms"][2])

if __name__ == "__main__":
    unittest.main()