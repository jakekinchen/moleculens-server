"""
Demonstration script for the atom label conversion system.
Shows how the conversion works in a practical example.
"""

import sys
import os
import json
from pprint import pprint

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    # Try to import the atom label conversion functions
    from agents.pubchem_agent_helper import (
        generate_atom_label_mapping,
        reverse_atom_label_mapping,
        validate_and_convert_script
    )
    from agent_model_config import create_agent_llm_service, AgentType
    from agent_factory import AgentFactory
    
    # Create a sample water molecule
    WATER_MOLECULE = {
        "name": "Water",
        "smiles": "O",
        "elements": ["O", "H", "H"],
        "sdf": """HETATM    1  O   HOH     1       0.000   0.000   0.000  1.00  0.00           O  
HETATM    2  H1  HOH     1       0.000   0.000   1.000  1.00  0.00           H  
HETATM    3  H2  HOH     1       0.943   0.000  -0.300  1.00  0.00           H  
CONECT    1    2    3
END"""
    }
    
    # Sample script with numeric indices
    SCRIPT_WITH_INDICES = {
        "title": "Water Molecule Structure",
        "content": [
            {
                "timecode": "00:00",
                "atoms": ["0"],
                "caption": "This is the oxygen atom at the center of the water molecule."
            },
            {
                "timecode": "00:05",
                "atoms": ["1", "2"],
                "caption": "These are the two hydrogen atoms bonded to the oxygen."
            },
            {
                "timecode": "00:10", 
                "atoms": ["0", "1", "2"],
                "caption": "Together, they form the water molecule with a bent structure."
            }
        ]
    }
    
    def demo_atom_label_conversion():
        """Demonstrate the atom label conversion system"""
        print("=" * 60)
        print("ATOM LABEL CONVERSION DEMONSTRATION")
        print("=" * 60)
        
        # Step 1: Generate mapping from indices to element labels
        print("\n1. GENERATING ELEMENT LABEL MAPPING")
        print("-" * 60)
        mapping = generate_atom_label_mapping(WATER_MOLECULE)
        print("Mapping from indices to element labels:")
        pprint(mapping)
        
        # Step 2: Generate reverse mapping
        print("\n2. GENERATING REVERSE MAPPING")
        print("-" * 60)
        reverse_mapping = reverse_atom_label_mapping(mapping)
        print("Mapping from element labels to indices:")
        pprint(reverse_mapping)
        
        # Step 3: Convert script from indices to element labels
        print("\n3. CONVERTING SCRIPT FROM INDICES TO ELEMENT LABELS")
        print("-" * 60)
        print("Original script with indices:")
        pprint(SCRIPT_WITH_INDICES)
        
        element_script = validate_and_convert_script(
            script=SCRIPT_WITH_INDICES,
            molecule_data=WATER_MOLECULE,
            use_element_labels=True
        )
        
        print("\nConverted script with element labels:")
        pprint(element_script)
        
        # Step 4: Convert back to indices
        print("\n4. CONVERTING SCRIPT FROM ELEMENT LABELS BACK TO INDICES")
        print("-" * 60)
        indices_script = validate_and_convert_script(
            script=element_script,
            molecule_data=WATER_MOLECULE,
            use_element_labels=False,
            convert_back_to_indices=True
        )
        
        print("Final script converted back to indices:")
        pprint(indices_script)
        
        # Step 5: Demonstrate configuration via AgentFactory
        print("\n5. CONFIGURING PUBCHEM AGENT VIA FACTORY")
        print("-" * 60)
        print("Agent with element labels enabled (better for LLM understanding):")
        print("agent = AgentFactory.create_pubchem_agent(")
        print("    use_element_labels=True,")
        print("    convert_back_to_indices=False")
        print(")")
        
        print("\nAgent with conversion back to indices (for visualization):")
        print("agent = AgentFactory.create_pubchem_agent(")
        print("    use_element_labels=True,")
        print("    convert_back_to_indices=True")
        print(")")
        
        print("\n" + "=" * 60)
        print("DEMONSTRATION COMPLETE")
        print("=" * 60)
    
    if __name__ == "__main__":
        demo_atom_label_conversion()
        
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("This demo requires the full environment with RDKit and other dependencies.")
    print("Please check the README_ATOM_LABELS.md file for more information.")