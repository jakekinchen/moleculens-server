#!/usr/bin/env python3
"""
Test script for the inject_data_into_template method in MoleculeVisualizer.

This script demonstrates how to use the new method to inject script data and PDB data
into the output.html template file.
"""

import os
import json
from molecule_visualizer import MoleculeVisualizer

# Example script data
script_data = {
    "title": "Methane: The Simplest Hydrocarbon",
    "content": [
        {
            "timecode": "00:00",
            "atoms": [],
            "caption": "Methane is the simplest hydrocarbon with the molecular formula CH₄."
        },
        {
            "timecode": "00:05",
            "atoms": ["0"],
            "caption": "Methane consists of a single carbon atom at the center."
        },
        {
            "timecode": "00:10",
            "atoms": ["1", "2", "3", "4"],
            "caption": "The carbon atom is bonded to four hydrogen atoms in a tetrahedral arrangement."
        },
        {
            "timecode": "00:15",
            "atoms": ["0", "1", "2", "3", "4"],
            "caption": "Methane is the main component of natural gas and a potent greenhouse gas."
        }
    ]
}

# Example PDB data
pdb_data = """COMPND    297
HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  
HETATM    2  H1  UNL     1       0.635   0.635   0.635  1.00  0.00           H  
HETATM    3  H2  UNL     1      -0.635   0.635  -0.635  1.00  0.00           H  
HETATM    4  H3  UNL     1      -0.635  -0.635   0.635  1.00  0.00           H  
HETATM    5  H4  UNL     1       0.635  -0.635  -0.635  1.00  0.00           H  
CONECT    1    2    3    4    5
END"""

def main():
    # Test direct data injection
    print("Testing direct data injection...")
    output_path = MoleculeVisualizer.inject_data_into_template(
        script_data=script_data,
        pdb_data=pdb_data
    )
    print(f"Output saved to: {output_path}")
    
    # Create example files for file-based injection
    current_dir = os.path.dirname(os.path.abspath(__file__))
    script_data_path = os.path.join(current_dir, 'example_script_data.json')
    pdb_data_path = os.path.join(current_dir, 'example_pdb_data.txt')
    
    # Save example script data to file
    with open(script_data_path, 'w') as f:
        json.dump(script_data, f, indent=2)
    
    # Save example PDB data to file
    with open(pdb_data_path, 'w') as f:
        f.write(pdb_data)
    
    print("\nFiles created for testing file-based injection:")
    print(f"- Script data: {script_data_path}")
    print(f"- PDB data: {pdb_data_path}")
    
    print("\nYou can now test file-based injection with:")
    print(f"MoleculeVisualizer.inject_data_into_template(")
    print(f"    script_data_path='{script_data_path}',")
    print(f"    pdb_data_path='{pdb_data_path}'")
    print(f")")

if __name__ == "__main__":
    main() 