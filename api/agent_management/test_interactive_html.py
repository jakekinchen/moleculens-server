#!/usr/bin/env python3
"""
Test script for the generate_interactive_html method in MoleculeVisualizer.

This script demonstrates how to use the generate_interactive_html method
to create an interactive HTML visualization with custom script and PDB data.
"""

import os
import json
from typing import Dict, Any, List
from molecule_visualizer import MoleculeVisualizer

# Example PDB data
pdb_data = """COMPND    297
HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  
HETATM    2  H1  UNL     1       0.635   0.635   0.635  1.00  0.00           H  
HETATM    3  H2  UNL     1      -0.635   0.635  -0.635  1.00  0.00           H  
HETATM    4  H3  UNL     1      -0.635  -0.635   0.635  1.00  0.00           H  
HETATM    5  H4  UNL     1       0.635  -0.635  -0.635  1.00  0.00           H  
CONECT    1    2    3    4    5
END"""

# Example script data - content must be an array for scriptData.content.map() to work
script_data: Dict[str, Any] = {
    "title": "Methane: The Simplest Hydrocarbon",
    "content": [  # This must be an array, not an object
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

def main():
    # Get the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Test 1: Generate interactive HTML with both script data and PDB data
    print("Test 1: Generating interactive HTML with both script data and PDB data...")
    output_path = os.path.join(current_dir, 'test_interactive_output.html')
    result = MoleculeVisualizer.generate_interactive_html(
        pdb_data=pdb_data,
        title="Methane",
        script_data=script_data,
        output_path=output_path
    )
    print(f"Output saved to: {result}")
    
    # Test 2: Generate interactive HTML with only PDB data and title
    print("\nTest 2: Generating interactive HTML with only PDB data and title...")
    output_path2 = os.path.join(current_dir, 'test_interactive_simple.html')
    result2 = MoleculeVisualizer.generate_interactive_html(
        pdb_data=pdb_data,
        title="Methane",
        output_path=output_path2
    )
    print(f"Output saved to: {result2}")
    
    # Test 3: Generate interactive HTML as a string
    print("\nTest 3: Generating interactive HTML as a string...")
    html_content = MoleculeVisualizer.generate_interactive_html(
        pdb_data=pdb_data,
        title="Methane"
    )
    print(f"Generated HTML content length: {len(html_content)} characters")
    
    # Save the string content to a file for verification
    string_output_path = os.path.join(current_dir, 'test_interactive_string.html')
    with open(string_output_path, 'w') as f:
        f.write(html_content)
    print(f"String content saved to: {string_output_path}")

if __name__ == "__main__":
    main() 