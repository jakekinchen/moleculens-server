#!/usr/bin/env python3
"""
Example usage of the molecule_visualizer.py script.

This script demonstrates how to use the molecule_visualizer.py script
to generate a visualization HTML file with custom script and PDB data.
"""

import os
from molecule_visualizer import generate_visualization

# Define paths
current_dir = os.path.dirname(os.path.abspath(__file__))
template_path = os.path.join(current_dir, 'visualization-40.html')
output_path = os.path.join(current_dir, 'ethane_visualization.html')

# Example 1: Using file paths
print("Example 1: Using file paths")
generate_visualization(
    template_path=template_path,
    output_path=output_path,
    script_data_path=os.path.join(current_dir, 'example_script_data.json'),
    pdb_data_path=os.path.join(current_dir, 'example_pdb_data.txt')
)

# Example 2: Using direct data
print("\nExample 2: Using direct data")
# Define script data directly
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

# Define PDB data directly
pdb_data = """COMPND    297
HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  
HETATM    2  H1  UNL     1       0.635   0.635   0.635  1.00  0.00           H  
HETATM    3  H2  UNL     1      -0.635   0.635  -0.635  1.00  0.00           H  
HETATM    4  H3  UNL     1      -0.635  -0.635   0.635  1.00  0.00           H  
HETATM    5  H4  UNL     1       0.635  -0.635  -0.635  1.00  0.00           H  
CONECT    1    2    3    4    5
END"""

# Generate visualization with direct data
generate_visualization(
    template_path=template_path,
    output_path=os.path.join(current_dir, 'methane_visualization.html'),
    script_data=script_data,
    pdb_data=pdb_data
)

print("\nDone! Check the generated HTML files.") 