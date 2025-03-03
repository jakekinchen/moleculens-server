#!/usr/bin/env python3
"""
Test script for file-based injection using the inject_data_into_template method.

This script demonstrates how to use the inject_data_into_template method
with file paths instead of direct data.
"""

import os
from molecule_visualizer import MoleculeVisualizer

def main():
    # Get the current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define file paths
    script_data_path = os.path.join(current_dir, 'example_script_data.json')
    pdb_data_path = os.path.join(current_dir, 'example_pdb_data.txt')
    
    # Check if the files exist
    if not os.path.exists(script_data_path):
        print(f"Error: Script data file not found: {script_data_path}")
        print("Please run test_inject_data.py first to create the example files.")
        return
    
    if not os.path.exists(pdb_data_path):
        print(f"Error: PDB data file not found: {pdb_data_path}")
        print("Please run test_inject_data.py first to create the example files.")
        return
    
    # Test file-based injection
    print("Testing file-based injection...")
    output_path = MoleculeVisualizer.inject_data_into_template(
        script_data_path=script_data_path,
        pdb_data_path=pdb_data_path
    )
    print(f"Output saved to: {output_path}")

if __name__ == "__main__":
    main() 