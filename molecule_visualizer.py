#!/usr/bin/env python3
"""
Molecule Visualizer Generator

This script takes an HTML template file and injects new values for scriptData and pdbData
without needing to store the entire HTML content as a Python string.
"""

import re
import json
import os
import argparse
from pathlib import Path


def generate_visualization(
    template_path,
    output_path,
    script_data=None,
    pdb_data=None,
    script_data_path=None,
    pdb_data_path=None
):
    """
    Generate a visualization HTML file by injecting new script and PDB data.
    
    Args:
        template_path (str): Path to the HTML template file
        output_path (str): Path where the generated HTML file will be saved
        script_data (dict, optional): New script data to inject
        pdb_data (str, optional): New PDB data to inject
        script_data_path (str, optional): Path to a JSON file containing script data
        pdb_data_path (str, optional): Path to a file containing PDB data
    """
    # Read the template file
    with open(template_path, 'r') as f:
        template_content = f.read()
    
    # Load script data from file if provided
    if script_data_path and not script_data:
        with open(script_data_path, 'r') as f:
            script_data = json.load(f)
    
    # Load PDB data from file if provided
    if pdb_data_path and not pdb_data:
        with open(pdb_data_path, 'r') as f:
            pdb_data = f.read()
    
    # Replace script data if provided
    if script_data:
        # Convert script data to JSON string with proper formatting
        script_json = json.dumps(script_data, indent=2)
        
        # Find the scriptData constant in the template
        script_pattern = r'const scriptData = \{[\s\S]*?\};'
        script_match = re.search(script_pattern, template_content)
        
        if script_match:
            # Replace the scriptData constant
            template_content = template_content[:script_match.start()] + \
                               f'const scriptData = {script_json};' + \
                               template_content[script_match.end():]
    
    # Replace PDB data if provided
    if pdb_data:
        # Find the pdbData constant in the template
        pdb_pattern = r'const pdbData = `[\s\S]*?`;'
        pdb_match = re.search(pdb_pattern, template_content)
        
        if pdb_match:
            # Replace the pdbData constant
            template_content = template_content[:pdb_match.start()] + \
                               f'const pdbData = `{pdb_data}`;' + \
                               template_content[pdb_match.end():]
    
    # Write the modified content to the output file
    with open(output_path, 'w') as f:
        f.write(template_content)
    
    print(f"Visualization generated successfully: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Generate molecule visualization HTML files')
    parser.add_argument('--template', required=True, help='Path to the HTML template file')
    parser.add_argument('--output', required=True, help='Path where the generated HTML file will be saved')
    parser.add_argument('--script-data', help='Path to a JSON file containing script data')
    parser.add_argument('--pdb-data', help='Path to a file containing PDB data')
    
    args = parser.parse_args()
    
    generate_visualization(
        template_path=args.template,
        output_path=args.output,
        script_data_path=args.script_data,
        pdb_data_path=args.pdb_data
    )


if __name__ == '__main__':
    main() 