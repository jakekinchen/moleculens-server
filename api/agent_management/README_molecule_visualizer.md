# Molecule Visualizer in Agent Management

This module provides functionality for visualizing molecular structures in the agent management system. It includes methods for converting SDF data to PDB format, generating JavaScript code for embedding in Three.js scenes, and injecting custom script and PDB data into the output.html template.

## Features

- Convert SDF data to PDB format in-memory
- Generate minimal JavaScript code for embedding in existing Three.js scenes
- Generate standalone HTML viewers for molecules
- Show toggleable atomic annotations that face the camera
- Inject custom script data (animation steps, captions, etc.) into the output.html template
- Inject custom PDB data (molecule structure) into the output.html template

## Usage

### Injecting Data into the Template

The `generate_interactive_html` method allows you to generate an interactive HTML visualization by injecting PDB data and optional script data into the output.html template:

```python
from molecule_visualizer import MoleculeVisualizer

# Example 1: Generate HTML content as a string
pdb_data = """COMPND    297
HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C
# ... more PDB data ...
END"""

script_data = {
    "title": "Methane: The Simplest Hydrocarbon",
    "content": [
        {
            "timecode": "00:00",
            "atoms": [],
            "caption": "Methane is the simplest hydrocarbon with the molecular formula CH₄."
        },
        # ... more steps ...
    ]
}

# Generate HTML content as a string
html_content = MoleculeVisualizer.generate_interactive_html(
    pdb_data=pdb_data,
    title="Methane",
    script_data=script_data
)

# Example 2: Save HTML content to a file
output_path = MoleculeVisualizer.generate_interactive_html(
    pdb_data=pdb_data,
    title="Methane",
    script_data=script_data,
    output_path="methane_visualization.html"
)

# Example 3: Generate HTML with only PDB data and title (auto-generates simple script)
simple_html = MoleculeVisualizer.generate_interactive_html(
    pdb_data=pdb_data,
    title="Methane"
)
```

### Other Visualization Methods

The module also provides methods for generating JavaScript code and HTML viewers from SDF data:

```python
from molecule_visualizer import MoleculeVisualizer

# Example SDF data
sdf_data = """
  Molecule Name

  10 10  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    # ... more SDF data ...
$$$$
"""

# Generate JavaScript code for embedding in an existing Three.js scene
js_code = MoleculeVisualizer.generate_js_code_from_sdf(sdf_data, name="Molecule Name")

# Generate a standalone HTML viewer
html_content = MoleculeVisualizer.generate_html_viewer_from_sdf(sdf_data, name="Molecule Name")
```

## File Formats

### Script Data (JSON)

The script data should be a JSON object with the following structure:

```json
{
  "title": "Molecule Name: Description",
  "content": [
    {
      "timecode": "00:00",
      "atoms": [],
      "caption": "Initial description of the molecule."
    },
    {
      "timecode": "00:05",
      "atoms": ["0", "1"],
      "caption": "Description of highlighted atoms."
    },
    // ... more steps ...
  ]
}
```

**Important Note**: The `content` property must be an array (not an object) for the visualization to work correctly. The JavaScript code calls `scriptData.content.map()` which requires an array.

- `timecode`: Time in "MM:SS" format when this step should be shown
- `atoms`: Array of atom indices to highlight (as strings)
- `caption`: Text description to display

### PDB Data

The PDB data should be in standard Protein Data Bank format:

```
COMPND    MOLECULE_ID
HETATM    1  ATOM_NAME  UNL     1       X_COORD  Y_COORD  Z_COORD  1.00  0.00           ELEMENT
HETATM    2  ATOM_NAME  UNL     1       X_COORD  Y_COORD  Z_COORD  1.00  0.00           ELEMENT
...
CONECT    1    2    3    4
...
END
```

## Example Files

- `example_script_data.json`: Example script data for methane
- `example_pdb_data.txt`: Example PDB data for methane
- `test_inject_data.py`: Example script demonstrating how to use the inject_data_into_template method

## How It Works

The `inject_data_into_template` method uses regular expressions to find and replace the `scriptData` and `pdbData` constants in the output.html template. This approach avoids having to store the entire HTML content as a Python string with escaped characters, making it easier to maintain and update the template.
