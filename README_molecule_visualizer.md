# Molecule Visualizer Generator

This tool allows you to generate molecule visualization HTML files by injecting custom script and PDB data into an HTML template. It avoids the need to store the entire HTML content as a Python string with escaped characters.

## Features

- Inject custom script data (animation steps, captions, etc.) into an HTML template
- Inject custom PDB data (molecule structure) into an HTML template
- Use either file paths or direct data in your Python code
- Simple command-line interface for quick generation

## Installation

No installation required. Just copy the `molecule_visualizer.py` script to your project directory.

## Usage

### Command Line

```bash
# Generate a visualization using file paths
python molecule_visualizer.py --template visualization-40.html --output ethane.html --script-data example_script_data.json --pdb-data example_pdb_data.txt
```

### Python API

```python
from molecule_visualizer import generate_visualization

# Example 1: Using file paths
generate_visualization(
    template_path='visualization-40.html',
    output_path='ethane.html',
    script_data_path='example_script_data.json',
    pdb_data_path='example_pdb_data.txt'
)

# Example 2: Using direct data
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

pdb_data = """COMPND    297
HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  
# ... more PDB data ...
END"""

generate_visualization(
    template_path='visualization-40.html',
    output_path='methane.html',
    script_data=script_data,
    pdb_data=pdb_data
)
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

- `example_script_data.json`: Example script data for ethane
- `example_pdb_data.txt`: Example PDB data for ethane
- `example_usage.py`: Example Python script demonstrating how to use the generator

## How It Works

The script uses regular expressions to find and replace the `scriptData` and `pdbData` constants in the HTML template. This approach avoids having to store the entire HTML content as a Python string with escaped characters, making it easier to maintain and update the template.

## License

MIT 