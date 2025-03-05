# Atom Label Conversion System

This document describes the bidirectional atom label conversion system implemented in the PubChem Agent to improve molecule script generation.

## Overview

The system supports two types of atom labels:

1. **Numeric Indices** (0, 1, 2...)
   - Used internally for visualization
   - More compact representation
   - Less human-readable

2. **Element-Based Labels** (C1, H2, O1...)
   - More human-readable
   - Better for LLM script generation
   - Easier for non-technical users to understand

## Key Components

### 1. Label Generation

The system uses three methods to generate element-based labels, with prioritization:

1. **Primary**: Extract labels from PDB/SDF data
   ```
   HETATM    1  C1  UNL     1       3.814  -1.599   0.000  1.00  0.00           C
   ```

2. **Secondary**: Generate from elements array
   ```python
   elements = ['C', 'C', 'O', 'H', 'H']
   # Generates: C1, C2, O1, H1, H2
   ```

3. **Fallback**: Parse from SMILES using RDKit
   ```python
   molecule = Chem.MolFromSmiles('CCO')
   # Generates: C1, C2, O1
   ```

### 2. Conversion Functions

The system provides the following functions:

- `generate_atom_label_mapping(molecule_data)` - Creates mapping from indices to element labels
- `reverse_atom_label_mapping(mapping)` - Creates inverse mapping from element labels to indices
- `validate_and_convert_script(script, molecule_data, use_element_labels, convert_back_to_indices)` - Handles script validation and bidirectional conversion

### 3. PubChemAgent Configuration

The PubChemAgent can be configured with two parameters:

```python
agent = AgentFactory.create_pubchem_agent(
    use_element_labels=True,        # Use element-based labels for LLM
    convert_back_to_indices=False   # Keep as element labels in the output
)
```

## Workflow

1. **LLM Script Generation**: 
   - Convert numeric indices to element labels (C1, H2)
   - Send to LLM for script generation
   - LLM understands element references better

2. **Output Options**:
   - Keep as element labels (for human reading)
   - Convert back to numeric indices (for visualization)

## Testing

The system includes comprehensive tests:
- Unit tests using mocked implementations
- Integration tests with RDKit
- Tests for all generation methods and bidirectional conversion

## Usage Example

```python
# Create PubChemAgent with element-based labels
pubchem_agent = AgentFactory.create_pubchem_agent(
    use_element_labels=True,
    convert_back_to_indices=False
)

# Generate a molecule package with element-based labels in the script
molecule_package = pubchem_agent.get_molecule_package("Show me water")

# The script will contain element-based labels like O1, H1, H2
# instead of numeric indices like 0, 1, 2
```