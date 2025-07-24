# Examples and Demos

This directory contains example usage, demo scripts, and sample outputs for the Moleculens server.

## Directory Structure

```
examples/
├── calvin_cycle/      # Calvin cycle photosynthesis examples
├── pipeline_demos/    # Complete pipeline demonstrations
├── compose_rendered_images.py  # Image composition utilities
└── README.md
```

## Sample Outputs

The `sample_outputs/` directory in the project root contains example generated files:
- SVG diagrams
- YAML specifications
- Transparent molecular images

## Running Examples

Most example files can be run directly:

```bash
python examples/compose_rendered_images.py
python tests/demos/test_calvin_cycle.py
```

Note: Some examples may require environment setup (API keys, dependencies).
