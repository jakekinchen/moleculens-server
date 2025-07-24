# Transparent Background Rendering with PyMOL

## Overview

The Moleculens server now supports transparent background rendering for creating molecular visualizations suitable for presentations, infographics, and overlays. This feature uses PyMOL's ray-tracing capabilities with alpha channel support.

## Key Features

### 1. Transparent Background Support
- **Command**: `set ray_opaque_background, 0`
- **Result**: PNG files with alpha channel transparency
- **Use case**: Overlay molecules on custom backgrounds

### 2. Enhanced Rendering Quality
- **Ray tracing**: High-quality software rendering
- **Antialiasing**: Smooth edges and lines
- **Shadows**: Realistic depth perception
- **Custom DPI**: Publication-quality resolution

### 3. LLM Integration
The system can understand natural language requests for transparency:
- "Show me caffeine with a transparent background"
- "Create a molecule with no background for presentation"
- "Generate transparent PNG of protein structure"

## API Usage

### Basic Transparent Rendering
```json
{
  "description": "Show me 1ubq protein with transparent background",
  "format": "image",
  "transparent_background": true,
  "ray_trace": true,
  "resolution": [1920, 1080],
  "dpi": 300
}
```

### Advanced Rendering Options
```json
{
  "description": "Publication quality caffeine molecule",
  "format": "image",
  "transparent_background": true,
  "ray_trace": true,
  "ray_trace_mode": "default",
  "resolution": [1920, 1080],
  "dpi": 300,
  "antialias": true,
  "ray_shadow": true,
  "depth_cue": true,
  "background_color": "white"
}
```

## Available Templates

### 1. Transparent Molecule Scene
```python
transparent_molecule_scene(structure_id="1ubq", style="cartoon")
```
- Basic molecule with transparent background
- Chainbow coloring
- Optimized for overlays

### 2. Publication Quality Scene
```python
publication_quality_scene(structure_id="1ubq", transparent=True)
```
- High-quality rendering settings
- Optional highlighting
- Professional appearance

### 3. Transparent Binding Site
```python
transparent_binding_site_scene(structure_id="1abc", selection="resi 45", transparent_bg=True)
```
- Focused binding site visualization
- Transparent background
- Suitable for presentations

## Essential PyMOL Commands

The system automatically generates these key commands for transparency:

```pymol
# Enable transparency
set ray_opaque_background, 0

# Quality settings
set antialias, 1
set ray_shadow, 1
set depth_cue, 1

# High-resolution rendering
ray 1920, 1080
png molecule.png, dpi=300, ray=1
```

## Use Cases

### 1. Presentation Slides
- Overlay molecules on custom backgrounds
- Maintain visual consistency with slide themes
- Professional appearance

### 2. Infographics
- Combine molecular structures with diagrams
- Create educational materials
- Scientific publications

### 3. Web Graphics
- Transparent PNGs for web overlays
- Interactive visualizations
- Dynamic content

## Technical Details

### Ray Tracing Settings
- **Mode 0**: Default ray tracing
- **Mode 1**: Poster-style rendering
- **Mode 2**: Black and white
- **Mode 3**: Cartoon outline

### Alpha Channel Support
- PNG format preserves transparency
- Compatible with image editing software
- Suitable for compositing

### Performance Considerations
- Ray tracing is CPU-intensive
- Higher resolutions increase render time
- Caching reduces repeated requests

## Example Workflow

1. **Request**: "Show me a transparent caffeine molecule for my presentation"
2. **LLM Processing**: Identifies transparent background requirement
3. **Template Selection**: Uses `transparent_molecule_scene`
4. **PyMOL Commands**: Generates commands with `ray_opaque_background, 0`
5. **Rendering**: Ray-traces at specified resolution with alpha channel
6. **Output**: PNG file with transparent background

## Security

All transparent rendering commands are validated against the security whitelist:
- ✅ `set ray_opaque_background, 0`
- ✅ `set antialias, 1`
- ✅ `ray 1920, 1080`
- ✅ `png output.png, dpi=300, ray=1`

## Integration with Molecule Diagrams

The transparent background feature is designed to work seamlessly with the molecule diagram pipeline, enabling:
- Transparent 3D structures in 2D infographics
- Layered visualizations
- Mixed media presentations

This enhancement makes the Moleculens server ideal for creating professional molecular visualizations suitable for presentations, publications, and educational materials.
