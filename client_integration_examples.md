# Client Integration Guide: 2D Transparent PNGs & 3D Data

This guide shows how to get both 2D transparent PNGs and 3D data from the MoleculeLens API for client-side rendering.

## 🎯 Quick Summary

| Data Type | Endpoint | Format | Best For |
|-----------|----------|--------|----------|
| **2D Transparent PNG** | `/render` | `format: "image"` | UI overlays, diagrams, thumbnails |
| **3D PDB Data** | `/render` | `format: "model"` | Three.js reconstruction |
| **3D PDB (Known)** | `/rcsb/fetch-structure/` | `format: "pdb"` | Known protein structures |

## 🖼️ Getting 2D Transparent PNGs

### Basic Usage
```javascript
const get2DTransparentPNG = async (moleculeName) => {
  const response = await fetch('/render', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      description: `Show ${moleculeName} molecule with transparent background`,
      format: "image",
      transparent_background: true,
      ray_trace: true,
      resolution: [512, 512],
      dpi: 150
    })
  });

  // Handle both direct PNG and large file URL responses
  if (response.headers.get('content-type')?.includes('image/png')) {
    return response.blob(); // Direct PNG blob
  } else {
    const data = await response.json();
    const imageResponse = await fetch(data.url);
    return imageResponse.blob(); // PNG from URL
  }
};

// Usage
const pngBlob = await get2DTransparentPNG('caffeine');
const imageUrl = URL.createObjectURL(pngBlob);
document.getElementById('molecule-img').src = imageUrl;
```

### High-Quality 2D Rendering
```javascript
const getHighQuality2D = async (moleculeName) => {
  return fetch('/render', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      description: `Show ${moleculeName} molecule structure with transparent background`,
      format: "image",
      transparent_background: true,
      ray_trace: true,
      resolution: [1024, 1024],
      dpi: 300,
      ray_trace_mode: "poster", // Publication quality
      antialias: true,
      ray_shadow: false, // Better for transparent backgrounds
      background_color: "white"
    })
  });
};
```

## 🧬 Getting 3D PDB Data for Three.js

### Method 1: Via Render Endpoint (Any Molecule)
```javascript
const get3DPDBData = async (moleculeName) => {
  const response = await fetch('/render', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      description: `Load ${moleculeName} molecule structure`,
      format: "model" // Returns PDB format
    })
  });

  if (response.headers.get('content-type')?.includes('chemical/x-pdb')) {
    return response.text(); // Direct PDB text
  } else {
    const data = await response.json();
    const pdbResponse = await fetch(data.url);
    return pdbResponse.text(); // PDB from URL
  }
};

// Usage with Three.js PDB loader
import { PDBLoader } from 'three/examples/jsm/loaders/PDBLoader.js';

const load3DMolecule = async (moleculeName) => {
  const pdbData = await get3DPDBData(moleculeName);

  const loader = new PDBLoader();
  const pdb = loader.parse(pdbData);

  // Add to scene
  scene.add(pdb);

  return pdb;
};
```

### Method 2: Via RCSB Endpoint (Known Structures)
```javascript
const getPDBFromRCSB = async (pdbId) => {
  const response = await fetch('/rcsb/fetch-structure/', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      identifier: pdbId, // e.g., "1ubq", "2itx"
      format: "pdb"
    })
  });

  const data = await response.json();
  return data.data; // PDB text content
};

// Usage
const ubiquitinPDB = await getPDBFromRCSB('1ubq');
```

## 🚀 Enhanced Endpoints (Recommended)

I've created enhanced endpoints that make this even easier:

### Single Molecule (Both 2D + 3D)
```javascript
const getMoleculeData = async (moleculeName) => {
  const response = await fetch('/render/molecule', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      molecule_name: moleculeName,
      render_type: "both", // "2d_transparent", "3d_pdb", or "both"
      size: "medium", // "small", "medium", "large"
      quality: "high" // "fast", "high", "publication"
    })
  });

  const data = await response.json();
  return {
    name: data.molecule_name,
    png: data.png_base64 ? `data:image/png;base64,${data.png_base64}` : data.png_url,
    pdb: data.pdb_data,
    metadata: data.metadata
  };
};

// Usage
const { png, pdb } = await getMoleculeData('caffeine');
document.getElementById('molecule-img').src = png;
// Load PDB into Three.js...
```

### Batch Processing
```javascript
const getBatchMolecules = async (molecules) => {
  const response = await fetch('/render/batch', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      molecules: molecules,
      render_type: "2d_transparent",
      size: "small",
      quality: "fast"
    })
  });

  return response.json();
};

// Usage
const results = await getBatchMolecules(['caffeine', 'aspirin', 'glucose']);
results.forEach(result => {
  if (result.png_base64) {
    // Display 2D image
    const img = document.createElement('img');
    img.src = `data:image/png;base64,${result.png_base64}`;
    document.body.appendChild(img);
  }
});
```

### Protein Structures
```javascript
const getProteinData = async (pdbId) => {
  const response = await fetch('/render/protein', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      pdb_id: pdbId,
      render_type: "both",
      representation: "cartoon", // "cartoon", "surface", "stick", "ribbon"
      size: "medium",
      quality: "high"
    })
  });

  return response.json();
};

// Usage
const protein = await getProteinData('1ubq');
```

## 🎨 Three.js Integration Example

```javascript
import * as THREE from 'three';
import { PDBLoader } from 'three/examples/jsm/loaders/PDBLoader.js';

class MoleculeViewer {
  constructor(containerId) {
    this.container = document.getElementById(containerId);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    this.renderer = new THREE.WebGLRenderer({ alpha: true });

    this.renderer.setSize(window.innerWidth, window.innerHeight);
    this.container.appendChild(this.renderer.domElement);

    this.loader = new PDBLoader();
  }

  async loadMolecule(moleculeName) {
    try {
      // Get PDB data from API
      const pdbData = await get3DPDBData(moleculeName);

      // Parse with Three.js
      const pdb = this.loader.parse(pdbData);

      // Style the molecule
      pdb.traverse((child) => {
        if (child.isMesh) {
          child.material.transparent = true;
          child.material.opacity = 0.8;
        }
      });

      // Add to scene
      this.scene.add(pdb);

      // Position camera
      const box = new THREE.Box3().setFromObject(pdb);
      const center = box.getCenter(new THREE.Vector3());
      const size = box.getSize(new THREE.Vector3());

      this.camera.position.copy(center);
      this.camera.position.z += size.length();
      this.camera.lookAt(center);

      return pdb;
    } catch (error) {
      console.error('Failed to load molecule:', error);
    }
  }

  render() {
    this.renderer.render(this.scene, this.camera);
  }
}

// Usage
const viewer = new MoleculeViewer('molecule-container');
await viewer.loadMolecule('caffeine');
viewer.render();
```

## 📊 Performance Considerations

### For 2D Images:
- Use `size: "small"` and `quality: "fast"` for thumbnails
- Use `size: "medium"` and `quality: "high"` for main displays
- Use `size: "large"` and `quality: "publication"` for detailed views

### For 3D Data:
- PDB files are typically small (< 1MB) and load quickly
- Cache PDB data on client side to avoid repeated requests
- Use the batch endpoint for multiple molecules

### Caching Strategy:
```javascript
class MoleculeCache {
  constructor() {
    this.cache = new Map();
  }

  async getMolecule(name, type = 'both') {
    const key = `${name}_${type}`;

    if (this.cache.has(key)) {
      return this.cache.get(key);
    }

    const data = await getMoleculeData(name);
    this.cache.set(key, data);

    return data;
  }
}

const cache = new MoleculeCache();
```

## 🔧 Error Handling

```javascript
const safeGetMolecule = async (moleculeName) => {
  try {
    const data = await getMoleculeData(moleculeName);
    return data;
  } catch (error) {
    console.error(`Failed to load ${moleculeName}:`, error);

    // Fallback to basic structure fetch
    try {
      const pdbData = await getPDBFromRCSB(moleculeName);
      return { pdb: pdbData, png: null };
    } catch (fallbackError) {
      console.error('Fallback also failed:', fallbackError);
      return null;
    }
  }
};
```

This setup gives you maximum flexibility - transparent 2D PNGs for UI elements and full 3D PDB data for interactive Three.js visualization!
