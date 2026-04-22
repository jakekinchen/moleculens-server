# Moleculens Server

Production-ready FastAPI microservice for Psi4 cubeprop orbital/density overlays for Three.js.

## Overview

Given a 3D SDF (with coordinates), this service computes molecular orbital overlays (HOMO/LUMO + optional density), converts them to meshes server-side, caches results, and returns data in a Three.js-friendly format.

### Key Features

- **Pre-meshed surfaces**: Returns ready-to-render mesh data (vertices, normals, indices)
- **Dual lobe support**: Separate meshes for positive and negative orbital lobes
- **Compressed output**: Gzipped base64-encoded typed arrays for efficient transfer
- **Caching**: Results cached by geometry + parameters for fast repeated queries
- **Postgres job queue**: Scalable async processing with de-duplication

## Quick Start

### Development

```bash
# Install dependencies
uv sync --dev

# Start services (includes Postgres)
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d

# API available at http://localhost:8001
```

### Production

```bash
# Set required environment variables
export DATABASE_URL=postgresql://user:pass@host:5432/moleculens
export ALLOWED_ORIGINS=https://your-app.com

# Build and run
docker compose up -d
```

## API Endpoints

### Health Check
```http
GET /health
```

### Submit Computation
```http
POST /v1/orbitals/compute
Content-Type: application/json

{
  "sdfContent": "...",
  "method": "scf",
  "basis": "sto-3g",
  "orbitals": ["homo", "lumo", "density"],
  "gridSpacing": 0.25,
  "isovalue": 0.05
}
```

### Poll Job Status
```http
GET /v1/orbitals/jobs/{jobId}
```

### Download Artifact
```http
GET /v1/orbitals/artifacts/{cacheKey}/{filename}
```

### Submit Conformer
```http
POST /v1/conformers/compute
Content-Type: application/json

{
  "molblock2d": "...",
  "params": {
    "method": "etkdg_v3",
    "opt": "uff",
    "maxAttempts": 10,
    "maxOptIters": 200,
    "addHs": false
  },
  "waitMs": 800,
  "geomVersion": "2026-01-07_v1"
}
```

### Poll Conformer Job Status
```http
GET /v1/conformers/jobs/{jobId}
```

### Download Conformer Artifact
```http
GET /v1/conformers/artifacts/{cacheKey}/structure.sdf
```

## Response Format

```json
{
  "cached": false,
  "jobId": "...",
  "status": "done",
  "cacheKey": "...",
  "result": {
    "orbitals": {
      "homo": {
        "positive": {
          "vertices": "<base64 gzipped Float32Array>",
          "normals": "<base64 gzipped Float32Array>",
          "indices": "<base64 gzipped Uint32Array>",
          "vertexCount": 1234,
          "triangleCount": 2345
        },
        "negative": { ... },
        "energyEv": -12.34,
        "isovalue": 0.05
      },
      "lumo": { ... }
    },
    "density": { ... },
    "meta": {
      "method": "scf",
      "basis": "sto-3g",
      "gridSpacingAngstrom": 0.25,
      "psi4NumThreads": 3,
      "computeTimeMs": 5432.1
    }
  }
}
```

## Three.js Integration

```javascript
// Decode mesh data
function decodeMesh(meshData) {
  const compressed = atob(meshData.vertices);
  const decompressed = pako.inflate(compressed);
  const vertices = new Float32Array(decompressed.buffer);

  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
  // ... similarly for normals and indices

  return geometry;
}
```

## Configuration

| Environment Variable | Default | Description |
|---------------------|---------|-------------|
| DATABASE_URL | required | PostgreSQL connection string |
| MOLECULE_CACHE_DIR | /var/molecule-cache | Artifact cache directory |
| PSI_SCRATCH | /var/psi4_scratch | Psi4 scratch directory |
| PSI4_NUM_THREADS | 3 | Psi4 thread count |
| PSI4_MEMORY_MB | 4096 | Psi4 memory limit |
| WORKER_POLL_SECONDS | 1 | Worker poll interval |
| ALLOWED_ORIGINS | * | CORS allowed origins |
| LOG_LEVEL | INFO | Logging level |
| METRICS_AUTH_USER | metrics | Basic auth username for `GET /metrics/cache` |
| METRICS_AUTH_PASSWORD_HASH | unset | BCrypt hash for the metrics route password |

The cache metrics endpoint is intended for operator use, not anonymous public access. In the single-box Caddy deployment it is exposed at `https://api.moleculens.com/metrics/cache` only after HTTP basic auth succeeds.

Generate a password hash with Caddy and add it to your `.env` before deploying:

```bash
docker exec molecule-caddy caddy hash-password --plaintext 'choose-a-strong-password'
```

Then query it with:

```bash
curl -u "$METRICS_AUTH_USER:your-plaintext-password" \
  https://api.moleculens.com/metrics/cache
```

## Architecture

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Client    │────▶│  API (FastAPI)  │────▶│  PostgreSQL │
└─────────────┘     └─────────────┘     └──────┬──────┘
                                               │
                    ┌─────────────┐            │
                    │   Worker    │◀───────────┘
                    │   (Psi4)    │
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │   Cache     │
                    │   (Disk)    │
                    └─────────────┘
```

## Development

```bash
# Install dev dependencies
uv sync --dev

# Run linter
uv run ruff check .

# Run formatter
uv run ruff format .

# Run unit tests
uv run pytest tests/unit -v

# Run integration tests (requires Docker)
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d
uv run pytest tests/integration -v
```

## License

MIT
