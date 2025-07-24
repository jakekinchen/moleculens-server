# Moleculens Server - Project Structure

## Overview
This document describes the organized structure of the Moleculens server codebase.

## Directory Structure

```
moleculens-server/
├── api/                          # Main application code
│   ├── agent_management/         # AI agents and LLM services
│   ├── dependencies/             # FastAPI dependencies
│   ├── diagram/                  # Diagram generation pipeline
│   ├── routers/                  # API route handlers
│   ├── static/                   # Static files and test outputs
│   ├── tests/                    # Legacy API tests (being migrated)
│   └── utils/                    # Shared utilities
├── tests/                        # Organized test suite
│   ├── integration/              # End-to-end pipeline tests
│   ├── unit/                     # Isolated component tests
│   ├── demos/                    # Demo and example tests
│   └── fixtures/                 # Test data (in api/tests/fixtures)
├── examples/                     # Example usage and demos
│   ├── calvin_cycle/             # Calvin cycle photosynthesis examples
│   ├── pipeline_demos/           # Complete pipeline demonstrations
│   └── compose_rendered_images.py
├── sample_outputs/               # Example generated files
├── docs/                         # Documentation files
│   ├── AGENTS.md                 # Agent development roadmap
│   ├── COMPLETE_PIPELINE_DEMO.md # Pipeline documentation
│   ├── MOLECULE_DIAGRAM_PLAN.md  # Diagram generation plan
│   └── TRANSPARENT_RENDERING_IMPLEMENTATION.md
├── .env                          # Environment variables
├── .gitignore                    # Git ignore rules
├── .pre-commit-config.yaml       # Pre-commit hooks
├── Dockerfile                    # Container configuration
├── README.md                     # Main project documentation
├── pytest.ini                   # Test configuration
└── setup.cfg                    # Project configuration
```

## Key Components

### API (`api/`)
- **agent_management/**: AI agents for molecular data processing
- **routers/**: FastAPI endpoints organized by feature
- **diagram/**: Molecular diagram generation pipeline
- **utils/**: Shared utilities (caching, security, rate limiting)

### Tests (`tests/`)
- **integration/**: Full pipeline and API endpoint tests
- **unit/**: Individual component tests
- **demos/**: Example usage and demonstration scripts

### Examples (`examples/`)
- **calvin_cycle/**: Photosynthesis pathway examples
- **pipeline_demos/**: Complete workflow demonstrations

## Running the Project

```bash
# Install dependencies
pip install -r api/requirements.txt

# Run the server
python api/main.py

# Run tests
pytest                    # All tests
pytest tests/unit/        # Unit tests only
pytest tests/integration/ # Integration tests only
```

## Development Workflow

1. **API Development**: Add new features in `api/` with appropriate tests
2. **Testing**: Write unit tests in `tests/unit/`, integration tests in `tests/integration/`
3. **Examples**: Add usage examples in `examples/` for new features
4. **Documentation**: Update relevant docs in the root directory

## Migration Notes

- Legacy test files have been moved from root to `tests/` subdirectories
- API tests in `api/tests/` are being gradually migrated to the main `tests/` directory
- Demo files and sample outputs are now properly organized
