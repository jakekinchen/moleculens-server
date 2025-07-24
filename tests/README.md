# Test Organization

This directory contains all tests for the Moleculens server, organized by type:

## Directory Structure

```
tests/
├── integration/     # End-to-end pipeline tests
├── unit/           # Isolated component tests
├── demos/          # Demo and example tests
└── fixtures/       # Test data and fixtures
```

## Running Tests

```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/unit/

# Run only integration tests
pytest tests/integration/

# Run tests with specific markers
pytest -m unit
pytest -m integration
pytest -m "not slow"
```

## Test Categories

- **Unit Tests**: Test individual components in isolation
- **Integration Tests**: Test complete pipelines and workflows
- **Demo Tests**: Example usage and demonstration scripts
