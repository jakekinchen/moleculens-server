# Moleculens Server

A FastAPI-based server for molecular visualization and analysis, featuring AI-powered molecule rendering and diagram generation.

## Project Structure

See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for detailed organization of the codebase.

```
moleculens-server/
├── api/              # Main application code
├── tests/            # Organized test suite
├── examples/         # Usage examples and demos
├── sample_outputs/   # Example generated files
├── docs/             # Documentation
└── README.md         # This file
```

## Features

- 3D molecular visualization using PyMOL
- 2D molecular diagram generation
- AI-powered prompt interpretation
- PubChem integration
- RCSB integration
- Real-time rendering and animation

## Prerequisites

- Docker and Docker Compose
- Python 3.9+
- OpenAI API key (for AI features)
- Git

## Quick Start

1. Clone the repository:
```bash
git clone https://github.com/jakekinchen/moleculens-server.git
cd moleculens-server
```

2. Create a `.env` file in the root directory:
```env
OPENAI_API_KEY=your-api-key
ENVIRONMENT=development  # or production
```

3. Build and start the containers:
```bash
docker-compose up --build
```

The server will be available at `http://localhost:8000`.

### Docker Commands

- **Stop containers:**
  ```bash
  docker-compose down
  ```

- **Rebuild without cache:**
  ```bash
  docker-compose build --no-cache
  ```

- **View logs:**
  ```bash
  docker-compose logs -f
  ```

## Project Structure

```
moleculens-server/
├── api/                    # Main application code
│   ├── agent_management/   # AI agent configuration
│   ├── dependencies/       # FastAPI dependencies
│   ├── routers/           # API routes
│   ├── static/            # Static files
│   ├── tests/             # Test suite
│   │   └── fixtures/      # Test data
│   ├── utils/             # Utility functions
│   └── main.py           # Application entry point
├── docs/                  # Documentation
└── docker-compose.yml    # Docker configuration
```

## Development

### Running Tests

The test suite is organized into categories:

```bash
# Run all tests
pytest

# Run by directory
pytest tests/unit/           # Unit tests only
pytest tests/integration/    # Integration tests only
pytest tests/demos/          # Demo/example tests

# Run by markers
pytest -m unit              # Unit tests
pytest -m integration       # Integration tests
pytest -m "not slow"        # Exclude slow tests
```

See [tests/README.md](tests/README.md) for detailed test organization.

### Code Style

The project follows PEP 8 guidelines. Use `black` for code formatting and `flake8` for linting.

### Hot Reload

The development server supports hot reload. Any code changes to `main.py` will reflect immediately without needing to rebuild.

## Troubleshooting

### Common Docker Issues

1. **ModuleNotFoundError: No module named 'openai'**
   Solution options:
   - Install directly in container: `docker exec -it fastapi_app pip install openai`
   - Specify version in requirements.txt: `openai==1.3.7`
   - Rebuild containers: `docker-compose build --no-cache`

2. **Container fails to start**
   - Check logs: `docker-compose logs api`
   - Verify environment variables
   - Check disk space: `df -h`

For more detailed troubleshooting and deployment instructions, see [docs/deployment.md](docs/deployment.md).

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests
5. Submit a pull request

## License

[Your License Here]
