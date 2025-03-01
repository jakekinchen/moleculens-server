# Free Throw Agentic System

A modular system for generating and visualizing 3D animations using LLM-generated keyframes. This system uses OpenAI's GPT models to generate keyframe descriptions and converts them into interactive Three.js visualizations.

## Features

- LLM-powered keyframe generation
- Dynamic 3D animations using Three.js
- Modular architecture with support for multiple LLM providers
- Automatic keyframe validation and error handling
- Real-time animation updates based on keyframe timing

## Installation

1. Clone the repository
2. Install dependencies:
```bash
pip install -r requirements.txt
```
3. Set up your OpenAI API key:
```bash
export OPENAI_API_KEY='your-api-key-here'
```

## Usage

Run the main script:
```bash
python freethrow-agent.py
```

The script will:
1. Prompt you for what you'd like to visualize
2. Generate keyframes using the LLM
3. Create an HTML file with the visualization
4. Save the result as `output.html`

Open `output.html` in your web browser to view the visualization.

## Architecture

The system consists of several key components:

- `llm_service.py`: Handles LLM provider integration
- `keyframe_system.py`: Manages keyframe generation and validation
- `freethrow-agent.py`: Main orchestrator and visualization generator

## Supported Animation Types

The system currently supports these animation types:
- Rotation
- Color changes
- Scaling

## Requirements

- Python 3.7+
- OpenAI API key
- Modern web browser with WebGL support

## License

MIT License 