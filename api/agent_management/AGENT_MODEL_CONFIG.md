# Agent-Model Configuration System

This document describes the agent-model configuration system that maps specific LLM models to different agents in the visualization pipeline.

## Architecture

The agent-model configuration system follows these principles:

1. **Separation of concerns**: Each agent performs a specific task in the visualization pipeline.
2. **Model Registry pattern**: Models are registered centrally and can be accessed by name.
3. **Agent-Model mapping**: Each agent is configured to use the optimal model for its specific task.
4. **Factory pattern**: Agents are created through factory methods that handle model selection.

## Key Components

### Model Registry

The `ModelRegistry` class in `models.py` serves as a central repository for all model configurations:

- Register models with `ModelRegistry.register(name, model_cls, factory)`
- Retrieve model classes with `ModelRegistry.get_model(name)`
- Create model instances with `ModelRegistry.create_instance(name, **kwargs)`
- List all registered models with `ModelRegistry.list_models()`

### Model Configuration

The `model_config.py` file defines:

- `ModelInfo` class for storing model metadata
- `ModelCategory` enum for categorizing models by capability (GENERAL, CODE, REASONING, VISION)
- Registration of all supported models with their capabilities
- Helper functions for working with models

### Agent-Model Configuration

The `agent_model_config.py` file defines:

- `AgentType` enum for all agent types in the pipeline
- `AgentModelConfig` class for defining which model an agent should use
- Default mappings of agents to their preferred models
- Helper functions for retrieving model configurations for agents

### Agent Factory

The `agent_factory.py` file implements:

- Factory methods for creating each type of agent with the appropriate model
- A method for creating all agents at once with consistent model selection
- Support for overriding model selection when needed

## Default Agent-Model Mappings

| Agent Type        | Preferred Model         | Fallback Models                      | Required Categories |
|-------------------|-------------------------|------------------------------------|-------------------|
| DOMAIN_VALIDATOR  | o3-mini                 | gpt-4o, claude-3-5-sonnet-latest   | GENERAL           |
| SCRIPT            | gpt-4o                  | claude-3-7-sonnet-latest, gpt-4.5-preview | REASONING  |
| ORCHESTRATION     | gpt-4o                  | claude-3-7-sonnet-latest, gpt-4.5-preview | REASONING  |
| GEOMETRY          | claude-3-7-sonnet-latest| gpt-4.5-preview, qwen-2.5-coder-32b | CODE             |
| ANIMATION         | claude-3-7-sonnet-latest| gpt-4.5-preview, gpt-4o             | CODE             |
| CAPTION           | o3-mini                 | gpt-4o, llama3-70b-8192             | GENERAL           |
| AGGREGATOR        | o3-mini                 | gpt-4o, llama3-70b-8192             | GENERAL           |

## Usage

### Basic Usage

```python
from agent_management.agent_factory import AgentFactory

# Create an agent with its default model
script_agent = AgentFactory.create_script_agent()

# Create an agent with a specific model
geometry_agent = AgentFactory.create_geometry_agent("claude-3-7-sonnet-latest")
```

### Model Selection Logic

For each agent, the model selection follows this priority:

1. Explicitly specified model (override_model parameter)
2. Preferred model for that agent type
3. Fallback models if the preferred model is not available
4. Any model that matches the required categories for the agent

### API Endpoints

The following endpoints support agent-model configuration:

- `GET /models/` - List all available models with their capabilities
- `GET /agent-models/` - List all agent types with their preferred models
- `POST /prompt/` - Process a prompt using the specified model

#### Query Parameters

Most endpoints in the API now support the following query parameters for model selection:

- `model`: Specific model name to use (e.g., "gpt-4o", "claude-3-7-sonnet-latest")
- `preferred_model_category`: Model category preference (GENERAL, CODE, REASONING, VISION)

#### Supported Endpoints

These endpoints now support model specification parameters:

- `/prompt/generate-geometry/` - Controls models for geometry generation
- `/prompt/generate-script/` - Controls models for script generation
- `/prompt/generate-orchestration/` - Controls models for orchestration plan generation
- `/prompt/generate-animation/` - Controls models for animation code generation
- `/prompt/generate-geometry-for-plan/` - Controls models for geometry generation from plans
- `/prompt/process/` - Controls models for the entire pipeline
- `/prompt/` - Controls models for the entire pipeline

#### Example API Usage

```javascript
// Frontend example: Override model for geometry generation
fetch('/prompt/generate-geometry/', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    prompt: "Create a DNA double helix",
    model: "claude-3-7-sonnet-latest"  // Use Claude instead of default
  })
});

// Frontend example: Use specific model for animation generation
fetch('/prompt/generate-animation/', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    // ... other parameters
    model: "gpt-4o"  // Override to use GPT-4o
  })
});
```

## Configuration

To modify agent-model mappings:

1. Edit `agent_model_config.py` to update the `DEFAULT_AGENT_MODELS` list
2. Register new models in `model_config.py` if needed

## Adding New Models

To add a new model:

1. Add the model to `model_config.py` with its capabilities
2. Update agent model preferences in `agent_model_config.py` if needed

## Adding New Agents

To add a new agent type:

1. Implement the agent class
2. Add the agent type to `AgentType` enum
3. Add the agent to `DEFAULT_AGENT_MODELS`
4. Add a factory method in `AgentFactory`
