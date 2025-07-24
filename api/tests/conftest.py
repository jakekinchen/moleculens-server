"""Test configuration and shared fixtures."""

import os
import sys
import types
from enum import Enum
from typing import Any, Dict, Generator, Generic, Optional, Type, TypeVar

import pytest
from pydantic import BaseModel

from fastapi import FastAPI
from fastapi.testclient import TestClient

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

import importlib.util
from pathlib import Path

# Load the RCSB agent module directly and expose it under the expected package path
_agent_path = (
    Path(__file__).resolve().parents[1]
    / "agent_management"
    / "agents"
    / "rcsb_agent.py"
)
spec_agent = importlib.util.spec_from_file_location(
    "agent_management.agents.rcsb_agent", _agent_path
)
if spec_agent is None:
    raise ImportError("Could not create module spec")
rcsb_agent_mod = importlib.util.module_from_spec(spec_agent)
assert spec_agent and spec_agent.loader
spec_agent.loader.exec_module(rcsb_agent_mod)

# Minimal model_config stub required by models and llm_service
model_config_stub = types.ModuleType("api.agent_management.model_config")


class _ProviderType(str, Enum):
    OPENAI = "openai"


class _LLMModelConfig(BaseModel):  # type: ignore[misc]
    provider: _ProviderType
    model_name: str


model_config_stub.ProviderType = _ProviderType  # type: ignore[attr-defined]
model_config_stub.LLMModelConfig = _LLMModelConfig  # type: ignore[attr-defined]
sys.modules["api.agent_management.model_config"] = model_config_stub

agent_pkg = types.ModuleType("agent_management.agents")
agent_pkg.rcsb_agent = rcsb_agent_mod  # type: ignore[attr-defined]
agent_root = types.ModuleType("agent_management")
agent_root.agents = agent_pkg  # type: ignore[attr-defined]
# sys.modules["agent_management"] = agent_root
# sys.modules["agent_management.agents"] = agent_pkg
# sys.modules["agent_management.agents.rcsb_agent"] = rcsb_agent_mod
api_mod = types.ModuleType("api")
# api_mod.agent_management = agent_root  # type: ignore[attr-defined]
# sys.modules["api"] = api_mod
# sys.modules["api.agent_management"] = agent_root
api_mod.__path__ = [str(Path(__file__).resolve().parents[2])]
# sys.modules["api.agent_management"] = agent_root
# Allow loading of real subpackages under api.agent_management
agent_root.__path__ = [str(Path(__file__).resolve().parents[1] / "agent_management")]
sys.modules["api.agent_management.agents"] = agent_pkg
sys.modules["api.agent_management.agents.rcsb_agent"] = rcsb_agent_mod

# Load real providers package and modules
providers_init = (
    Path(__file__).resolve().parents[1]
    / "agent_management"
    / "providers"
    / "__init__.py"
)
spec_providers = importlib.util.spec_from_file_location(
    "api.agent_management.providers", providers_init
)
assert spec_providers is not None
assert spec_providers.loader is not None
providers_pkg = importlib.util.module_from_spec(spec_providers)
spec_providers.loader.exec_module(providers_pkg)
sys.modules["api.agent_management.providers"] = providers_pkg

openai_file = (
    Path(__file__).resolve().parents[1]
    / "agent_management"
    / "providers"
    / "openai_provider.py"
)
spec_openai = importlib.util.spec_from_file_location(
    "api.agent_management.providers.openai_provider", openai_file
)
assert spec_openai is not None
assert spec_openai.loader is not None
openai_mod = importlib.util.module_from_spec(spec_openai)
spec_openai.loader.exec_module(openai_mod)
sys.modules["api.agent_management.providers.openai_provider"] = openai_mod

# Provide lightweight stubs for models and the LLM service used in unit tests
models_stub = types.ModuleType("agent_management.models")


class LLMResponse(BaseModel):
    content: str = ""
    model: str = ""
    usage: Dict[str, int] = {}


T = TypeVar("T")


class StructuredLLMRequest(BaseModel, Generic[T]):
    user_prompt: str
    output_schema: Type[T]
    llm_config: Optional[_LLMModelConfig] = None


models_stub.LLMResponse = LLMResponse  # type: ignore[attr-defined]
models_stub.StructuredLLMRequest = StructuredLLMRequest  # type: ignore[attr-defined]
sys.modules["agent_management.models"] = models_stub
agent_root.models = models_stub  # type: ignore[attr-defined]
sys.modules["api.agent_management.models"] = models_stub


llm_stub = types.ModuleType("agent_management.llm_service")


class _DummyProvider:
    def generate(self, request: Any) -> LLMResponse:
        return LLMResponse(content="ok", model="test", usage={})

    def generate_structured(self, request: Any) -> Dict[str, Any]:
        return {"key": "value"}


class LLMService:
    def __init__(self, config: _LLMModelConfig):
        self.config = config
        self._provider = _DummyProvider()

    def generate(self, request: "LLMRequest") -> LLMResponse:
        return self._provider.generate(request)

    def generate_structured(self, request: Any) -> Dict[str, Any]:
        return self._provider.generate_structured(request)


class LLMRequest(BaseModel):
    user_prompt: str


llm_stub.LLMService = LLMService  # type: ignore[attr-defined]
llm_stub.LLMRequest = LLMRequest  # type: ignore[attr-defined]
llm_stub.LLMResponse = LLMResponse  # type: ignore[attr-defined]
sys.modules["agent_management.llm_service"] = llm_stub
agent_root.llm_service = llm_stub  # type: ignore[attr-defined]
sys.modules["api.agent_management.llm_service"] = llm_stub

# Stub heavy optional dependencies before loading router module
chem_mod = types.ModuleType("rdkit.Chem")
chem_mod.Fragments = object  # type: ignore[attr-defined]
chem_mod.Descriptors = object  # type: ignore[attr-defined]
chem_mod.AllChem = object  # type: ignore[attr-defined]
rdkit_mod = types.ModuleType("rdkit")
rdkit_mod.Chem = chem_mod  # type: ignore[attr-defined]
sys.modules["rdkit"] = rdkit_mod
sys.modules["rdkit.Chem"] = chem_mod

openai_mod = types.ModuleType("openai")
openai_mod.OpenAI = object  # type: ignore[attr-defined]
types_mod = types.ModuleType("openai.types")
chat_mod = types.ModuleType("openai.types.chat")
chat_completion_mod = types.ModuleType("openai.types.chat.chat_completion")
setattr(types_mod, "Completion", object)
setattr(chat_mod, "ChatCompletion", object)
setattr(chat_mod, "ChatCompletionMessage", object)
setattr(chat_mod, "ChatCompletionMessageParam", dict)
setattr(chat_completion_mod, "Choice", object)
sys.modules["openai"] = openai_mod
sys.modules["openai.types"] = types_mod
sys.modules["openai.types.chat"] = chat_mod
sys.modules["openai.types.chat.chat_completion"] = chat_completion_mod
# Stub openai.types.response_format for providers
response_format_mod = types.ModuleType("openai.types.response_format")


class ResponseFormatJSONSchema:
    def __init__(self, *args, **kwargs):
        pass


setattr(response_format_mod, "ResponseFormatJSONSchema", ResponseFormatJSONSchema)
sys.modules["openai.types.response_format"] = response_format_mod

import importlib.util
from pathlib import Path

_routes_path = Path(__file__).resolve().parents[1] / "routers" / "rcsb" / "routes.py"
spec = importlib.util.spec_from_file_location("rcsb_routes", _routes_path)
if spec is None:
    raise ImportError("Could not create module spec")
rcsb_module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(rcsb_module)
rcsb_router = rcsb_module.router

# Create a minimal FastAPI app for unit tests to avoid importing the full server
app = FastAPI()
app.include_router(rcsb_router)


@app.post("/prompt/")
def submit_prompt(prompt: Dict[str, str]) -> Dict[str, str]:
    return {"job_id": "job"}


@app.get("/prompt/process/{job_id}")
def job_status(job_id: str) -> Dict[str, str]:
    return {"status": "completed"}


@app.post("/prompt/fetch-molecule-data/")
def fetch_molecule_data(query: Dict[str, str]) -> Dict[str, Any]:
    return {"molecule_data": {}}


@app.get("/health")
def health() -> Dict[str, str]:
    return {"status": "healthy"}


# Set test environment
os.environ["ENVIRONMENT"] = "test"

# Stub the heavy routers package so importing api modules doesn't pull full FastAPI routes
sys.modules.setdefault("api.routers", types.ModuleType("api.routers"))
sys.modules.setdefault("routers", sys.modules["api.routers"])


@pytest.fixture
def test_client() -> Generator[TestClient, None, None]:
    """Create a test client for the FastAPI app."""
    with TestClient(app) as client:
        yield client


@pytest.fixture
def mock_openai_response() -> Dict[str, Any]:
    """Mock response from OpenAI API."""
    return {
        "choices": [
            {
                "message": {"content": "Test response", "role": "assistant"},
                "finish_reason": "stop",
                "index": 0,
            }
        ],
        "model": "gpt-3.5-turbo",
        "object": "chat.completion",
        "usage": {"completion_tokens": 10, "prompt_tokens": 20, "total_tokens": 30},
    }


@pytest.fixture
def mock_pubchem_response() -> Dict[str, Any]:
    """Mock response from PubChem API."""
    return {
        "PC_Compounds": [
            {
                "id": {"id": {"cid": 962}},
                "props": [{"urn": {"label": "IUPAC Name"}, "value": {"sval": "water"}}],
            }
        ]
    }


@pytest.fixture
def test_env_vars(monkeypatch: pytest.MonkeyPatch) -> None:
    """Set up test environment variables."""
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setenv("ENVIRONMENT", "test")
    monkeypatch.setenv("REDIS_HOST", "localhost")
    monkeypatch.setenv("REDIS_PORT", "6379")
