import json
from types import SimpleNamespace
import requests
import pytest

from agent_management.agents.pubchem_agent import PubChemAgent
from agent_management.llm_service import LLMService, LLMModelConfig, ProviderType


def make_agent():
    config = LLMModelConfig(provider=ProviderType.OPENAI, model_name="gpt-3.5-turbo")
    llm_service = LLMService(config)
    return PubChemAgent(llm_service)


@pytest.fixture()
def sample_sdf():
    path = "api/compound_data/cid_4140276/compound_details.json"
    with open(path) as f:
        data = json.load(f)
    return data["structure"]["sdf_2d"]


def test_get_molecule_2d_info(monkeypatch, sample_sdf):
    agent = make_agent()

    class FakeCompound(SimpleNamespace):
        pass

    fake_compound = FakeCompound(
        cid=4140276,
        iupac_name="Water",
        molecular_formula="H2O",
    )

    def fake_interpret(query: str):
        return "water"

    def fake_search(query: str):
        return [fake_compound]

    class FakeResponse:
        status_code = 200
        text = sample_sdf

    def fake_get(url, timeout=30):
        return FakeResponse()

    monkeypatch.setattr(agent, "interpret_user_query", fake_interpret)
    monkeypatch.setattr(agent, "_search_with_fallbacks", fake_search)
    monkeypatch.setattr(requests, "get", fake_get)

    result = agent.get_molecule_2d_info("water")
    assert result["cid"] == 4140276
    assert len(result["atoms"]) > 0
    assert len(result["bonds"]) > 0
