from fastapi.testclient import TestClient

from api.main import app
from agent_management.llm_service import LLMService
from api.routers.prompt.routes import DiagramPlan, MoleculePlacement

client = TestClient(app)


def test_generate_molecule_diagram(monkeypatch):
    dummy_plan = DiagramPlan(
        plan="simple plan",
        molecule_list=[MoleculePlacement(molecule="water", x=0, y=0, width=100, height=100)],
        arrows=None,
    )

    def fake_structured(self, request):
        return dummy_plan

    def fake_layout(self, queries):
        return [
            {
                "atoms": [
                    {"element": "O", "x": 0.0, "y": 0.0},
                    {"element": "H", "x": 1.0, "y": 0.0},
                    {"element": "H", "x": 0.0, "y": 1.0},
                ],
                "bonds": [
                    {"start": 0, "end": 1, "order": 1},
                    {"start": 0, "end": 2, "order": 1},
                ],
                "name": "water",
                "cid": 1,
                "formula": "H2O",
                "box": queries[0]["box"],
                "query": "water",
            }
        ]

    monkeypatch.setattr(LLMService, "generate_structured", fake_structured)
    monkeypatch.setattr(
        "agent_management.agents.pubchem_agent.PubChemAgent.get_molecules_2d_layout",
        fake_layout,
    )

    response = client.post(
        "/prompt/generate-molecule-diagram/",
        json={"prompt": "a diagram"},
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["diagram_plan"]["plan"] == "simple plan"
    assert "<svg" in payload["diagram_image"]
