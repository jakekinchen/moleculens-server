from fastapi.testclient import TestClient

from api.main import app

# ---------------------------------------------------------------------------
# Stub requests.get in the patch module so no real network traffic occurs
# ---------------------------------------------------------------------------
sample_json = {
    "coordset": [
        {"seq_id": 1, "cartn_x": 1.1, "cartn_y": 2.2, "cartn_z": 3.3},
        {"seq_id": 2, "cartn_x": 4.4, "cartn_y": 5.5, "cartn_z": 6.6},
    ]
}


class _MockResponse:
    def __init__(self, json_data, status_code=200):
        self._json_data = json_data
        self.status_code = status_code

    def json(self):
        return self._json_data


def _mock_get(url, timeout=30):
    assert url.endswith("/sequence/coordinates/1ABC")
    return _MockResponse(sample_json)


def test_sequence_coordinates_endpoint(monkeypatch):
    # Patch the requests.get used inside the monkey‑patched method
    monkeypatch.setattr(
        "api.agent_management.agents.rcsb_sequence_patch.requests.get",
        _mock_get,
    )

    client = TestClient(app)
    resp = client.get("/rcsb/sequence-coordinates/1ABC")
    assert resp.status_code == 200
    data = resp.json()

    assert data["1"] == [1.1, 2.2, 3.3]
    assert data["2"] == [4.4, 5.5, 6.6]
