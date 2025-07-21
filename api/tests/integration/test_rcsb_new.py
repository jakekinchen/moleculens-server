import pytest
from fastapi.testclient import TestClient

from agent_management.agents.rcsb_agent import RCSBAgent


@pytest.mark.integration
def test_pairwise_alignment_endpoint(test_client: TestClient):
    original = RCSBAgent.fetch_pairwise_alignment
    def fake(self, id_one: str, id_two: str) -> dict:
        assert id_one == "1AAA"
        assert id_two == "2BBB"
        return {"score": 42}
    RCSBAgent.fetch_pairwise_alignment = fake
    try:
        resp = test_client.post("/rcsb/align/", json={"identifier1": "1AAA", "identifier2": "2BBB"})
        assert resp.status_code == 200
        assert resp.json()["metadata"]["score"] == 42
    finally:
        RCSBAgent.fetch_pairwise_alignment = original


@pytest.mark.integration
def test_group_entries_endpoint(test_client: TestClient):
    original = RCSBAgent.fetch_group_entries
    def fake(self, group_id: str) -> dict:
        assert group_id == "G123"
        return {"entries": ["1AAA", "2BBB"]}
    RCSBAgent.fetch_group_entries = fake
    try:
        resp = test_client.get("/rcsb/group/G123")
        assert resp.status_code == 200
        assert resp.json()["metadata"]["entries"] == ["1AAA", "2BBB"]
    finally:
        RCSBAgent.fetch_group_entries = original


@pytest.mark.integration
def test_feature_annotations_endpoint(test_client: TestClient):
    original = RCSBAgent.fetch_feature_annotations
    def fake(self, identifier: str) -> dict:
        assert identifier == "1AAA"
        return {"features": []}
    RCSBAgent.fetch_feature_annotations = fake
    try:
        resp = test_client.get("/rcsb/feature-annotations/1AAA")
        assert resp.status_code == 200
        assert resp.json()["annotations"]["features"] == []
    finally:
        RCSBAgent.fetch_feature_annotations = original
