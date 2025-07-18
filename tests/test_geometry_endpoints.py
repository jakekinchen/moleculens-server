from fastapi.testclient import TestClient
from api.main import app

client = TestClient(app)


def test_cube_endpoint_returns_jsx():
    resp = client.get("/geometry/cube")
    assert resp.status_code == 200
    data = resp.json()
    assert "jsx" in data and "<Canvas" in data["jsx"]


def test_html_page_returns_html():
    resp = client.get("/geometry/html-test-page")
    assert resp.status_code == 200
    data = resp.json()
    assert "html" in data and "<html" in data["html"].lower() 