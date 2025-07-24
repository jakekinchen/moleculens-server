import pytest
from pydantic import BaseModel

import api.agent_management.providers.openai_provider as provider


class DummyModel(BaseModel):
    foo: str
    bar: int


class DummyMessage:
    def __init__(self, content):
        self.content = content


class DummyChoice:
    def __init__(self, message):
        self.message = message


class DummyCompletion:
    def __init__(self, content):
        self.choices = [DummyChoice(DummyMessage(content))]


def test_generate_structured_success(monkeypatch):
    # Arrange
    dummy_json = '{"foo":"hello","bar":123}'

    def fake_create(model, messages, response_format):
        return DummyCompletion(dummy_json)

    monkeypatch.setattr(provider.client.chat.completions, "create", fake_create)

    # Act
    result = provider.generate_structured(
        "Test prompt",
        DummyModel,
        model_name="test-model",
        system_prompt="Test system",
    )

    # Assert
    assert isinstance(result, DummyModel)
    assert result.foo == "hello"
    assert result.bar == 123


def test_generate_structured_no_content(monkeypatch):
    # Arrange: create returns no content
    def fake_create(model, messages, response_format):
        return DummyCompletion(None)

    monkeypatch.setattr(provider.client.chat.completions, "create", fake_create)

    # Act & Assert
    with pytest.raises(ValueError) as excinfo:
        provider.generate_structured("Test prompt", DummyModel)
    assert "No content received" in str(excinfo.value)
