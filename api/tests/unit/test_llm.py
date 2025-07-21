# mypy: ignore-errors

"""Unit tests for LLM service."""
from typing import Any, Dict
from unittest.mock import MagicMock, patch
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))

import pytest
from pydantic import BaseModel

from api.agent_management.llm_service import LLMRequest, LLMResponse, LLMService
from api.agent_management.model_config import LLMModelConfig, ProviderType
from api.agent_management.models import StructuredLLMRequest


@pytest.fixture
def llm_service() -> LLMService:
    """Create a test LLM service instance."""
    config = LLMModelConfig(provider=ProviderType.OPENAI, model_name="gpt-3.5-turbo")
    return LLMService(config)


@pytest.mark.unit
class TestLLMService:
    """Test suite for LLMService."""

    def test_init(self, llm_service: LLMService) -> None:
        """Test LLMService initialization."""
        assert llm_service.config.provider == ProviderType.OPENAI
        assert llm_service.config.model_name == "gpt-3.5-turbo"

    def test_generate(self, llm_service: LLMService) -> None:
        """Test generate method (mock provider)."""
        llm_service._provider.generate = MagicMock(return_value=LLMResponse(content="ok", model="test", usage={}))  # type: ignore
        response = llm_service.generate(LLMRequest(user_prompt="Test"))
        assert response.content == "ok"

    def test_generate_structured(self, llm_service: LLMService) -> None:
        """Test generate_structured method with provider mocked via attribute."""
        expected_result = {"key": "value"}

        # Mock the provider's generate_structured directly
        llm_service._provider.generate_structured = MagicMock(return_value=expected_result)  # type: ignore

        class TestSchema(BaseModel):
            key: str = "value"

        request = StructuredLLMRequest[TestSchema](  # type: ignore
            user_prompt="Test prompt",
            output_schema=TestSchema,  # type: ignore
        )

        response: Dict[str, str] = llm_service.generate_structured(request)
        assert response == expected_result
        llm_service._provider.generate_structured.assert_called_once()  # type: ignore
