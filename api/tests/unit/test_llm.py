# mypy: ignore-errors
"""Unit tests for LLM service."""
import os
import sys
from typing import Dict, Type, TypeVar
from unittest.mock import MagicMock

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
)

import pytest
from pydantic import BaseModel

from api.agent_management.llm_service import LLMService
from api.agent_management.model_config import LLMModelConfig, ProviderType


class TestSchema(BaseModel):
    key: str = "value"


T = TypeVar("T", bound=BaseModel)


@pytest.fixture
def llm_service() -> LLMService[TestSchema]:
    """Create a test LLM service instance."""
    config = LLMModelConfig(provider=ProviderType.OPENAI, model_name="gpt-3.5-turbo")
    return LLMService[TestSchema](config)


@pytest.mark.unit
class TestLLMService:
    """Test suite for LLMService."""

    def test_init(self, llm_service: LLMService[TestSchema]) -> None:
        """Test LLMService initialization."""
        assert llm_service.config.provider == ProviderType.OPENAI
        assert llm_service.config.model_name == "gpt-3.5-turbo"

    def test_generate_structured(self, llm_service: LLMService[TestSchema]) -> None:
        """Test generate_structured method with provider mocked."""
        expected_result = TestSchema(key="value")

        # Mock the provider's generate_structured directly
        llm_service._provider.generate_structured = MagicMock(
            return_value=expected_result
        )

        response = llm_service.generate_structured(
            user_prompt="Test prompt",
            response_model=TestSchema,
            system_prompt="Test system prompt",
        )
        assert response == expected_result
        llm_service._provider.generate_structured.assert_called_once_with(
            user_prompt="Test prompt",
            response_model=TestSchema,
            system_prompt="Test system prompt",
            model_name="gpt-3.5-turbo",
        )
