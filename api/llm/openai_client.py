"""Centralized OpenAI client helper to ensure consistent initialization."""

import os
import threading
from typing import Optional

from openai import OpenAI

_client_instance: Optional[OpenAI] = None
_lock = threading.Lock()


def get_client(api_key: Optional[str] = None) -> OpenAI:
    """Thread-safe singleton OpenAI client.

    Args:
        api_key: Optional API key. If not provided, will use OPENAI_API_KEY env var.

    Returns:
        OpenAI client instance

    Raises:
        ValueError: If no API key is available
    """
    global _client_instance

    if _client_instance is None:
        with _lock:
            if _client_instance is None:  # double-check
                effective_api_key = api_key or os.environ.get("OPENAI_API_KEY")
                if not effective_api_key:
                    raise ValueError("OpenAI API key is required. Set OPENAI_API_KEY env var or pass api_key.")
                _client_instance = OpenAI(api_key=effective_api_key)

    return _client_instance


def reset_client() -> None:
    """Reset the singleton client instance. Useful for testing."""
    global _client_instance
    with _lock:
        _client_instance = None
