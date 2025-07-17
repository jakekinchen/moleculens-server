from api.utils.security import validate_commands
import pytest


def test_allowed_commands():
    validate_commands(["fetch('1abc')", "hide everything", "show sticks"])


def test_disallowed_command():
    with pytest.raises(ValueError):
        validate_commands(["delete all"])
