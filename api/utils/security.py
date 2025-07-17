ALLOWED_COMMANDS = {
    "fetch",
    "load",
    "hide",
    "show",
    "color",
    "bg_color",
    "orient",
    "turn",
    "rotate",
    "png",
    "save",
}


def _command_name(command: str) -> str:
    """Return the base command name from a PyMOL command string.

    Supports both `cmd.fetch("1abc")` and `fetch 1abc` syntaxes.
    """
    # Strip optional cmd. prefix
    if command.startswith("cmd."):
        command = command[4:]

    command = command.strip()

    # Handle parentheses style e.g. fetch("1abc")
    before_paren = command.split("(", 1)[0]
    name = before_paren.split()[0]  # Also handles whitespace style
    return name


def validate_commands(commands: list[str]) -> None:
    """Validate PyMOL commands against the allowed whitelist."""

    for cmd in commands:
        name = _command_name(cmd)
        if name not in ALLOWED_COMMANDS:
            raise ValueError(f"Command '{name}' not allowed")
