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
    "label",
    "select",  # Essential for creating selections
    "set",  # Essential for setting properties like transparency
    "zoom",  # Essential for focusing on specific regions
    "surface",  # For surface representations
    # NEW: Advanced rendering commands
    "ray",  # Ray-tracing
    "antialias",  # Edge smoothing
    "depth_cue",  # Depth emphasis
    "clip",  # Viewport clipping
    "field_of_view",  # Camera controls
    "spectrum",  # Color gradients
    "ramp_new",  # Custom color ramps
    "distance",  # Measurement objects
    "angle",  # Angle measurements
    "dihedral",  # Dihedral measurements
    "cgo_arrow",  # Custom graphics objects
    "h_add",  # Hydrogen bond detection
    "find_pairs",  # Interaction detection
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
