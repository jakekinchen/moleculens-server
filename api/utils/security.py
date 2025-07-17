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


def validate_commands(commands: list[str]):
    for cmd in commands:
        name = cmd.split()[0]
        if name not in ALLOWED_COMMANDS:
            raise ValueError(f"Command '{name}' not allowed")
