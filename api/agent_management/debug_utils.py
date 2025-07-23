"""Debug utilities for agent management."""

import os
from pathlib import Path

# Debug flags
DEBUG_PUBCHEM = os.getenv("DEBUG_PUBCHEM", "false").lower() == "true"


def write_debug_file(
    filename: str, content: str, debug_dir: str = "debug_output"
) -> None:
    """Write debug content to a file.

    Parameters
    ----------
    filename : str
        Name of the debug file
    content : str
        Content to write
    debug_dir : str
        Directory to write debug files to
    """
    if not DEBUG_PUBCHEM:
        return

    debug_path = Path(debug_dir)
    debug_path.mkdir(exist_ok=True)

    file_path = debug_path / filename
    with open(file_path, "w") as f:
        f.write(content)
