"""
Debug utilities for the agent management system.
"""

import os

# Debug flag - set to False to disable debug logging
DEBUG_PUBCHEM = True


def write_debug_file(filename: str, content: str) -> None:
    """Write debug content to a file if DEBUG_PUBCHEM is True."""
    if not DEBUG_PUBCHEM:
        return

    # Create debug directory if it doesn't exist
    current_file = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file)
    parent_dir = os.path.dirname(current_dir)
    workspace_dir = os.path.dirname(parent_dir)
    debug_dir = os.path.join(workspace_dir, "debug")

    print(f"[DEBUG] Current file: {current_file}")
    print(f"[DEBUG] Current directory: {current_dir}")
    print(f"[DEBUG] Parent directory: {parent_dir}")
    print(f"[DEBUG] Workspace directory: {workspace_dir}")
    print(f"[DEBUG] Debug directory: {debug_dir}")

    os.makedirs(debug_dir, exist_ok=True)

    # Write content to file (overwriting if exists)
    filepath = os.path.join(debug_dir, filename)
    try:
        with open(filepath, "w") as f:
            f.write(content)
        print(f"[DEBUG] Debug file written: {filepath}")
    except Exception as e:
        print(f"[ERROR] Error writing debug file {filepath}: {e}")

        # Try writing to a different location as a fallback
        fallback_dir = os.path.join(current_dir, "debug_fallback")
        os.makedirs(fallback_dir, exist_ok=True)
        fallback_path = os.path.join(fallback_dir, filename)
        try:
            with open(fallback_path, "w") as f:
                f.write(content)
            print(f"[DEBUG] Fallback debug file written: {fallback_path}")
        except Exception as e2:
            print(f"[ERROR] Error writing fallback debug file {fallback_path}: {e2}")
