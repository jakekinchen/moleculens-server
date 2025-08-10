# Clean module exports for refactored pymol package
__all__ = [
    "services",
    "clients",
    "visualizers",
    "pymol_security",
    "pymol_templates",
    "pymol_translator",
    "scene_packager",
    "scene_spec",
]

# Import system PyMOL for 3D rendering
def get_system_pymol():
    """Import and return the system PyMOL module."""
    import sys
    import importlib
    
    # Temporarily remove current directory from path
    original_path = sys.path[:]
    sys.path = [p for p in sys.path if not p.endswith('/api')]
    
    try:
        # Force reload of pymol module
        if 'pymol' in sys.modules:
            del sys.modules['pymol']
        if 'pymol.cmd' in sys.modules:
            del sys.modules['pymol.cmd']
            
        import pymol as system_pymol
        return system_pymol
    finally:
        sys.path = original_path
