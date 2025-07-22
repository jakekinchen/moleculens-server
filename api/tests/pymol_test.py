import os

os.environ["PYMOL_QUIET"] = "1"
os.environ["PYMOL_HEADLESS"] = "1"

try:
    import pymol

    print("PyMOL imported successfully!")
    pymol.finish_launching(["pymol", "-cq"])  # Launch in headless mode
    print("PyMOL launched successfully!")

    # Try to create a simple molecule
    pymol.cmd.fragment("ala")
    print("Created alanine fragment successfully!")

    # Clean up
    pymol.cmd.quit()
except ImportError as e:
    print(f"Failed to import PyMOL: {e}")
except Exception as e:
    print(f"Error while using PyMOL: {e}")
