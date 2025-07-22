import importlib.util
import os
import sys


def test_imports():
    print("Python path:", sys.path)
    print("Current directory:", os.getcwd())
    print("Directory contents:", os.listdir())

    try:
        # Try to load the module directly
        spec = importlib.util.spec_from_file_location(
            "routes", "/app/api/routers/render/routes.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        print("Successfully loaded routes module")
        assert hasattr(module, "router")
        print("Router found in module")

    except Exception as e:
        print(f"Error: {e}")
        raise
