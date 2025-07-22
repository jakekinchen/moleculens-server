import importlib.util
import os
import sys
from pathlib import Path


def test_imports():
    print("Python path:", sys.path)
    print("Current directory:", os.getcwd())
    print("Directory contents:", os.listdir())

    try:
        routes_path = (
            Path(__file__).resolve().parents[1] / "routers" / "render" / "routes.py"
        )
        spec = importlib.util.spec_from_file_location("routes", routes_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        print("Successfully loaded routes module")
        assert hasattr(module, "router")
        print("Router found in module")

    except Exception as e:
        print(f"Error: {e}")
        raise
