#!/usr/bin/env python3
"""Test runner script for MoleculeLens Server API tests."""

import os
import subprocess
import sys
from pathlib import Path


def main():
    """Run all tests and generate coverage report."""
    # Set up environment
    os.environ["PYTHONPATH"] = str(Path(__file__).parent.parent / "api")
    os.environ["ENVIRONMENT"] = "development"

    # Change to project root
    project_root = Path(__file__).parent.parent
    os.chdir(project_root)

    print("🧪 Running MoleculeLens Server API Tests")
    print("=" * 50)

    # Run tests with pytest
    cmd = [sys.executable, "-m", "pytest", "tests/", "-v", "--tb=short", "--color=yes", "--durations=10"]

    try:
        subprocess.run(cmd, check=True)
        print("\n✅ All tests passed!")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Tests failed with exit code {e.returncode}")
        return e.returncode
    except KeyboardInterrupt:
        print("\n⚠️  Tests interrupted by user")
        return 1


if __name__ == "__main__":
    sys.exit(main())
