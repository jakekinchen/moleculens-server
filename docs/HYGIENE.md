# Code Hygiene & Type Safety Guide

This document outlines the code quality standards, type hygiene practices, and tooling used in the Moleculens Server codebase.

## Overview

We use **Ruff** as our single tool for linting, formatting, and type hygiene enforcement. This replaces the traditional Black + isort + Flake8 trio with a single, blazing-fast Rust-powered tool that's 30-120x faster.

## Tools & Configuration

### Ruff Configuration

Our Ruff configuration is defined in `pyproject.toml`:

```toml
[tool.ruff]
# Universal settings
target-version = "py39"              # Match project Python version
line-length = 120                    # Match former Black config
extend-exclude = ["migrations", "build", ".venv", "__pycache__"]
fix = true                           # Auto-apply fixes on save / CI

[tool.ruff.lint]
select = ["E", "F", "I", "UP", "N", "B"]  # PEP8, flakes, isort, pyupgrade, pep8-naming, bugbear
extend-select = [
    "TCH",   # type-checking completeness
    "PYI",   # stub-style purity
    "UP",    # pyupgrade (already on, but explicit)
]
ignore = ["E501"]                         # Long lines already enforced by formatter
extend-fixable = ["TCH", "PYI"]          # Allow auto-fixing type issues

[tool.ruff.format]                        # Black parity
indent-style = "space"
quote-style = "double"
docstring-code-format = true
```

### Dependency Management

We use **uv** for fast, modern Python dependency management:

```bash
# Install dependencies
uv sync

# Add new dependency
uv add package-name

# Add development dependency
uv add --group lint ruff
```

## Type Hygiene Standards

### Modern Type Annotations

We enforce modern Python type annotations with these rules:

#### ✅ DO: Use modern built-in types
```python
from __future__ import annotations

def process_data(items: list[dict[str, Any]]) -> tuple[str, int]:
    """Use built-in types instead of typing module equivalents."""
    pass
```

#### ❌ DON'T: Use deprecated typing module types
```python
from typing import List, Dict, Tuple  # Deprecated in Python 3.9+

def process_data(items: List[Dict[str, Any]]) -> Tuple[str, int]:
    pass
```

#### ✅ DO: Use TYPE_CHECKING for runtime imports
```python
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from openai.types.chat import ChatCompletionMessageParam

def generate_response(messages: list[ChatCompletionMessageParam]) -> str:
    """Import only needed for type checking, not runtime."""
    pass
```

#### ✅ DO: Use Optional for Python 3.9 compatibility
```python
from typing import Optional

# Python 3.9 compatible
def process_name(name: Optional[str] = None) -> str:
    pass

# Python 3.10+ only (avoid for now)
def process_name(name: str | None = None) -> str:
    pass
```

### Future Annotations

All new files should include:
```python
from __future__ import annotations
```

This enables:
- Forward references without quotes
- Better performance (deferred evaluation)
- Modern type syntax compatibility

## Commands & Workflows

### Daily Development

```bash
# Format code (replaces black + isort)
uv run ruff format .

# Lint and auto-fix issues
uv run ruff check . --fix

# Check specific rule categories
uv run ruff check . --select TCH,PYI,UP --fix
```

### Pre-commit Integration

We use pre-commit hooks that automatically run Ruff:

```yaml
# .pre-commit-config.yaml
- repo: https://github.com/astral-sh/ruff-pre-commit
  rev: v0.12.5
  hooks:
  - id: ruff
    args: [--fix]
    stages: [commit]
  - id: ruff-format
    stages: [commit]
```

Install pre-commit hooks:
```bash
pip install pre-commit
pre-commit install
```

### CI/CD Integration

Our GitHub Actions workflow enforces code quality:

```yaml
# .github/workflows/lint.yml
- name: Type hygiene check
  run: |
    uv sync --group lint
    uv run ruff check . --select TCH,PYI,UP --output-format=github
```

## Rule Categories Explained

### Core Rules (Always Enabled)
- **E, W**: PEP 8 style violations
- **F**: Pyflakes (undefined names, unused imports)
- **I**: Import sorting (isort replacement)
- **UP**: Pyupgrade (modernize Python syntax)
- **N**: PEP 8 naming conventions
- **B**: Bugbear (common Python gotchas)

### Type Hygiene Rules (Strict Mode)
- **TCH**: Type-checking completeness
  - Move runtime-only imports to `TYPE_CHECKING` blocks
  - Detect missing type annotations
- **PYI**: Stub-style purity
  - Enforce clean type stub patterns
  - Detect unnecessary runtime code in type definitions
- **UP**: Pyupgrade (enhanced)
  - Replace `typing.List` → `list`
  - Replace `typing.Dict` → `dict`
  - Replace `typing.Tuple` → `tuple`
  - Replace `typing.Optional[X]` → `X | None` (Python 3.10+)

## Migration Guide

### From Black + isort + Flake8

If you have muscle memory for the old tools, use these aliases:

```bash
# Add to ~/.zshrc or ~/.bashrc
alias black="uv run ruff format"
alias flake8="uv run ruff check"
alias isort="echo 'Import sorting is now handled by ruff format'"
```

### Fixing Type Issues

Common patterns for fixing type hygiene violations:

#### 1. Modernize type imports
```python
# Before
from typing import List, Dict, Optional

def process(data: List[Dict[str, Optional[str]]]) -> None:
    pass

# After
from __future__ import annotations
from typing import Optional

def process(data: list[dict[str, Optional[str]]]) -> None:
    pass
```

#### 2. Move runtime imports to TYPE_CHECKING
```python
# Before
from expensive_module import SomeClass

def func(obj: SomeClass) -> None:
    pass

# After
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from expensive_module import SomeClass

def func(obj: SomeClass) -> None:
    pass
```

#### 3. Add missing type annotations
```python
# Before
def calculate_total(items):
    return sum(item.price for item in items)

# After
def calculate_total(items: list[Item]) -> float:
    return sum(item.price for item in items)
```

## Architecture-Specific Guidelines

### Pydantic Models

```python
from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, Field

class MoleculeData(BaseModel):
    """Use modern type annotations in Pydantic models."""

    name: str
    formula: str
    atoms: list[dict[str, Any]] = Field(default_factory=list)
    bonds: Optional[list[dict[str, Any]]] = None

    # Use Field(default_factory=...) for mutable defaults
    metadata: dict[str, Any] = Field(default_factory=dict)
```

### FastAPI Routes

```python
from __future__ import annotations
from typing import TYPE_CHECKING
from fastapi import APIRouter

if TYPE_CHECKING:
    from pydantic import BaseModel

router = APIRouter()

@router.post("/molecules/")
async def create_molecule(data: MoleculeData) -> dict[str, Any]:
    """Type annotations help with OpenAPI generation."""
    return {"status": "created", "id": data.name}
```

### Service Classes

```python
from __future__ import annotations
from typing import TYPE_CHECKING, Optional
import logging

if TYPE_CHECKING:
    from external_lib import SomeType

logger = logging.getLogger(__name__)

class MoleculeService:
    """Service classes should have comprehensive type annotations."""

    def __init__(self, timeout: float = 30.0) -> None:
        self.timeout = timeout

    async def fetch_data(self, query: str) -> Optional[dict[str, Any]]:
        """Return type hints help with IDE support."""
        try:
            # Implementation here
            return {"data": query}
        except Exception as e:
            logger.error(f"Failed to fetch data: {e}")
            return None
```

## Current Status

As of the latest migration:

- ✅ **Ruff configured** with modern type hygiene rules
- ✅ **Pre-commit hooks** updated to use Ruff
- ✅ **CI pipeline** enforces type hygiene
- ✅ **Key files modernized** (diagram models, LLM provider)
- 📊 **~900 type issues remaining** (down from 1000+)

### Progress Tracking

Check current type hygiene status:
```bash
# Count remaining issues
uv run ruff check . --select TCH,PYI,UP --exit-zero | wc -l

# See specific issues
uv run ruff check . --select TCH,PYI,UP
```

## Best Practices

### 1. Incremental Improvement
- Fix type issues in files you're already modifying
- Don't create massive "type fixing" PRs
- Focus on high-impact files first (models, services, APIs)

### 2. Python Version Compatibility
- Use `Optional[X]` instead of `X | None` for Python 3.9 compatibility
- Use `list[X]` instead of `List[X]` (supported in Python 3.9+)
- Always include `from __future__ import annotations`

### 3. Performance Considerations
- Use `TYPE_CHECKING` blocks for expensive imports
- Defer type evaluation with future annotations
- Avoid runtime type checking unless necessary

### 4. Documentation
- Type annotations serve as documentation
- Use descriptive parameter names with types
- Add docstrings for complex type relationships

## Troubleshooting

### Common Issues

#### "unsupported operand type(s) for |"
```python
# Problem: Using Python 3.10+ syntax on Python 3.9
def func(x: str | None) -> None:
    pass

# Solution: Use Optional for Python 3.9 compatibility
from typing import Optional

def func(x: Optional[str]) -> None:
    pass
```

#### "Move import into TYPE_CHECKING block"
```python
# Problem: Runtime import only used for typing
from expensive_module import SomeClass

def func(obj: SomeClass) -> None:
    pass

# Solution: Move to TYPE_CHECKING
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from expensive_module import SomeClass

def func(obj: SomeClass) -> None:
    pass
```

#### Ruff not fixing issues automatically
```bash
# Some issues require manual fixes
uv run ruff check . --select TCH,PYI,UP --fix

# Check what can't be auto-fixed
uv run ruff check . --select TCH,PYI,UP --diff
```

### Getting Help

1. **Check Ruff documentation**: https://docs.astral.sh/ruff/
2. **Rule explanations**: `uv run ruff rule TCH001`
3. **Project-specific questions**: Check existing code patterns in `/api/diagram/models.py` and `/api/llm/openai_provider.py`

## Contributing

When contributing to this codebase:

1. **Run pre-commit hooks**: `pre-commit run --all-files`
2. **Check type hygiene**: `uv run ruff check . --select TCH,PYI,UP`
3. **Format code**: `uv run ruff format .`
4. **Follow patterns**: Use existing modernized files as examples
5. **Test imports**: Ensure your changes don't break Python 3.9 compatibility

Remember: The goal is gradual, sustainable improvement of code quality and type safety! 🚀
