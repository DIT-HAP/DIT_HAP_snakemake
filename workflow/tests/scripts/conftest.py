"""Fixtures for testing CLI scripts that are not importable as modules.

workflow/scripts/figures/ has no __init__.py, and workflow/src/figures.py shadows
the directory name on sys.path, so these scripts can only be loaded by file path.
"""

# =============================================================================
# IMPORTS
# =============================================================================
import importlib.util
from collections.abc import Callable
from pathlib import Path
from types import ModuleType

import pytest

# =============================================================================
# CONSTANTS
# =============================================================================
FIGURES_DIR = Path(__file__).resolve().parents[3] / "workflow" / "scripts" / "figures"

# =============================================================================
# FIXTURES
# =============================================================================
@pytest.fixture(scope="session")
def load_script() -> Callable[[str], ModuleType]:
    """Return a loader that imports a figure CLI script by filename stem."""
    def _load(stem: str) -> ModuleType:
        path = FIGURES_DIR / f"{stem}.py"
        if not path.exists():
            raise FileNotFoundError(f"No such script: {path}")
        spec = importlib.util.spec_from_file_location(f"_script_{stem}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    return _load
