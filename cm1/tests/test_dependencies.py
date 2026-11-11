"""Ensure all dependencies declared in pyproject.toml are importable."""

import importlib
import importlib.metadata
from pathlib import Path
import tomllib
import pytest


def _pyproject_deps() -> list[str]:
    pyproject = Path(__file__).parent.parent.parent / "pyproject.toml"
    with pyproject.open("rb") as f:
        data = tomllib.load(f)
    return data["project"]["dependencies"]


def _canonical(dep: str) -> str:
    """Strip version specifiers and extras, return importable module name."""
    import re

    # grab just the package name before any version specifier
    name = re.split(r"[><=!;[]", dep)[0].strip()
    # pip name → import name (hyphens to underscores, known special cases)
    _IMPORT_NAMES = {
        "ipython": "IPython",
        "ipykernel": "ipykernel",
        "scikit-learn": "sklearn",
        "h5netcdf": "h5netcdf",
        "f90nml": "f90nml",
        "nwsspc-sharplib": "nwsspc",
    }
    return _IMPORT_NAMES.get(name, name.replace("-", "_"))


deps = _pyproject_deps()


@pytest.mark.parametrize("dep", deps)
def test_dependency_importable(dep):
    module = _canonical(dep)
    try:
        importlib.import_module(module)
    except ImportError:
        pytest.fail(
            f"'{dep}' declared in pyproject.toml but '{module}' is not importable"
        )
