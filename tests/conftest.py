"""
Pytest configuration and shared fixtures for DNA breathing dynamics tests.

This module provides common test fixtures, utilities, and configuration
for the entire test suite.
"""

import random
from pathlib import Path
from typing import Generator

import numpy as np
import pytest
from mpmath import mp


def pytest_load_initial_conftests(early_config, parser, args):
    """
    Keep default addopts usable even when pytest-cov isn't installed.

    The project-level addopts include coverage flags; stripping them here
    prevents failures that force callers to bypass config with -c /dev/null,
    which in turn avoids marker and cache-dir warnings.
    """
    try:
        import pytest_cov  # type: ignore  # noqa: F401
    except ImportError:
        cov_flags = {"--cov=src", "--cov-report=term-missing", "--cov-report=html"}
        args[:] = [arg for arg in args if arg not in cov_flags]


def pytest_addoption(parser):
    """
    Define no-op coverage options so pytest doesn't error when pytest-cov
    is absent. If pytest-cov is installed it will override these options.
    """
    parser.addoption(
        "--cov",
        action="append",
        default=[],
        help="Dummy coverage option; real behavior comes from pytest-cov if installed.",
    )
    parser.addoption(
        "--cov-report",
        action="append",
        default=[],
        help="Dummy coverage option; real behavior comes from pytest-cov if installed.",
    )


from src.core.params import (
    BREATHING_FREQUENCY_DEFAULT,
    DEFAULT_MPMATH_DPS,
    TEMPERATURE_DEFAULT,
)


@pytest.fixture(scope="session", autouse=True)
def set_random_seeds() -> None:
    """Set random seeds for reproducibility across all tests."""
    random.seed(42)
    np.random.seed(42)


@pytest.fixture(scope="session", autouse=True)
def configure_mpmath() -> None:
    """Configure mpmath precision for all tests."""
    mp.dps = DEFAULT_MPMATH_DPS


def pytest_configure(config):
    """Register markers explicitly for environments that skip pyproject parsing."""
    markers = [
        "unit: Unit tests",
        "integration: Integration tests",
        "performance: Performance benchmark tests",
        "validation: Scientific validation tests",
        "smoke: Quick smoke tests (<5s)",
        "slow: Slow tests (>5s)",
        "apple_silicon: Tests requiring Apple Silicon hardware",
    ]
    for marker in markers:
        config.addinivalue_line("markers", marker)

    # If config rootdir is forced to /dev (e.g., via -c /dev/null), redirect
    # pytest's cache to a writable path inside the repository to avoid warnings.
    if str(config.rootpath) == "/dev":
        repo_root = Path(__file__).resolve().parents[1]
        cache_dir = repo_root / ".pytest_cache"
        cache_dir.mkdir(exist_ok=True)
        # Update both option and cache object so cacheprovider writes inside repo
        config.option.cache_dir = str(cache_dir)
        if hasattr(config, "cache"):
            config.cache._cachedir = cache_dir


@pytest.fixture
def sample_dna_sequence() -> str:
    """Provide a sample DNA sequence for testing."""
    return "ATCGATCGATCGATCG"


@pytest.fixture
def sample_rna_sequence() -> str:
    """Provide a sample RNA sequence for testing."""
    return "AUCGAUCGAUCGAUCG"


@pytest.fixture
def long_dna_sequence() -> str:
    """Provide a longer DNA sequence for integration tests."""
    bases = "ACGT"
    return "".join(random.choice(bases) for _ in range(1000))


@pytest.fixture
def gc_rich_sequence() -> str:
    """Provide a GC-rich DNA sequence for testing."""
    return "GCGCGCGCGCGCGCGC"


@pytest.fixture
def at_rich_sequence() -> str:
    """Provide an AT-rich DNA sequence for testing."""
    return "ATATATATATATAT"


@pytest.fixture
def default_temperature() -> float:
    """Provide default physiological temperature."""
    return TEMPERATURE_DEFAULT


@pytest.fixture
def temperature_range() -> Generator[float, None, None]:
    """Provide a range of temperatures for testing."""
    for temp in np.linspace(273.15, 373.15, 10):
        yield temp


@pytest.fixture
def default_breathing_frequency() -> float:
    """Provide default breathing frequency."""
    return BREATHING_FREQUENCY_DEFAULT


@pytest.fixture
def numpy_array_1d() -> np.ndarray:
    """Provide a 1D numpy array for testing."""
    return np.array([1.0, 2.0, 3.0, 4.0, 5.0])


@pytest.fixture
def numpy_array_2d() -> np.ndarray:
    """Provide a 2D numpy array for testing."""
    return np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
