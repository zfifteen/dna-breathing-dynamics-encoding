"""
Pytest configuration and shared fixtures for DNA breathing dynamics tests.

This module provides common test fixtures, utilities, and configuration
for the entire test suite.
"""

import random
from typing import Generator

import numpy as np
import pytest
from mpmath import mp

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


# Markers for categorizing tests
pytest_markers = {
    "unit": "Unit tests",
    "integration": "Integration tests",
    "performance": "Performance benchmark tests",
    "validation": "Scientific validation tests",
    "smoke": "Quick smoke tests (<5s)",
    "slow": "Slow tests (>5s)",
    "apple_silicon": "Tests requiring Apple Silicon hardware",
}
