"""
Core framework components for DNA breathing dynamics analysis.

This module provides fundamental utilities, parameter management,
and common functionality used across the framework.
"""

from src.core.params import (
    BREATHING_FREQUENCY_DEFAULT,
    DEFAULT_MPFR_PRECISION,
    DEFAULT_MPMATH_DPS,
    KAPPA_GEO_DEFAULT,
    KAPPA_STAR_DEFAULT,
    TEMPERATURE_DEFAULT,
    validate_dna_sequence,
    validate_kappa_geo,
    validate_temperature,
)

__all__ = [
    "DEFAULT_MPMATH_DPS",
    "DEFAULT_MPFR_PRECISION",
    "BREATHING_FREQUENCY_DEFAULT",
    "TEMPERATURE_DEFAULT",
    "KAPPA_GEO_DEFAULT",
    "KAPPA_STAR_DEFAULT",
    "validate_dna_sequence",
    "validate_kappa_geo",
    "validate_temperature",
]
