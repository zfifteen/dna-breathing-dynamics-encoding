"""
Centralized parameter management for DNA breathing dynamics framework.

This module provides standardized parameter definitions, validation functions,
and configuration management following scientific computing best practices.
"""

from typing import Optional, Union
import re

# =============================================================================
# Precision Settings
# =============================================================================

# High-precision arithmetic settings (mpmath)
DEFAULT_MPMATH_DPS = 50  # Decimal places for mpmath calculations

# MPFR precision for C extensions (bits)
DEFAULT_MPFR_PRECISION = 256  # 256-bit precision for ultra-high scale

# =============================================================================
# DNA Bioinformatics Parameters
# =============================================================================

# Valid nucleotide characters
DNA_NUCLEOTIDES = frozenset("ACGTN")  # N = any nucleotide
RNA_NUCLEOTIDES = frozenset("ACGUN")

# Default biological parameters
TEMPERATURE_DEFAULT = 310.15  # K (37°C, physiological temperature)
TEMPERATURE_MIN = 273.15  # K (0°C)
TEMPERATURE_MAX = 373.15  # K (100°C)

# DNA breathing dynamics defaults
BREATHING_FREQUENCY_DEFAULT = 1e6  # Hz (approximate base pair opening frequency)
BREATHING_AMPLITUDE_DEFAULT = 0.1  # nm (typical breathing amplitude)

# =============================================================================
# Mathematical Framework Parameters
# =============================================================================

# Geodesic mapping exponent (from Z Framework)
KAPPA_GEO_DEFAULT = 0.3
KAPPA_GEO_MIN = 0.0
KAPPA_GEO_MAX = 1.0

# Z_5D calibration parameter (from unified-framework)
KAPPA_STAR_DEFAULT = 0.04449
KAPPA_STAR_MIN = 0.0
KAPPA_STAR_MAX = 0.1

# Numerical stability thresholds
EPSILON_DEFAULT = 1e-15  # Machine epsilon threshold
CONVERGENCE_TOLERANCE = 1e-12  # Iterative solver tolerance

# =============================================================================
# Apple Silicon Optimization Parameters
# =============================================================================

# AMX tile dimensions for matrix operations
AMX_TILE_ROWS = 16
AMX_TILE_COLS = 64

# Memory alignment for SIMD operations (bytes)
MEMORY_ALIGNMENT = 64

# =============================================================================
# Validation Functions
# =============================================================================


def validate_dna_sequence(
    sequence: str,
    allow_rna: bool = False,
    allow_lowercase: bool = True
) -> str:
    """
    Validate DNA/RNA sequence composition.

    Args:
        sequence: Nucleotide sequence string
        allow_rna: If True, allow U (uracil) for RNA sequences
        allow_lowercase: If True, automatically convert to uppercase

    Returns:
        Validated and normalized sequence (uppercase)

    Raises:
        ValueError: If sequence contains invalid characters
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")

    # Convert to uppercase if allowed
    if allow_lowercase:
        sequence = sequence.upper()

    # Determine valid nucleotide set
    valid_nucleotides = RNA_NUCLEOTIDES if allow_rna else DNA_NUCLEOTIDES

    # Check for invalid characters
    invalid_chars = set(sequence) - valid_nucleotides
    if invalid_chars:
        raise ValueError(
            f"Invalid nucleotide characters: {sorted(invalid_chars)}. "
            f"Valid characters: {sorted(valid_nucleotides)}"
        )

    return sequence


def validate_temperature(temp: float) -> float:
    """
    Validate temperature parameter.

    Args:
        temp: Temperature in Kelvin

    Returns:
        Validated temperature

    Raises:
        ValueError: If temperature is outside valid range
    """
    if not TEMPERATURE_MIN <= temp <= TEMPERATURE_MAX:
        raise ValueError(
            f"Temperature {temp}K outside valid range "
            f"[{TEMPERATURE_MIN}, {TEMPERATURE_MAX}]K"
        )
    return temp


def validate_kappa_geo(kappa_geo: float) -> float:
    """
    Validate geodesic exponent parameter.

    Args:
        kappa_geo: Geodesic mapping exponent

    Returns:
        Validated kappa_geo

    Raises:
        ValueError: If kappa_geo is outside valid range
    """
    if not KAPPA_GEO_MIN <= kappa_geo <= KAPPA_GEO_MAX:
        raise ValueError(
            f"kappa_geo {kappa_geo} outside valid range "
            f"[{KAPPA_GEO_MIN}, {KAPPA_GEO_MAX}]"
        )
    return kappa_geo


def validate_kappa_star(kappa_star: float) -> float:
    """
    Validate Z_5D calibration parameter.

    Args:
        kappa_star: Z_5D calibration parameter

    Returns:
        Validated kappa_star

    Raises:
        ValueError: If kappa_star is outside valid range
    """
    if not KAPPA_STAR_MIN <= kappa_star <= KAPPA_STAR_MAX:
        raise ValueError(
            f"kappa_star {kappa_star} outside valid range "
            f"[{KAPPA_STAR_MIN}, {KAPPA_STAR_MAX}]"
        )
    return kappa_star


def validate_positive(value: float, name: str = "value") -> float:
    """
    Validate that a numeric parameter is positive.

    Args:
        value: Numeric value to validate
        name: Parameter name for error messages

    Returns:
        Validated value

    Raises:
        ValueError: If value is not positive
    """
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value}")
    return value


def validate_probability(prob: float, name: str = "probability") -> float:
    """
    Validate probability parameter is in [0, 1].

    Args:
        prob: Probability value
        name: Parameter name for error messages

    Returns:
        Validated probability

    Raises:
        ValueError: If probability is outside [0, 1]
    """
    if not 0.0 <= prob <= 1.0:
        raise ValueError(f"{name} must be in [0, 1], got {prob}")
    return prob
