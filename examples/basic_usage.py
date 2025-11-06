#!/usr/bin/env python3
"""
Basic usage examples for DNA Breathing Dynamics framework.

This script demonstrates fundamental operations including:
- DNA sequence validation
- Parameter usage
- High-precision calculations
- Statistical validation
"""

import numpy as np
from mpmath import mp

from src.core.params import (
    DEFAULT_MPMATH_DPS,
    KAPPA_GEO_DEFAULT,
    TEMPERATURE_DEFAULT,
    validate_dna_sequence,
    validate_temperature,
)


def example_sequence_validation() -> None:
    """Demonstrate DNA sequence validation."""
    print("\n" + "=" * 60)
    print("Example 1: DNA Sequence Validation")
    print("=" * 60)

    # Valid DNA sequence
    sequence = "ATCGATCGATCG"
    print(f"\nOriginal sequence: {sequence}")

    validated = validate_dna_sequence(sequence)
    print(f"Validated sequence: {validated}")

    # Lowercase conversion
    lowercase_seq = "atcg"
    validated_lower = validate_dna_sequence(lowercase_seq)
    print(f"\nLowercase '{lowercase_seq}' -> '{validated_lower}'")

    # RNA sequence
    rna_seq = "AUCGAUCG"
    validated_rna = validate_dna_sequence(rna_seq, allow_rna=True)
    print(f"RNA sequence '{rna_seq}' -> '{validated_rna}'")

    # Invalid sequence (will raise error)
    try:
        invalid_seq = "ATCGXYZ"
        validate_dna_sequence(invalid_seq)
    except ValueError as e:
        print(f"\nInvalid sequence 'ATCGXYZ' raises: {e}")


def example_parameter_usage() -> None:
    """Demonstrate parameter usage."""
    print("\n" + "=" * 60)
    print("Example 2: Parameter Usage")
    print("=" * 60)

    # Default parameters
    temp_celsius = TEMPERATURE_DEFAULT - 273.15
    print(f"\nDefault temperature: {TEMPERATURE_DEFAULT}K ({temp_celsius}°C)")
    print(f"Default kappa_geo: {KAPPA_GEO_DEFAULT}")
    print(f"Default mpmath precision: {DEFAULT_MPMATH_DPS} decimal places")

    # Validate custom temperature
    custom_temp = 298.15  # Room temperature
    validated_temp = validate_temperature(custom_temp)
    validated_celsius = validated_temp - 273.15
    print(f"\nValidated temperature: {validated_temp}K ({validated_celsius}°C)")


def example_high_precision() -> None:
    """Demonstrate high-precision calculations."""
    print("\n" + "=" * 60)
    print("Example 3: High-Precision Arithmetic")
    print("=" * 60)

    # Set precision
    mp.dps = DEFAULT_MPMATH_DPS
    print(f"\nPrecision: {mp.dps} decimal places")

    # High-precision calculation
    from mpmath import exp, pi, sqrt

    # Calculate e^π
    result = exp(pi)
    print(f"\ne^π = {result}")

    # Golden ratio
    phi = (1 + sqrt(5)) / 2
    print(f"φ (golden ratio) = {phi}")

    # Verify precision
    print(f"\nφ² - φ - 1 = {phi**2 - phi - 1}")  # Should be ~0


def example_statistical_analysis() -> None:
    """Demonstrate statistical analysis."""
    print("\n" + "=" * 60)
    print("Example 4: Statistical Analysis")
    print("=" * 60)

    # Generate sample data
    np.random.seed(42)
    gc_content_data = np.random.beta(a=2, b=2, size=100)

    print(f"\nGC content analysis (n={len(gc_content_data)})")
    print(f"Mean: {np.mean(gc_content_data):.4f}")
    print(f"Std: {np.std(gc_content_data):.4f}")
    print(f"Range: [{np.min(gc_content_data):.4f}, {np.max(gc_content_data):.4f}]")

    # Calculate confidence interval (manual)
    from scipy import stats

    ci = stats.t.interval(
        0.95,
        len(gc_content_data) - 1,
        loc=np.mean(gc_content_data),
        scale=stats.sem(gc_content_data),
    )
    print(f"95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]")


def main() -> None:
    """Run all examples."""
    print("\n" + "=" * 60)
    print("DNA Breathing Dynamics Framework - Basic Usage Examples")
    print("=" * 60)

    example_sequence_validation()
    example_parameter_usage()
    example_high_precision()
    example_statistical_analysis()

    print("\n" + "=" * 60)
    print("Examples completed successfully!")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
