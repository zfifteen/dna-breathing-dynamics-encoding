"""
Scientific validation proof pack for DNA breathing dynamics.

This package provides rigorous statistical validation, empirical testing,
and reproducibility verification for research claims.
"""

from proof_pack.validation_framework import (
    ValidationResult,
    ValidationSuite,
    bootstrap_confidence_interval,
    validate_convergence,
    validate_numerical_accuracy,
    validate_statistical_hypothesis,
)

__all__ = [
    "ValidationResult",
    "ValidationSuite",
    "bootstrap_confidence_interval",
    "validate_convergence",
    "validate_numerical_accuracy",
    "validate_statistical_hypothesis",
]
