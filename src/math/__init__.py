"""
Mathematical modeling and analysis components.

This module provides high-precision numerical computation, statistical
analysis, and mathematical framework implementations for DNA dynamics.
"""

from typing import List

from src.math.kappa_weighted import (
    kappa,
    theta_prime,
    dna_to_complex,
    apply_phase_shift,
    compute_spectral_entropy,
    disruption_score,
    bootstrap_ci,
    compute_correlation_metrics,
    compute_auc_metrics,
    generate_validation_report,
)

__all__: List[str] = [
    "kappa",
    "theta_prime",
    "dna_to_complex",
    "apply_phase_shift",
    "compute_spectral_entropy",
    "disruption_score",
    "bootstrap_ci",
    "compute_correlation_metrics",
    "compute_auc_metrics",
    "generate_validation_report",
]
