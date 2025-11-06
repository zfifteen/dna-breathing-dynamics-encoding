#!/usr/bin/env python3
"""
Quick validation demonstration (~2 min runtime).

Demonstrates scientific validation framework with example tests.
"""

import numpy as np

from proof_pack.validation_framework import (
    ValidationSuite,
    bootstrap_confidence_interval,
    validate_convergence,
    validate_numerical_accuracy,
    validate_statistical_hypothesis,
)


def main() -> None:
    """Run quick validation demonstration."""
    print("\nDNA Breathing Dynamics - Quick Validation Demo")
    print("=" * 60)

    suite = ValidationSuite("Quick Validation Demo")

    # Test 1: Statistical hypothesis testing
    print("\n[1/4] Statistical hypothesis testing...")
    observed = np.random.normal(loc=1.0, scale=0.5, size=100)
    expected = np.random.normal(loc=0.0, scale=0.5, size=100)

    result = validate_statistical_hypothesis(
        observed, expected, test_type="t-test", alpha=0.05
    )
    suite.add_result(result)
    print(f"    Result: {'PASS' if result.passed else 'FAIL'} (p={result.p_value:.6f})")

    # Test 2: Bootstrap confidence intervals
    print("\n[2/4] Bootstrap confidence interval...")
    data = np.random.exponential(scale=2.0, size=200)
    ci = bootstrap_confidence_interval(
        data, statistic_func=np.mean, n_resamples=1000, random_seed=42
    )
    print(f"    95% CI for mean: [{ci[0]:.4f}, {ci[1]:.4f}]")
    print(f"    Sample mean: {np.mean(data):.4f}")

    # Test 3: Numerical accuracy
    print("\n[3/4] Numerical accuracy validation...")
    computed = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    reference = np.array([1.0000001, 2.0000001, 3.0000001, 4.0000001, 5.0000001])

    result = validate_numerical_accuracy(computed, reference, rtol=1e-6, atol=1e-6)
    suite.add_result(result)
    print(f"    Result: {'PASS' if result.passed else 'FAIL'}")
    print(f"    Max error: {result.metadata['max_absolute_error']:.2e}")

    # Test 4: Convergence validation
    print("\n[4/4] Convergence validation...")
    # Simulate converging sequence
    iterations = np.array([1 / (i + 1) for i in range(100)])

    result = validate_convergence(iterations, tolerance=1e-2, window=10)
    suite.add_result(result)
    print(f"    Result: {'PASS' if result.passed else 'FAIL'}")
    print(f"    Final value: {iterations[-1]:.6f}")

    # Summary
    print("\n" + suite.summary())

    if suite.run_all():
        print("\n✓ All validations passed!")
        return

    print("\n✗ Some validations failed. See details above.")


if __name__ == "__main__":
    main()
