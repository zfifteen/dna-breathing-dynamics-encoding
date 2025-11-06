"""
Scientific validation framework for DNA breathing dynamics analysis.

This module provides statistical hypothesis testing, bootstrap confidence
intervals, and empirical validation following rigorous scientific standards.
"""

from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats


class ValidationResult:
    """Container for validation test results."""

    def __init__(
        self,
        test_name: str,
        passed: bool,
        p_value: Optional[float] = None,
        statistic: Optional[float] = None,
        confidence_interval: Optional[Tuple[float, float]] = None,
        effect_size: Optional[float] = None,
        metadata: Optional[Dict] = None,
    ):
        self.test_name = test_name
        self.passed = passed
        self.p_value = p_value
        self.statistic = statistic
        self.confidence_interval = confidence_interval
        self.effect_size = effect_size
        self.metadata = metadata or {}

    def __repr__(self) -> str:
        status = "PASS" if self.passed else "FAIL"
        return f"ValidationResult({self.test_name}: {status})"

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [f"\n{'='*60}", f"Test: {self.test_name}", f"Status: {'PASS' if self.passed else 'FAIL'}"]

        if self.p_value is not None:
            lines.append(f"p-value: {self.p_value:.6f}")

        if self.statistic is not None:
            lines.append(f"Test statistic: {self.statistic:.6f}")

        if self.confidence_interval is not None:
            ci_low, ci_high = self.confidence_interval
            lines.append(f"95% CI: [{ci_low:.6f}, {ci_high:.6f}]")

        if self.effect_size is not None:
            lines.append(f"Effect size: {self.effect_size:.6f}")

        lines.append("="*60)
        return "\n".join(lines)


def bootstrap_confidence_interval(
    data: np.ndarray,
    statistic_func: Callable[[np.ndarray], float],
    n_resamples: int = 1000,
    confidence_level: float = 0.95,
    random_seed: Optional[int] = None,
) -> Tuple[float, float]:
    """
    Calculate bootstrap confidence interval for a statistic.

    Args:
        data: Input data array
        statistic_func: Function to calculate statistic from data
        n_resamples: Number of bootstrap resamples (default: 1000)
        confidence_level: Confidence level (default: 0.95)
        random_seed: Random seed for reproducibility

    Returns:
        Tuple of (lower_bound, upper_bound) for confidence interval
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    n = len(data)
    bootstrap_statistics = np.zeros(n_resamples)

    for i in range(n_resamples):
        # Resample with replacement
        resample = np.random.choice(data, size=n, replace=True)
        bootstrap_statistics[i] = statistic_func(resample)

    # Calculate percentiles
    alpha = 1 - confidence_level
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100

    lower_bound = np.percentile(bootstrap_statistics, lower_percentile)
    upper_bound = np.percentile(bootstrap_statistics, upper_percentile)

    return (lower_bound, upper_bound)


def validate_statistical_hypothesis(
    observed_data: np.ndarray,
    expected_data: Optional[np.ndarray] = None,
    test_type: str = "t-test",
    alpha: float = 0.05,
    alternative: str = "two-sided",
) -> ValidationResult:
    """
    Perform statistical hypothesis testing.

    Args:
        observed_data: Observed data array
        expected_data: Expected/reference data (for paired/two-sample tests)
        test_type: Type of test ('t-test', 'wilcoxon', 'ks', 'chi2')
        alpha: Significance level (default: 0.05)
        alternative: Alternative hypothesis ('two-sided', 'less', 'greater')

    Returns:
        ValidationResult with test outcomes
    """
    if test_type == "t-test":
        if expected_data is None:
            # One-sample t-test against zero
            statistic, p_value = stats.ttest_1samp(observed_data, 0, alternative=alternative)
        else:
            # Two-sample t-test
            statistic, p_value = stats.ttest_ind(
                observed_data, expected_data, alternative=alternative
            )
    elif test_type == "wilcoxon":
        # Non-parametric Wilcoxon test
        if expected_data is None:
            statistic, p_value = stats.wilcoxon(observed_data, alternative=alternative)
        else:
            statistic, p_value = stats.ranksums(
                observed_data, expected_data, alternative=alternative
            )
    elif test_type == "ks":
        # Kolmogorov-Smirnov test
        if expected_data is None:
            # KS test against normal distribution
            statistic, p_value = stats.kstest(observed_data, "norm")
        else:
            statistic, p_value = stats.ks_2samp(observed_data, expected_data)
    elif test_type == "chi2":
        # Chi-squared goodness of fit
        statistic, p_value = stats.chisquare(observed_data, expected_data)
    else:
        raise ValueError(f"Unknown test type: {test_type}")

    passed = p_value < alpha

    # Calculate effect size (Cohen's d for t-tests)
    effect_size = None
    if test_type == "t-test" and expected_data is not None:
        pooled_std = np.sqrt(
            (np.std(observed_data, ddof=1) ** 2 + np.std(expected_data, ddof=1) ** 2) / 2
        )
        if pooled_std > 0:
            effect_size = (np.mean(observed_data) - np.mean(expected_data)) / pooled_std

    return ValidationResult(
        test_name=f"Statistical Test ({test_type})",
        passed=passed,
        p_value=p_value,
        statistic=statistic,
        effect_size=effect_size,
        metadata={"test_type": test_type, "alpha": alpha, "alternative": alternative},
    )


def validate_numerical_accuracy(
    computed_values: np.ndarray,
    reference_values: np.ndarray,
    rtol: float = 1e-9,
    atol: float = 1e-12,
) -> ValidationResult:
    """
    Validate numerical accuracy against reference values.

    Args:
        computed_values: Computed values to validate
        reference_values: Reference/ground truth values
        rtol: Relative tolerance
        atol: Absolute tolerance

    Returns:
        ValidationResult with accuracy metrics
    """
    passed = np.allclose(computed_values, reference_values, rtol=rtol, atol=atol)

    max_error = np.max(np.abs(computed_values - reference_values))
    mean_error = np.mean(np.abs(computed_values - reference_values))
    relative_error = np.max(
        np.abs((computed_values - reference_values) / (reference_values + 1e-15))
    )

    return ValidationResult(
        test_name="Numerical Accuracy",
        passed=passed,
        metadata={
            "max_absolute_error": max_error,
            "mean_absolute_error": mean_error,
            "max_relative_error": relative_error,
            "rtol": rtol,
            "atol": atol,
        },
    )


def validate_convergence(
    sequence: np.ndarray, tolerance: float = 1e-10, window: int = 10
) -> ValidationResult:
    """
    Validate convergence of iterative algorithm.

    Args:
        sequence: Sequence of values from iterative algorithm
        tolerance: Convergence tolerance
        window: Window size for checking convergence

    Returns:
        ValidationResult with convergence status
    """
    if len(sequence) < window:
        return ValidationResult(
            test_name="Convergence",
            passed=False,
            metadata={"reason": "Insufficient iterations"},
        )

    # Check last window for convergence
    recent_values = sequence[-window:]
    max_change = np.max(np.abs(np.diff(recent_values)))

    passed = max_change < tolerance

    return ValidationResult(
        test_name="Convergence",
        passed=passed,
        metadata={
            "max_change": max_change,
            "tolerance": tolerance,
            "iterations": len(sequence),
        },
    )


class ValidationSuite:
    """Orchestrate multiple validation tests."""

    def __init__(self, name: str):
        self.name = name
        self.results: List[ValidationResult] = []

    def add_result(self, result: ValidationResult) -> None:
        """Add validation result to suite."""
        self.results.append(result)

    def run_all(self) -> bool:
        """Check if all validations passed."""
        return all(result.passed for result in self.results)

    def summary(self) -> str:
        """Generate summary report."""
        total = len(self.results)
        passed = sum(1 for r in self.results if r.passed)
        failed = total - passed

        lines = [
            f"\n{'='*60}",
            f"Validation Suite: {self.name}",
            f"Total tests: {total}",
            f"Passed: {passed}",
            f"Failed: {failed}",
            f"Success rate: {100*passed/total:.1f}%",
            "="*60,
        ]

        for result in self.results:
            status = "✓" if result.passed else "✗"
            lines.append(f"{status} {result.test_name}")

        lines.append("="*60)
        return "\n".join(lines)

    def detailed_report(self) -> str:
        """Generate detailed report with all test results."""
        lines = [self.summary()]

        for result in self.results:
            if not result.passed:  # Show details for failed tests
                lines.append(result.summary())

        return "\n".join(lines)
