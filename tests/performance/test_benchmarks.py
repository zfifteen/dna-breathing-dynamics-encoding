"""
Performance benchmark tests.

Tests performance characteristics and regression monitoring.
"""

import random

import pytest

pytest.importorskip(
    "pytest_benchmark", reason="pytest-benchmark plugin is optional for performance tests"
)

from src.core.params import validate_dna_sequence


@pytest.fixture(scope="module")
def long_sequence() -> str:
    """Pre-generated 10k bp DNA sequence for benchmarking."""
    rng = random.Random(0)
    return "".join(rng.choice("ACGT") for _ in range(10_000))


@pytest.fixture(scope="module")
def very_long_sequence() -> str:
    """Pre-generated 250k bp DNA sequence for benchmarking."""
    rng = random.Random(1)
    return "".join(rng.choice("ACGT") for _ in range(250_000))


@pytest.mark.performance
class TestSequenceValidationPerformance:
    """Benchmark sequence validation performance."""

    @pytest.mark.benchmark(group="validation")
    def test_short_sequence_validation(
        self, benchmark, sample_dna_sequence: str
    ) -> None:
        """Benchmark validation of short sequences."""
        result = benchmark(validate_dna_sequence, sample_dna_sequence)
        assert result is not None

    @pytest.mark.benchmark(group="validation")
    def test_long_sequence_validation(self, benchmark, long_sequence: str) -> None:
        """Benchmark validation of long sequences (10k bp)."""
        result = benchmark(validate_dna_sequence, long_sequence)
        assert result is not None

    @pytest.mark.benchmark(group="validation")
    @pytest.mark.slow
    def test_very_long_sequence_validation(
        self, benchmark, very_long_sequence: str
    ) -> None:
        """Benchmark validation of very long sequences (250k bp)."""
        result = benchmark(validate_dna_sequence, very_long_sequence)
        assert result is not None


@pytest.mark.performance
@pytest.mark.apple_silicon
class TestAppleSiliconOptimization:
    """Tests for Apple Silicon-specific optimizations."""

    def test_simd_compatibility(self) -> None:
        """Test SIMD operation compatibility (placeholder)."""
        # This would test actual SIMD operations when implemented
        assert True
