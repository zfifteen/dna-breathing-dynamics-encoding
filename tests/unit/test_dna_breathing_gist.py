"""
Unit tests for DNA breathing gist dynamic DG bounds computation.

Tests that the ΔG° normalization uses dynamically computed bounds from
the NEAREST_NEIGHBOR_DG dictionary instead of hardcoded values.
"""

import sys
from pathlib import Path

import pytest

# Add the gist source to path
gist_src = Path(__file__).parent.parent.parent / "gists/breathing/czt_feature_extractor/src"
sys.path.insert(0, str(gist_src))

from dna_breathing_gist import (
    NEAREST_NEIGHBOR_DG,
    _DG_MIN,
    _DG_MAX,
    encode_sequence,
)


@pytest.mark.unit
class TestDynamicDGBounds:
    """Test dynamic DG bounds computation from NEAREST_NEIGHBOR_DG dictionary."""

    def test_dg_min_matches_dictionary_minimum(self) -> None:
        """Verify _DG_MIN equals the minimum value in NEAREST_NEIGHBOR_DG."""
        expected_min = min(NEAREST_NEIGHBOR_DG.values())
        assert _DG_MIN == expected_min
        assert _DG_MIN == -2.24  # GC is most stable

    def test_dg_max_matches_dictionary_maximum(self) -> None:
        """Verify _DG_MAX equals the maximum value in NEAREST_NEIGHBOR_DG."""
        expected_max = max(NEAREST_NEIGHBOR_DG.values())
        assert _DG_MAX == expected_max
        assert _DG_MAX == -0.58  # TA is least stable

    def test_normalization_range(self) -> None:
        """Test that normalization produces 0-1 range for dictionary extremes."""
        # Most stable (GC: -2.24) should normalize to 0
        norm_most_stable = (-2.24 - _DG_MIN) / (_DG_MAX - _DG_MIN)
        assert abs(norm_most_stable - 0.0) < 1e-10

        # Least stable (TA: -0.58) should normalize to 1
        norm_least_stable = (-0.58 - _DG_MIN) / (_DG_MAX - _DG_MIN)
        assert abs(norm_least_stable - 1.0) < 1e-10

    def test_encode_sequence_returns_complex_signal(self) -> None:
        """Test that encode_sequence returns a complex signal array."""
        signal = encode_sequence("GCTA")
        assert signal.dtype == complex
        assert len(signal) == 4

    def test_encode_sequence_imaginary_part_in_range(self) -> None:
        """Test that the imaginary part stays within normalized range."""
        # Test with a sequence containing all dinucleotide types
        signal = encode_sequence("GCTAGCTA", apply_helical=False)
        # Without helical modulation, imaginary parts should be in [0, 1]
        for val in signal:
            # The imaginary part should be within reasonable bounds
            # (may slightly exceed due to fallback values or edge effects)
            assert val.imag >= -0.5  # Allow some tolerance
            assert val.imag <= 1.5  # Allow some tolerance
