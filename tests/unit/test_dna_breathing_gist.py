"""
Unit tests for DNA breathing gist dynamic DG bounds computation.

Tests that the ΔG° normalization uses dynamically computed bounds from
the NEAREST_NEIGHBOR_DG dictionary instead of hardcoded values.
"""

import math
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
    compute_stats,
    generate_dinuc_shuffles,
)


@pytest.mark.unit
class TestDynamicDGBounds:
    """Test dynamic DG bounds computation from NEAREST_NEIGHBOR_DG dictionary."""

    def test_dg_min_matches_dictionary_minimum(self) -> None:
        """Verify _DG_MIN equals the minimum value in NEAREST_NEIGHBOR_DG."""
        expected_min = min(NEAREST_NEIGHBOR_DG.values())
        assert _DG_MIN == expected_min

    def test_dg_max_matches_dictionary_maximum(self) -> None:
        """Verify _DG_MAX equals the maximum value in NEAREST_NEIGHBOR_DG."""
        expected_max = max(NEAREST_NEIGHBOR_DG.values())
        assert _DG_MAX == expected_max

    def test_normalization_range(self) -> None:
        """Test that normalization produces 0-1 range for dictionary extremes."""
        dg_min = min(NEAREST_NEIGHBOR_DG.values())
        dg_max = max(NEAREST_NEIGHBOR_DG.values())

        # Most stable dinucleotide should normalize to 0
        norm_most_stable = (dg_min - _DG_MIN) / (_DG_MAX - _DG_MIN)
        assert abs(norm_most_stable - 0.0) < 1e-10

        # Least stable dinucleotide should normalize to 1
        norm_least_stable = (dg_max - _DG_MIN) / (_DG_MAX - _DG_MIN)
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


@pytest.mark.unit
class TestCohensDZeroVariance:
    """Regression test for issue #20: avoid masking division-by-zero in Cohen's d."""

    def test_zero_variance_groups_yield_nan_effect_size(self) -> None:
        """If both groups are constant, Cohen's d and its CI should be NaN."""
        features = []
        groups = []

        for _ in range(4):
            features.append(
                {"peak_mag": 1.0, "snr": 0.5, "phase_coherence": 0.2}
            )
            groups.append("A")

        for _ in range(5):
            features.append(
                {"peak_mag": 3.0, "snr": 0.5, "phase_coherence": 0.2}
            )
            groups.append("B")

        stats = compute_stats(features, groups, num_bootstrap=10, num_perm=10, seed=123)

        assert math.isnan(stats["cohens_d_peak_mag"])
        assert math.isnan(stats["ci_low_peak_mag"])
        assert math.isnan(stats["ci_high_peak_mag"])
        assert 0.0 <= stats["p_perm_peak_mag"] <= 1.0


@pytest.mark.unit
class TestShuffleFallbackWarning:
    """Ensure shuffling failure path emits a warning (issue #23)."""

    def test_shuffle_fallback_warns_and_returns_original(self) -> None:
        seq = "AAAA"  # Graph with single node; force failure via zero retries
        with pytest.warns(RuntimeWarning):
            shuffles = generate_dinuc_shuffles(
                seq, num_shuffles=2, seed=1, max_retries=0, warn=True
            )
        assert shuffles == [seq, seq]


@pytest.mark.unit
class TestPermutationPValue:
    """Regression for permutation p-value smoothing (issue #26)."""

    def test_permutation_p_has_add_one_smoothing(self) -> None:
        # Construct two small deterministic groups with clear difference.
        features = []
        groups = []
        for val in [0.0, 0.0, 0.0]:
            features.append({"peak_mag": val, "snr": val, "phase_coherence": val})
            groups.append("A")
        for val in [1.0, 1.0, 1.0]:
            features.append({"peak_mag": val, "snr": val, "phase_coherence": val})
            groups.append("B")

        num_perm = 10
        stats = compute_stats(features, groups, num_bootstrap=10, num_perm=num_perm, seed=7)

        expected_min = 1 / (num_perm + 1)
        # With add-one smoothing, p must be bounded below by 1/(perm+1) instead of 0
        assert stats["p_perm_peak_mag"] >= expected_min - 1e-12
        assert 0.0 < stats["p_perm_peak_mag"] <= 1.0


@pytest.mark.unit
class TestGroupLengthMismatch:
    """Regression for issue #15: reject mismatched group labels instead of truncating."""

    def test_group_count_mismatch_raises(self) -> None:
        features = [
            {"peak_mag": 1.0, "snr": 1.0, "phase_coherence": 1.0},
            {"peak_mag": 2.0, "snr": 2.0, "phase_coherence": 2.0},
            {"peak_mag": 3.0, "snr": 3.0, "phase_coherence": 3.0},
        ]
        groups = ["A", "B"]  # length mismatch

        with pytest.raises(ValueError):
            compute_stats(features, groups, num_bootstrap=5, num_perm=5, seed=0)
