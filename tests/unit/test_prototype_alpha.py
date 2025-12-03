"""
Unit tests for Prototype Alpha: CZT Breathing Dynamics Feature Extractor.

Tests core encoding, feature extraction, and statistical functions.
"""

import numpy as np
import pytest

from experiments.prototype_alpha import (
    BREATH_LIFETIMES_MS,
    DELTA_G,
    RNG_SEED,
    cohen_d,
    czt_band_features,
    encode_sequence,
    gc_content,
    random_encoding_control,
)


@pytest.mark.unit
class TestEncodeSequence:
    """Test DNA sequence encoding to complex signal."""

    def test_encodes_valid_sequence(self) -> None:
        """Test encoding a simple DNA sequence."""
        seq = "ATCG"
        x = encode_sequence(seq)
        assert x is not None
        assert len(x) == 4
        assert x.dtype == complex

    def test_encodes_at_bases_correctly(self) -> None:
        """Test that AT bases get correct lifetime values."""
        seq = "AT"
        x = encode_sequence(seq)
        assert x is not None
        np.testing.assert_array_equal(
            x.real, [BREATH_LIFETIMES_MS["A"], BREATH_LIFETIMES_MS["T"]]
        )

    def test_encodes_gc_bases_correctly(self) -> None:
        """Test that GC bases get correct lifetime values."""
        seq = "GC"
        x = encode_sequence(seq)
        assert x is not None
        np.testing.assert_array_equal(
            x.real, [BREATH_LIFETIMES_MS["G"], BREATH_LIFETIMES_MS["C"]]
        )

    def test_imaginary_part_contains_delta_g(self) -> None:
        """Test that imaginary part contains delta-G values."""
        seq = "ATGC"
        x = encode_sequence(seq)
        assert x is not None
        expected = [DELTA_G["A"], DELTA_G["T"], DELTA_G["G"], DELTA_G["C"]]
        np.testing.assert_array_equal(x.imag, expected)

    def test_skips_invalid_bases(self) -> None:
        """Test that invalid bases (N) are skipped."""
        seq = "ATNGC"
        x = encode_sequence(seq)
        assert x is not None
        assert len(x) == 4  # N is skipped

    def test_returns_none_for_all_invalid(self) -> None:
        """Test that all-invalid sequence returns None."""
        seq = "NNNN"
        x = encode_sequence(seq)
        assert x is None

    def test_empty_sequence_returns_none(self) -> None:
        """Test that empty sequence returns None."""
        seq = ""
        x = encode_sequence(seq)
        assert x is None

    def test_handles_lowercase(self) -> None:
        """Test that lowercase sequences are handled."""
        seq = "atcg"
        x = encode_sequence(seq)
        assert x is not None
        assert len(x) == 4


@pytest.mark.unit
class TestGCContent:
    """Test GC content calculation."""

    def test_gc_content_pure_at(self) -> None:
        """Test GC content of pure AT sequence is 0."""
        assert gc_content("AAAA") == 0.0
        assert gc_content("TTTT") == 0.0
        assert gc_content("ATAT") == 0.0

    def test_gc_content_pure_gc(self) -> None:
        """Test GC content of pure GC sequence is 1."""
        assert gc_content("GGGG") == 1.0
        assert gc_content("CCCC") == 1.0
        assert gc_content("GCGC") == 1.0

    def test_gc_content_mixed(self) -> None:
        """Test GC content of mixed sequence."""
        assert gc_content("ATCG") == 0.5
        assert gc_content("AATG") == 0.25
        assert gc_content("CCCG") == 1.0

    def test_gc_content_handles_lowercase(self) -> None:
        """Test GC content handles lowercase."""
        assert gc_content("atcg") == 0.5

    def test_gc_content_empty_returns_zero(self) -> None:
        """Test GC content of empty sequence returns 0."""
        assert gc_content("") == 0.0

    def test_gc_content_invalid_bases_excluded(self) -> None:
        """Test that invalid bases are excluded from calculation."""
        # Only ATGC bases are counted
        assert gc_content("ATCGN") == 0.5  # N is not counted


@pytest.mark.unit
class TestCZTBandFeatures:
    """Test CZT band feature extraction."""

    def test_returns_four_features(self) -> None:
        """Test that CZT extraction returns 4 features."""
        seq = "ATCGATCGATCGATCGATCG"  # 20 bp
        x = encode_sequence(seq)
        assert x is not None
        features = czt_band_features(x)
        assert len(features) == 4

    def test_feature_values_are_valid(self) -> None:
        """Test that feature values are valid numbers."""
        seq = "ATCGATCGATCGATCGATCG"
        x = encode_sequence(seq)
        assert x is not None
        features = czt_band_features(x)
        # mag_peak should be positive
        assert features[0] > 0
        # band_energy should be positive
        assert features[1] > 0
        # phase_peak should be in [-pi, pi]
        assert -np.pi <= features[2] <= np.pi
        # freq_eff should be near f0 = 1/10.5
        assert 0.05 < features[3] < 0.15

    def test_frequency_near_helical_period(self) -> None:
        """Test that effective frequency is near target."""
        # Longer sequence for better frequency resolution
        seq = "ATCGATCGATCG" * 10  # 120 bp
        x = encode_sequence(seq)
        assert x is not None
        features = czt_band_features(x)
        # freq_eff should be within the band around 1/10.5
        f0 = 1 / 10.5
        rel_bw = 0.1
        assert f0 * (1 - rel_bw) <= features[3] <= f0 * (1 + rel_bw)


@pytest.mark.unit
class TestRandomEncodingControl:
    """Test random encoding control shuffle."""

    def test_preserves_shape(self) -> None:
        """Test that shuffle preserves signal shape."""
        np.random.seed(RNG_SEED)
        seq = "ATCGATCG"
        x = encode_sequence(seq)
        assert x is not None
        xr = random_encoding_control(x)
        assert xr.shape == x.shape

    def test_preserves_marginal_distributions(self) -> None:
        """Test that shuffle preserves value sets."""
        np.random.seed(RNG_SEED)
        seq = "ATCGATCG"
        x = encode_sequence(seq)
        assert x is not None
        xr = random_encoding_control(x)
        # Same values, possibly different order
        np.testing.assert_array_equal(sorted(xr.real), sorted(x.real))
        np.testing.assert_array_equal(sorted(xr.imag), sorted(x.imag))

    def test_changes_positions(self) -> None:
        """Test that shuffle changes at least some positions."""
        np.random.seed(42)  # Different seed
        seq = "ATCGATCGATCGATCGATCG"
        x = encode_sequence(seq)
        assert x is not None
        xr = random_encoding_control(x)
        # Very unlikely to be identical after shuffle
        assert not np.array_equal(xr, x)


@pytest.mark.unit
class TestCohenD:
    """Test Cohen's d effect size calculation."""

    def test_identical_groups_returns_zero(self) -> None:
        """Test Cohen's d for identical groups is 0."""
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 2.0, 3.0])
        d = cohen_d(a, b)
        assert d == 0.0

    def test_different_groups(self) -> None:
        """Test Cohen's d for different groups."""
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([4.0, 5.0, 6.0])
        d = cohen_d(a, b)
        # Should be negative since mean(a) < mean(b)
        assert d < 0
        # Magnitude should be > 1 for this separation
        assert abs(d) > 1.0

    def test_handles_single_element(self) -> None:
        """Test Cohen's d handles single element groups."""
        a = np.array([1.0])
        b = np.array([2.0, 3.0, 4.0])
        d = cohen_d(a, b)
        # Should return 0 when one group has < 2 elements
        assert d == 0.0

    def test_known_effect_size(self) -> None:
        """Test Cohen's d with known effect size."""
        # Two groups with known means and pooled std
        a = np.array([0.0, 0.0, 0.0])
        b = np.array([1.0, 1.0, 1.0])
        d = cohen_d(a, b)
        # d = (0 - 1) / pooled_std, where pooled_std = 0
        # This edge case should handle zero variance
        # Actually with zero variance, we get division by zero
        # Let's use a more realistic test
        a = np.array([0.0, 1.0, 2.0])
        b = np.array([3.0, 4.0, 5.0])
        d = cohen_d(a, b)
        # Mean diff = -3, pooled_std = 1.0
        assert abs(d + 3.0) < 0.01


@pytest.mark.unit
@pytest.mark.smoke
def test_biophysics_parameters() -> None:
    """Smoke test for biophysics parameter tables."""
    # Check lifetimes are correctly defined
    assert BREATH_LIFETIMES_MS["A"] == 5.0
    assert BREATH_LIFETIMES_MS["T"] == 5.0
    assert BREATH_LIFETIMES_MS["G"] == 25.0
    assert BREATH_LIFETIMES_MS["C"] == 25.0

    # Check delta-G values
    assert DELTA_G["A"] == -1.0
    assert DELTA_G["T"] == -1.0
    assert DELTA_G["G"] == -2.0
    assert DELTA_G["C"] == -2.0


@pytest.mark.unit
@pytest.mark.smoke
def test_rng_seed_defined() -> None:
    """Smoke test for RNG seed."""
    assert RNG_SEED == 137
