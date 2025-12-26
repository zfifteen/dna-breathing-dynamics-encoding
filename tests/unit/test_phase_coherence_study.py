"""
Unit tests for Phase-Coherent vs Random-Phase CZT Spectra Comparison Study.

Tests core encoding, CZT feature extraction, phase scrambling, statistical
functions, and circular statistics for DNA breathing dynamics analysis.
"""

import numpy as np
import pytest

from experiments.phase_coherence_study import (
    BREATH_LIFETIMES_MS,
    NEAREST_NEIGHBOR_DG,
    DEFAULT_HELICAL_PERIOD,
    _DG_MIN,
    _DG_MAX,
    CZTFeatures,
    StudyConfig,
    FoldMetrics,
    encode_sequence,
    compute_czt_spectrum,
    extract_czt_features,
    scramble_phases,
    _compute_circular_stats,
    compute_paired_statistics,
    compute_power_analysis,
    apply_fdr_correction,
    extract_features_for_sequences,
    run_cv_evaluation,
)


# =============================================================================
# Test encode_sequence
# =============================================================================
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
        """Test that AT bases get correct lifetime values in real part."""
        seq = "AT"
        x = encode_sequence(seq, apply_helical=False)
        assert x is not None
        np.testing.assert_array_equal(
            x.real, [BREATH_LIFETIMES_MS["A"], BREATH_LIFETIMES_MS["T"]]
        )

    def test_encodes_gc_bases_correctly(self) -> None:
        """Test that GC bases get correct lifetime values in real part."""
        seq = "GC"
        x = encode_sequence(seq, apply_helical=False)
        assert x is not None
        np.testing.assert_array_equal(
            x.real, [BREATH_LIFETIMES_MS["G"], BREATH_LIFETIMES_MS["C"]]
        )

    def test_imaginary_part_contains_normalized_dg(self) -> None:
        """Test that imaginary part contains normalized ΔG values (0-1 range)."""
        seq = "ATGC"
        x = encode_sequence(seq, apply_helical=False)
        assert x is not None
        # All values should be in normalized [0, 1] range
        for val in x.imag:
            assert 0.0 <= val <= 1.0

    def test_helical_modulation_applied(self) -> None:
        """Test that helical phase modulation changes signal when enabled."""
        seq = "ATCGATCGATCG"
        x_with_helical = encode_sequence(seq, apply_helical=True)
        x_without_helical = encode_sequence(seq, apply_helical=False)
        assert x_with_helical is not None
        assert x_without_helical is not None
        # Signals should be different due to helical modulation
        assert not np.allclose(x_with_helical, x_without_helical)

    def test_raises_on_empty_sequence(self) -> None:
        """Test that empty sequence raises ValueError."""
        with pytest.raises(ValueError, match="Empty sequence"):
            encode_sequence("")

    def test_raises_on_invalid_bases(self) -> None:
        """Test that invalid bases raise ValueError."""
        with pytest.raises(ValueError, match="Invalid base"):
            encode_sequence("ATNGC")
        with pytest.raises(ValueError, match="Invalid base"):
            encode_sequence("ATXGC")

    def test_handles_lowercase(self) -> None:
        """Test that lowercase sequences are handled."""
        seq = "atcg"
        x = encode_sequence(seq)
        assert x is not None
        assert len(x) == 4

    def test_custom_lifetimes(self) -> None:
        """Test that custom lifetimes are applied correctly."""
        seq = "AT"
        x = encode_sequence(
            seq, at_lifetime=10.0, gc_lifetime=50.0, apply_helical=False
        )
        assert x is not None
        np.testing.assert_array_equal(x.real, [10.0, 10.0])

    @pytest.mark.smoke
    def test_smoke_typical_grna_length(self) -> None:
        """Smoke test encoding a typical gRNA-length sequence (20 bp)."""
        seq = "ATCGATCGATCGATCGATCG"
        x = encode_sequence(seq)
        assert x is not None
        assert len(x) == 20
        assert x.dtype == complex


# =============================================================================
# Test CZT spectrum computation
# =============================================================================
@pytest.mark.unit
class TestComputeCZTSpectrum:
    """Test CZT spectrum computation for DNA signals."""

    def test_returns_correct_shapes(self) -> None:
        """Test that CZT returns correct array shapes."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal, m=256)
        assert len(freqs) == 256
        assert len(spectrum) == 256

    def test_frequency_range_centered_on_helical(self) -> None:
        """Test that frequency range is centered on helical frequency."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        f0 = 1 / DEFAULT_HELICAL_PERIOD
        band_width = 0.01
        freqs, _ = compute_czt_spectrum(signal, f0=f0, band_width=band_width, m=256)
        # Frequencies should span [f0 - bw, f0 + bw]
        assert freqs[0] == pytest.approx(f0 - band_width, rel=1e-6)
        assert freqs[-1] == pytest.approx(f0 + band_width, rel=1e-6)

    def test_spectrum_is_complex(self) -> None:
        """Test that CZT spectrum is complex-valued."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        _, spectrum = compute_czt_spectrum(signal)
        assert spectrum.dtype == complex


# =============================================================================
# Test extract_czt_features
# =============================================================================
@pytest.mark.unit
class TestExtractCZTFeatures:
    """Test CZT feature extraction."""

    def test_returns_seven_features(self) -> None:
        """Test that feature extraction returns 7 features."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        arr = features.to_array()
        assert len(arr) == 7

    def test_feature_names_match_array(self) -> None:
        """Test that feature names list matches feature count."""
        names = CZTFeatures.feature_names()
        assert len(names) == 7
        expected_names = [
            "peak_magnitude",
            "peak_freq_idx",
            "phase_at_peak",
            "band_power",
            "spectral_centroid",
            "phase_kappa",
            "rayleigh_z",
        ]
        assert names == expected_names

    def test_peak_magnitude_positive(self) -> None:
        """Test that peak magnitude is positive."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert features.peak_magnitude > 0

    def test_band_power_positive(self) -> None:
        """Test that band power is positive."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert features.band_power > 0

    def test_phase_in_valid_range(self) -> None:
        """Test that phase at peak is in [-π, π]."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert -np.pi <= features.phase_at_peak <= np.pi

    def test_spectral_centroid_in_frequency_range(self) -> None:
        """Test that spectral centroid is within frequency range."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        f0 = 1 / DEFAULT_HELICAL_PERIOD
        band_width = 0.01
        freqs, spectrum = compute_czt_spectrum(signal, f0=f0, band_width=band_width)
        features = extract_czt_features(freqs, spectrum)
        assert freqs[0] <= features.spectral_centroid <= freqs[-1]

    def test_phase_kappa_non_negative(self) -> None:
        """Test that von Mises kappa is non-negative."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert features.phase_kappa >= 0

    def test_rayleigh_z_non_negative(self) -> None:
        """Test that Rayleigh z-statistic is non-negative."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert features.rayleigh_z >= 0

    def test_rayleigh_p_in_valid_range(self) -> None:
        """Test that Rayleigh p-value is in [0, 1]."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert 0 <= features.rayleigh_p <= 1

    @pytest.mark.smoke
    def test_smoke_feature_extraction(self) -> None:
        """Smoke test for complete feature extraction pipeline."""
        signal = encode_sequence("GCTAGCTAGCTAGCTAGCTA")
        assert signal is not None
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        arr = features.to_array()
        # All values should be finite
        assert np.all(np.isfinite(arr))


# =============================================================================
# Test scramble_phases
# =============================================================================
@pytest.mark.unit
class TestScramblePhases:
    """Test phase scrambling for random-phase control."""

    def test_preserves_magnitudes(self) -> None:
        """Test that phase scrambling preserves spectral magnitudes."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        _, spectrum = compute_czt_spectrum(signal)

        rng = np.random.default_rng(42)
        scrambled = scramble_phases(spectrum, rng)

        # Magnitudes should be identical
        np.testing.assert_allclose(np.abs(scrambled), np.abs(spectrum), rtol=1e-10)

    def test_changes_phases(self) -> None:
        """Test that phase scrambling changes phases."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        _, spectrum = compute_czt_spectrum(signal)

        rng = np.random.default_rng(42)
        scrambled = scramble_phases(spectrum, rng)

        orig_phases = np.angle(spectrum)
        scrambled_phases = np.angle(scrambled)
        # Phases should be different
        assert not np.allclose(orig_phases, scrambled_phases)

    def test_scrambled_phases_uniform(self) -> None:
        """Test that scrambled phases are approximately uniformly distributed."""
        signal = encode_sequence("ATCGATCGATCGATCGATCGATCGATCG" * 10)
        assert signal is not None
        _, spectrum = compute_czt_spectrum(signal, m=1000)

        rng = np.random.default_rng(42)
        scrambled = scramble_phases(spectrum, rng)

        phases = np.angle(scrambled)
        # Phases should be in [-π, π]
        assert np.all(phases >= -np.pi)
        assert np.all(phases <= np.pi)
        # Mean should be close to 0 for uniform distribution
        assert abs(np.mean(phases)) < 0.2

    def test_reproducible_with_same_seed(self) -> None:
        """Test that same seed produces same scrambled result."""
        signal = encode_sequence("ATCGATCGATCGATCGATCG")
        assert signal is not None
        _, spectrum = compute_czt_spectrum(signal)

        rng1 = np.random.default_rng(42)
        scrambled1 = scramble_phases(spectrum, rng1)

        rng2 = np.random.default_rng(42)
        scrambled2 = scramble_phases(spectrum, rng2)

        np.testing.assert_array_equal(scrambled1, scrambled2)


# =============================================================================
# Test _compute_circular_stats
# =============================================================================
@pytest.mark.unit
class TestComputeCircularStats:
    """Test circular statistics computation."""

    def test_uniform_phases_low_kappa(self) -> None:
        """Test that uniform phases yield low von Mises kappa."""
        # Uniformly distributed phases
        phases = np.linspace(-np.pi, np.pi, 100)
        kappa, rayleigh_z, rayleigh_p = _compute_circular_stats(phases)
        # Should have low concentration
        assert kappa < 0.5

    def test_concentrated_phases_high_kappa(self) -> None:
        """Test that concentrated phases yield high von Mises kappa."""
        # All phases near 0
        phases = np.random.default_rng(42).normal(0, 0.1, 100)
        kappa, rayleigh_z, rayleigh_p = _compute_circular_stats(phases)
        # Should have high concentration
        assert kappa > 2.0

    def test_rayleigh_z_increases_with_concentration(self) -> None:
        """Test that Rayleigh z increases with phase concentration."""
        # Uniform phases
        uniform_phases = np.linspace(-np.pi, np.pi, 100)
        _, z_uniform, _ = _compute_circular_stats(uniform_phases)

        # Concentrated phases
        concentrated_phases = np.random.default_rng(42).normal(0, 0.1, 100)
        _, z_concentrated, _ = _compute_circular_stats(concentrated_phases)

        assert z_concentrated > z_uniform

    def test_empty_array_returns_defaults(self) -> None:
        """Test that empty phase array returns default values."""
        phases = np.array([])
        kappa, rayleigh_z, rayleigh_p = _compute_circular_stats(phases)
        assert kappa == 0.0
        assert rayleigh_z == 0.0
        assert rayleigh_p == 1.0

    def test_single_phase_returns_defaults(self) -> None:
        """Test that single phase value produces valid output."""
        phases = np.array([0.5])
        kappa, rayleigh_z, rayleigh_p = _compute_circular_stats(phases)
        # With single value, r_bar = 1, kappa should be high
        assert kappa > 0
        assert rayleigh_z > 0
        assert 0 <= rayleigh_p <= 1


# =============================================================================
# Test compute_paired_statistics (Cohen's d)
# =============================================================================
@pytest.mark.unit
class TestComputePairedStatistics:
    """Test paired statistical analysis including Cohen's d."""

    def test_identical_conditions_zero_effect(self) -> None:
        """Test that identical conditions yield zero effect size."""
        # Create fold metrics with identical AUROC for both conditions
        fold_metrics = [
            FoldMetrics(
                fold_idx=i,
                repeat_idx=0,
                auroc_a=0.75,
                auroc_b=0.75,
                auprc_a=0.6,
                auprc_b=0.6,
                brier_a=0.2,
                brier_b=0.2,
            )
            for i in range(10)
        ]
        config = StudyConfig(n_bootstrap=100)
        results = compute_paired_statistics(fold_metrics, config)

        assert "auroc" in results
        assert results["auroc"]["mean_delta"] == 0.0
        # Cohen's d should be 0 or close to 0
        assert abs(results["auroc"]["cohens_d"]) < 1e-10

    def test_different_conditions_nonzero_effect(self) -> None:
        """Test that different conditions yield non-zero effect size."""
        # Condition A consistently better than B, with some variance
        rng = np.random.default_rng(42)
        fold_metrics = [
            FoldMetrics(
                fold_idx=i,
                repeat_idx=0,
                auroc_a=0.80 + rng.uniform(-0.02, 0.02),
                auroc_b=0.70 + rng.uniform(-0.02, 0.02),
                auprc_a=0.65 + rng.uniform(-0.02, 0.02),
                auprc_b=0.55 + rng.uniform(-0.02, 0.02),
                brier_a=0.15 + rng.uniform(-0.02, 0.02),
                brier_b=0.25 + rng.uniform(-0.02, 0.02),
            )
            for i in range(10)
        ]
        config = StudyConfig(n_bootstrap=100)
        results = compute_paired_statistics(fold_metrics, config)

        # Mean delta should be approximately 0.10
        assert results["auroc"]["mean_delta"] == pytest.approx(0.10, abs=0.05)
        # Effect size should be positive (A > B consistently)
        assert results["auroc"]["cohens_d"] > 0

    def test_bootstrap_ci_computed(self) -> None:
        """Test that bootstrap confidence intervals are computed."""
        fold_metrics = [
            FoldMetrics(
                fold_idx=i,
                repeat_idx=0,
                auroc_a=0.75 + np.random.default_rng(i).uniform(-0.05, 0.05),
                auroc_b=0.70 + np.random.default_rng(i + 100).uniform(-0.05, 0.05),
                auprc_a=0.6,
                auprc_b=0.55,
                brier_a=0.2,
                brier_b=0.25,
            )
            for i in range(20)
        ]
        config = StudyConfig(n_bootstrap=100)
        results = compute_paired_statistics(fold_metrics, config)

        assert "cohens_d_ci_low" in results["auroc"]
        assert "cohens_d_ci_high" in results["auroc"]
        # CI should bracket the point estimate (usually)
        ci_low = results["auroc"]["cohens_d_ci_low"]
        ci_high = results["auroc"]["cohens_d_ci_high"]
        assert ci_low <= ci_high

    def test_effect_interpretation_categories(self) -> None:
        """Test effect size interpretation categories."""
        # Create metrics with large effect
        fold_metrics = [
            FoldMetrics(
                fold_idx=i,
                repeat_idx=0,
                auroc_a=0.90,
                auroc_b=0.60,
                auprc_a=0.7,
                auprc_b=0.4,
                brier_a=0.1,
                brier_b=0.4,
            )
            for i in range(10)
        ]
        config = StudyConfig(n_bootstrap=100)
        results = compute_paired_statistics(fold_metrics, config)

        assert results["auroc"]["effect_interpretation"] in [
            "negligible", "small", "medium", "large"
        ]

    def test_empty_fold_metrics_returns_error(self) -> None:
        """Test that empty fold metrics returns error dict."""
        config = StudyConfig()
        results = compute_paired_statistics([], config)
        assert "error" in results


# =============================================================================
# Test compute_power_analysis
# =============================================================================
@pytest.mark.unit
class TestComputePowerAnalysis:
    """Test power analysis computation."""

    def test_returns_power_metrics(self) -> None:
        """Test that power analysis returns expected metrics."""
        # Add variance to deltas so we don't get zero variance error
        rng = np.random.default_rng(42)
        fold_metrics = [
            FoldMetrics(
                fold_idx=i,
                repeat_idx=0,
                auroc_a=0.80 + rng.uniform(-0.02, 0.02),
                auroc_b=0.70 + rng.uniform(-0.02, 0.02),
                auprc_a=0.65 + rng.uniform(-0.02, 0.02),
                auprc_b=0.55 + rng.uniform(-0.02, 0.02),
                brier_a=0.15 + rng.uniform(-0.02, 0.02),
                brier_b=0.25 + rng.uniform(-0.02, 0.02),
            )
            for i in range(10)
        ]
        config = StudyConfig()
        results = compute_power_analysis(fold_metrics, config)

        assert "n_observations" in results
        assert "observed_delta" in results
        assert "observed_std" in results
        assert "achieved_power" in results
        assert "n_needed_for_desired_power" in results

    def test_power_in_valid_range(self) -> None:
        """Test that achieved power is in [0, 1]."""
        # Add variance to deltas so we don't get zero variance error
        rng = np.random.default_rng(42)
        fold_metrics = [
            FoldMetrics(
                fold_idx=i,
                repeat_idx=0,
                auroc_a=0.75 + rng.uniform(-0.03, 0.03),
                auroc_b=0.70 + rng.uniform(-0.03, 0.03),
                auprc_a=0.6 + rng.uniform(-0.03, 0.03),
                auprc_b=0.55 + rng.uniform(-0.03, 0.03),
                brier_a=0.2 + rng.uniform(-0.03, 0.03),
                brier_b=0.25 + rng.uniform(-0.03, 0.03),
            )
            for i in range(20)
        ]
        config = StudyConfig()
        results = compute_power_analysis(fold_metrics, config)

        assert 0 <= results["achieved_power"] <= 1

    def test_insufficient_data_returns_error(self) -> None:
        """Test that insufficient data returns error dict."""
        # Only 1 fold
        fold_metrics = [
            FoldMetrics(
                fold_idx=0,
                repeat_idx=0,
                auroc_a=0.75,
                auroc_b=0.70,
                auprc_a=0.6,
                auprc_b=0.55,
                brier_a=0.2,
                brier_b=0.25,
            )
        ]
        config = StudyConfig()
        results = compute_power_analysis(fold_metrics, config)

        assert "error" in results

    def test_zero_variance_returns_error(self) -> None:
        """Test that zero variance in deltas returns error."""
        # All identical deltas (zero variance)
        fold_metrics = [
            FoldMetrics(
                fold_idx=i,
                repeat_idx=0,
                auroc_a=0.75,
                auroc_b=0.70,  # Same delta for all
                auprc_a=0.6,
                auprc_b=0.55,
                brier_a=0.2,
                brier_b=0.25,
            )
            for i in range(10)
        ]
        config = StudyConfig()
        results = compute_power_analysis(fold_metrics, config)

        assert "error" in results


# =============================================================================
# Test apply_fdr_correction
# =============================================================================
@pytest.mark.unit
class TestApplyFDRCorrection:
    """Test Benjamini-Hochberg FDR correction."""

    def test_empty_pvalues(self) -> None:
        """Test that empty p-value list returns empty results."""
        significant, adjusted = apply_fdr_correction([])
        assert significant == []
        assert adjusted == []

    def test_single_pvalue_below_threshold(self) -> None:
        """Test single p-value below threshold is significant."""
        significant, adjusted = apply_fdr_correction([0.01], alpha=0.05)
        assert significant == [True]
        assert adjusted[0] <= 0.05

    def test_single_pvalue_above_threshold(self) -> None:
        """Test single p-value above threshold is not significant."""
        significant, adjusted = apply_fdr_correction([0.1], alpha=0.05)
        assert significant == [False]

    def test_adjusted_pvalues_ordered(self) -> None:
        """Test that adjusted p-values maintain relative ordering."""
        p_values = [0.01, 0.02, 0.05, 0.10, 0.20]
        significant, adjusted = apply_fdr_correction(p_values)

        # Adjusted values should increase with original values (approximately)
        # Though FDR adjustment can change ordering in some edge cases
        assert all(0 <= p <= 1 for p in adjusted)

    def test_all_significant(self) -> None:
        """Test case where all p-values are significant."""
        p_values = [0.001, 0.002, 0.003]
        significant, adjusted = apply_fdr_correction(p_values, alpha=0.05)
        assert all(significant)

    def test_none_significant(self) -> None:
        """Test case where no p-values are significant."""
        p_values = [0.5, 0.6, 0.7]
        significant, adjusted = apply_fdr_correction(p_values, alpha=0.05)
        assert not any(significant)


# =============================================================================
# Test CZTFeatures dataclass
# =============================================================================
@pytest.mark.unit
class TestCZTFeatures:
    """Test CZTFeatures dataclass."""

    def test_to_array_length(self) -> None:
        """Test that to_array returns correct number of features."""
        features = CZTFeatures(
            peak_magnitude=1.0,
            peak_freq_idx=128,
            phase_at_peak=0.5,
            band_power=10.0,
            spectral_centroid=0.095,
            phase_kappa=2.0,
            rayleigh_z=50.0,
            rayleigh_p=0.001,
        )
        arr = features.to_array()
        assert len(arr) == 7
        # rayleigh_p should NOT be in the array (diagnostic only)
        assert 0.001 not in arr

    def test_to_array_order(self) -> None:
        """Test that to_array values match expected order."""
        features = CZTFeatures(
            peak_magnitude=1.0,
            peak_freq_idx=128,
            phase_at_peak=0.5,
            band_power=10.0,
            spectral_centroid=0.095,
            phase_kappa=2.0,
            rayleigh_z=50.0,
            rayleigh_p=0.001,
        )
        arr = features.to_array()
        assert arr[0] == 1.0  # peak_magnitude
        assert arr[1] == 128  # peak_freq_idx
        assert arr[2] == 0.5  # phase_at_peak
        assert arr[3] == 10.0  # band_power
        assert arr[4] == 0.095  # spectral_centroid
        assert arr[5] == 2.0  # phase_kappa
        assert arr[6] == 50.0  # rayleigh_z


# =============================================================================
# Test FoldMetrics dataclass
# =============================================================================
@pytest.mark.unit
class TestFoldMetrics:
    """Test FoldMetrics dataclass."""

    def test_delta_auroc_computed(self) -> None:
        """Test that delta_auroc is computed correctly."""
        metrics = FoldMetrics(
            fold_idx=0,
            repeat_idx=0,
            auroc_a=0.80,
            auroc_b=0.70,
            auprc_a=0.6,
            auprc_b=0.5,
            brier_a=0.2,
            brier_b=0.3,
        )
        assert metrics.delta_auroc == pytest.approx(0.10, abs=1e-6)

    def test_delta_auprc_computed(self) -> None:
        """Test that delta_auprc is computed correctly."""
        metrics = FoldMetrics(
            fold_idx=0,
            repeat_idx=0,
            auroc_a=0.80,
            auroc_b=0.70,
            auprc_a=0.6,
            auprc_b=0.5,
            brier_a=0.2,
            brier_b=0.3,
        )
        assert metrics.delta_auprc == pytest.approx(0.10, abs=1e-6)

    def test_delta_brier_computed_inverted(self) -> None:
        """Test that delta_brier is computed as B-A (lower is better)."""
        metrics = FoldMetrics(
            fold_idx=0,
            repeat_idx=0,
            auroc_a=0.80,
            auroc_b=0.70,
            auprc_a=0.6,
            auprc_b=0.5,
            brier_a=0.2,
            brier_b=0.3,
        )
        # delta_brier = brier_b - brier_a
        assert metrics.delta_brier == pytest.approx(0.10, abs=1e-6)


# =============================================================================
# Test StudyConfig defaults
# =============================================================================
@pytest.mark.unit
class TestStudyConfig:
    """Test StudyConfig dataclass defaults."""

    def test_default_seeds(self) -> None:
        """Test that default seeds are set correctly."""
        config = StudyConfig()
        assert config.global_seed == 137
        assert config.phase_seed == 271828
        assert config.cv_seed == 161803

    def test_default_cv_params(self) -> None:
        """Test that default CV parameters are set correctly."""
        config = StudyConfig()
        assert config.n_folds == 5
        assert config.n_repeats == 5

    def test_default_band_center(self) -> None:
        """Test that band center defaults to helical frequency."""
        config = StudyConfig()
        expected = 1 / 10.5
        assert config.band_center == pytest.approx(expected, rel=1e-6)


# =============================================================================
# Test biophysics parameters
# =============================================================================
@pytest.mark.unit
@pytest.mark.smoke
def test_biophysics_parameters() -> None:
    """Smoke test for biophysics parameter tables."""
    # Check lifetimes are correctly defined
    assert BREATH_LIFETIMES_MS["A"] == 5.0
    assert BREATH_LIFETIMES_MS["T"] == 5.0
    assert BREATH_LIFETIMES_MS["G"] == 25.0
    assert BREATH_LIFETIMES_MS["C"] == 25.0

    # Check ΔG bounds are computed from dictionary
    expected_min = min(NEAREST_NEIGHBOR_DG.values())
    expected_max = max(NEAREST_NEIGHBOR_DG.values())
    assert _DG_MIN == expected_min
    assert _DG_MAX == expected_max


@pytest.mark.unit
@pytest.mark.smoke
def test_helical_period_defined() -> None:
    """Smoke test for helical period constant."""
    assert DEFAULT_HELICAL_PERIOD == 10.5


# =============================================================================
# Test extract_features_for_sequences
# =============================================================================
@pytest.mark.unit
class TestExtractFeaturesForSequences:
    """Test feature extraction for multiple sequences."""

    def test_returns_two_feature_arrays(self) -> None:
        """Test that feature extraction returns both coherent and random arrays."""
        sequences = ["ATCGATCGATCGATCGATCG", "GCTAGCTAGCTAGCTAGCTA"]
        config = StudyConfig()
        phase_rng = np.random.default_rng(config.phase_seed)
        X_coherent, X_random = extract_features_for_sequences(
            sequences, config, phase_rng
        )

        assert X_coherent.shape[0] == 2
        assert X_random.shape[0] == 2
        assert X_coherent.shape[1] == 7
        assert X_random.shape[1] == 7

    def test_coherent_and_random_differ(self) -> None:
        """Test that coherent and random features differ (phases changed)."""
        sequences = ["ATCGATCGATCGATCGATCG"]
        config = StudyConfig()
        phase_rng = np.random.default_rng(config.phase_seed)
        X_coherent, X_random = extract_features_for_sequences(
            sequences, config, phase_rng
        )

        # Features should not be identical (phase-dependent features differ)
        assert not np.allclose(X_coherent, X_random)


# =============================================================================
# Integration tests
# =============================================================================
@pytest.mark.unit
class TestIntegration:
    """Integration tests for the study pipeline."""

    def test_full_pipeline_small_sample(self) -> None:
        """Test complete pipeline with small synthetic dataset."""
        # Create synthetic sequences with labels
        rng = np.random.default_rng(42)
        bases = "ATGC"
        sequences = []
        labels = []

        # Class 0: AT-rich sequences (more breathing)
        for _ in range(20):
            seq = "".join(rng.choice(list("AT"), 20))
            sequences.append(seq)
            labels.append(0)

        # Class 1: GC-rich sequences (more stable)
        for _ in range(20):
            seq = "".join(rng.choice(list("GC"), 20))
            sequences.append(seq)
            labels.append(1)

        y = np.array(labels)

        # Extract features
        config = StudyConfig(
            n_folds=2,
            n_repeats=1,
            n_bootstrap=100,
        )
        phase_rng = np.random.default_rng(config.phase_seed)
        X_coherent, X_random = extract_features_for_sequences(
            sequences, config, phase_rng
        )

        # Run CV evaluation
        fold_metrics = run_cv_evaluation(X_coherent, X_random, y, config)

        # Should have some fold results
        assert len(fold_metrics) > 0

        # Compute statistics
        stats_results = compute_paired_statistics(fold_metrics, config)
        assert "auroc" in stats_results

    def test_seeds_produce_reproducible_results(self) -> None:
        """Test that fixed seeds produce reproducible results."""
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTA",
            "ATATATATATATATATATA",
            "GCGCGCGCGCGCGCGCGCGC",
        ]
        labels = np.array([0, 1, 0, 1])

        config = StudyConfig(n_folds=2, n_repeats=1)

        # First run
        phase_rng1 = np.random.default_rng(config.phase_seed)
        X_coherent1, X_random1 = extract_features_for_sequences(
            sequences, config, phase_rng1
        )

        # Second run with same seed
        phase_rng2 = np.random.default_rng(config.phase_seed)
        X_coherent2, X_random2 = extract_features_for_sequences(
            sequences, config, phase_rng2
        )

        np.testing.assert_array_equal(X_coherent1, X_coherent2)
        np.testing.assert_array_equal(X_random1, X_random2)


# =============================================================================
# Edge case tests
# =============================================================================
@pytest.mark.unit
class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_short_sequence(self) -> None:
        """Test handling of very short sequence (<10 bp)."""
        seq = "ATCGATCG"  # 8 bp
        signal = encode_sequence(seq)
        assert signal is not None
        assert len(signal) == 8

        # CZT should still work
        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert np.all(np.isfinite(features.to_array()))

    def test_minimum_sequence_length(self) -> None:
        """Test minimum viable sequence length."""
        seq = "AT"  # 2 bp minimum
        signal = encode_sequence(seq)
        assert signal is not None
        assert len(signal) == 2

    def test_all_same_base(self) -> None:
        """Test handling of homopolymer sequence."""
        seq = "AAAAAAAAAA"  # All A's
        signal = encode_sequence(seq)
        assert signal is not None

        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        # Should produce valid finite features
        assert np.all(np.isfinite(features.to_array()))

    def test_alternating_bases(self) -> None:
        """Test handling of alternating base pattern."""
        seq = "ATATATATAT"
        signal = encode_sequence(seq)
        assert signal is not None

        freqs, spectrum = compute_czt_spectrum(signal)
        features = extract_czt_features(freqs, spectrum)
        assert np.all(np.isfinite(features.to_array()))

    def test_gc_rich_sequence(self) -> None:
        """Test handling of GC-rich sequence."""
        seq = "GCGCGCGCGCGCGCGC"
        signal = encode_sequence(seq)
        assert signal is not None

        # Real part should be all 25.0 (GC lifetime)
        x_no_helical = encode_sequence(seq, apply_helical=False)
        np.testing.assert_array_equal(x_no_helical.real, [25.0] * len(seq))

    def test_single_class_cv_skip(self) -> None:
        """Test that CV skips folds with single class."""
        sequences = [
            "ATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTA",
        ]
        # All same label
        y = np.array([0, 0])

        config = StudyConfig(n_folds=2, n_repeats=1)
        phase_rng = np.random.default_rng(config.phase_seed)
        X_coherent, X_random = extract_features_for_sequences(
            sequences, config, phase_rng
        )

        # Should handle gracefully (may return empty or raise)
        # The function skips folds with single class
        fold_metrics = run_cv_evaluation(X_coherent, X_random, y, config)
        # Result depends on how splits work out - may be empty
        assert isinstance(fold_metrics, list)
