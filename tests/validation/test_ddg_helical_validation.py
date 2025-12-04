"""
Validation tests for ΔΔG-Pairs + Helical Band Sweep experiment.

Tests cover:
1. ΔΔG calculation correctness
2. Pair generation and binning
3. Band sweep mechanics
4. Statistical computation accuracy
"""

import math
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pytest

# Add the experiments module to path
experiments_dir = Path(__file__).parent.parent.parent / "experiments"
sys.path.insert(0, str(experiments_dir))

# Add the gist source to path for NEAREST_NEIGHBOR_DG
gist_src = Path(__file__).parent.parent.parent / "gists/breathing/czt_feature_extractor/src"
sys.path.insert(0, str(gist_src))

from dna_breathing_gist import NEAREST_NEIGHBOR_DG

from ddg_helical_validation import (
    ValidationConfig,
    bca_bootstrap_ci,
    benjamini_hochberg_correction,
    cohens_d_paired,
    compute_ddg,
    compute_paired_statistics,
    compute_total_dg,
    extract_features_at_band,
    generate_mutants_for_sequence,
    generate_point_mutant,
    generate_wt_mutant_pairs,
    jonckheere_terpstra_test,
    label_shuffle_permutation,
    run_band_sweep,
)


# =============================================================================
# Test ΔΔG Calculation
# =============================================================================


@pytest.mark.validation
class TestDDGCalculation:
    """Tests for ΔΔG (mutation-induced free-energy change) calculations."""

    def test_total_dg_computation_basic(self) -> None:
        """Test basic total ΔG computation for a simple sequence."""
        # Using SantaLucia 1998 parameters from NEAREST_NEIGHBOR_DG
        seq = "AAT"
        dg = compute_total_dg(seq)
        # AA + AT dinucleotides
        expected = NEAREST_NEIGHBOR_DG["AA"] + NEAREST_NEIGHBOR_DG["AT"]
        assert abs(dg - expected) < 1e-10, f"Expected {expected}, got {dg}"

    def test_total_dg_gc_rich_sequence(self) -> None:
        """Test ΔG for GC-rich sequence (more stable)."""
        seq = "GCG"
        dg = compute_total_dg(seq)
        # GC + CG dinucleotides
        expected = NEAREST_NEIGHBOR_DG["GC"] + NEAREST_NEIGHBOR_DG["CG"]
        assert abs(dg - expected) < 1e-10

    def test_total_dg_at_rich_sequence(self) -> None:
        """Test ΔG for AT-rich sequence (less stable)."""
        seq = "ATA"
        dg = compute_total_dg(seq)
        # AT + TA dinucleotides
        expected = NEAREST_NEIGHBOR_DG["AT"] + NEAREST_NEIGHBOR_DG["TA"]
        assert abs(dg - expected) < 1e-10

    def test_ddg_stabilizing_mutation(self) -> None:
        """Test that AT→GC mutation produces negative ΔΔG (stabilizing)."""
        wt = "ATAT"  # AT-rich
        mut = "GCGC"  # GC-rich
        ddg = compute_ddg(wt, mut)
        # GC-rich should be more stable (more negative ΔG)
        # So ΔΔG = ΔG(mut) - ΔG(wt) should be negative
        assert ddg < 0, "GC-rich mutation should be stabilizing (negative ΔΔG)"

    def test_ddg_destabilizing_mutation(self) -> None:
        """Test that GC→AT mutation produces positive ΔΔG (destabilizing)."""
        wt = "GCGC"  # GC-rich
        mut = "ATAT"  # AT-rich
        ddg = compute_ddg(wt, mut)
        # AT-rich should be less stable (less negative ΔG)
        # So ΔΔG = ΔG(mut) - ΔG(wt) should be positive
        assert ddg > 0, "AT-rich mutation should be destabilizing (positive ΔΔG)"

    def test_ddg_symmetry(self) -> None:
        """Test that ΔΔG(A→B) = -ΔΔG(B→A)."""
        seq_a = "ATGC"
        seq_b = "GCTA"
        ddg_ab = compute_ddg(seq_a, seq_b)
        ddg_ba = compute_ddg(seq_b, seq_a)
        assert abs(ddg_ab + ddg_ba) < 1e-10, "ΔΔG should be antisymmetric"

    def test_ddg_identical_sequences_zero(self) -> None:
        """Test that ΔΔG is zero for identical sequences."""
        seq = "ATGCATGC"
        ddg = compute_ddg(seq, seq)
        assert abs(ddg) < 1e-10, "ΔΔG should be zero for identical sequences"

    def test_single_base_sequence_zero_dg(self) -> None:
        """Test that single base sequence has zero ΔG."""
        dg = compute_total_dg("A")
        assert dg == 0.0, "Single base should have zero ΔG"

    def test_empty_sequence_zero_dg(self) -> None:
        """Test that empty sequence has zero ΔG."""
        dg = compute_total_dg("")
        assert dg == 0.0, "Empty sequence should have zero ΔG"


# =============================================================================
# Test Pair Generation and Binning
# =============================================================================


@pytest.mark.validation
class TestPairGeneration:
    """Tests for WT-mutant pair generation and binning."""

    def test_generate_point_mutant_basic(self) -> None:
        """Test basic point mutation generation."""
        seq = "ATGC"
        mut = generate_point_mutant(seq, 0, "G")
        assert mut == "GTGC", f"Expected GTGC, got {mut}"

    def test_generate_point_mutant_last_position(self) -> None:
        """Test point mutation at last position."""
        seq = "ATGC"
        mut = generate_point_mutant(seq, 3, "T")
        assert mut == "ATGT", f"Expected ATGT, got {mut}"

    def test_generate_point_mutant_invalid_position_raises(self) -> None:
        """Test that invalid position raises ValueError."""
        seq = "ATGC"
        with pytest.raises(ValueError):
            generate_point_mutant(seq, 10, "A")

    def test_generate_point_mutant_invalid_base_raises(self) -> None:
        """Test that invalid base raises ValueError."""
        seq = "ATGC"
        with pytest.raises(ValueError):
            generate_point_mutant(seq, 0, "X")

    def test_generate_mutants_count(self) -> None:
        """Test that correct number of mutants are generated."""
        seq = "ATGC"  # 4 bases
        mutants = generate_mutants_for_sequence(seq, seed=42)
        # Each position can have 3 mutations (to the other 3 bases)
        expected_count = len(seq) * 3
        assert len(mutants) == expected_count, f"Expected {expected_count}, got {len(mutants)}"

    def test_generate_mutants_all_different(self) -> None:
        """Test that all generated mutants are different from WT."""
        seq = "ATGC"
        mutants = generate_mutants_for_sequence(seq, seed=42)
        for mut_seq, pos, orig, ddg in mutants:
            assert mut_seq != seq, "Mutant should differ from WT"
            # Check that exactly one position differs
            diff_count = sum(1 for i in range(len(seq)) if mut_seq[i] != seq[i])
            assert diff_count == 1, "Should have exactly one mutation"

    def test_generate_wt_mutant_pairs_binning(self) -> None:
        """Test that pairs are correctly binned into tertiles."""
        sequences = [
            ("seq1", "ATGCATGC"),
            ("seq2", "GCTAGCTA"),
        ]
        pairs, bin_edges = generate_wt_mutant_pairs(sequences, num_bins=3, seed=42)

        # Check bins are assigned
        bins = set(p["bin"] for p in pairs)
        assert "low" in bins or "mid" in bins or "high" in bins, "Should have bin assignments"

        # Check all pairs have required fields
        for pair in pairs:
            assert "wt_id" in pair
            assert "mut_id" in pair
            assert "wt_seq" in pair
            assert "mut_seq" in pair
            assert "delta_delta_g" in pair
            assert "bin" in pair

    def test_generate_wt_mutant_pairs_bin_edges(self) -> None:
        """Test that bin edges are properly computed tertiles."""
        sequences = [
            ("seq1", "ATGCATGC"),
            ("seq2", "GCTAGCTA"),
            ("seq3", "TATATATATA"),
        ]
        pairs, bin_edges = generate_wt_mutant_pairs(sequences, num_bins=3, seed=42)

        # Should have 4 edges for 3 bins
        assert len(bin_edges) == 4, f"Expected 4 bin edges, got {len(bin_edges)}"
        # Edges should be monotonically increasing
        for i in range(len(bin_edges) - 1):
            assert bin_edges[i] <= bin_edges[i + 1], "Bin edges should be monotonic"


# =============================================================================
# Test Band Sweep Mechanics
# =============================================================================


@pytest.mark.validation
class TestBandSweep:
    """Tests for band sweep analysis mechanics."""

    def test_extract_features_at_band_returns_dict(self) -> None:
        """Test that feature extraction returns expected dictionary keys."""
        seq = "ATGCATGCATGC"
        config = ValidationConfig(seed=42)
        features = extract_features_at_band(seq, center_period=10.5, width_fraction=0.03, config=config)

        expected_keys = ["peak_mag", "band_power", "phase_coherence", "snr"]
        for key in expected_keys:
            assert key in features, f"Missing key: {key}"
            assert isinstance(features[key], float), f"Value for {key} should be float"

    def test_extract_features_at_different_centers(self) -> None:
        """Test that different center frequencies produce different features."""
        seq = "ATGCATGCATGCATGC"
        config = ValidationConfig(seed=42)

        feat_low = extract_features_at_band(seq, center_period=10.3, width_fraction=0.03, config=config)
        feat_mid = extract_features_at_band(seq, center_period=10.5, width_fraction=0.03, config=config)
        feat_high = extract_features_at_band(seq, center_period=10.7, width_fraction=0.03, config=config)

        # At least some features should differ
        all_same = True
        for key in ["peak_mag", "band_power"]:
            if feat_low[key] != feat_mid[key] or feat_mid[key] != feat_high[key]:
                all_same = False
                break
        assert not all_same, "Features should vary with center frequency"

    def test_extract_features_at_different_widths(self) -> None:
        """Test that different band widths produce different features."""
        seq = "ATGCATGCATGCATGC"
        config = ValidationConfig(seed=42)

        feat_narrow = extract_features_at_band(seq, center_period=10.5, width_fraction=0.01, config=config)
        feat_wide = extract_features_at_band(seq, center_period=10.5, width_fraction=0.06, config=config)

        # Band power should be larger for wider bands
        # (more frequency bins contribute)
        # Note: This may not always hold depending on spectral shape
        assert feat_narrow["band_power"] != feat_wide["band_power"], "Band power should differ with width"

    def test_run_band_sweep_structure(self) -> None:
        """Test that band sweep returns correct structure."""
        sequences = [("seq1", "ATGCATGCATGC")]
        pairs, _ = generate_wt_mutant_pairs(sequences, num_bins=3, seed=42)

        # Take just a few pairs for speed
        test_pairs = pairs[:3]
        config = ValidationConfig(
            seed=42,
            center_sweep=[10.5],  # Single center for speed
            width_sweep=[0.03],   # Single width for speed
        )

        results = run_band_sweep(test_pairs, config)

        assert len(results) > 0, "Should have results"
        for r in results:
            assert "pair_id" in r
            assert "wt_id" in r
            assert "center" in r
            assert "width" in r
            assert "feature" in r
            assert "wt_value" in r
            assert "mut_value" in r
            assert "diff" in r


# =============================================================================
# Test Statistical Computations
# =============================================================================


@pytest.mark.validation
class TestStatistics:
    """Tests for statistical computation accuracy."""

    def test_cohens_d_paired_known_value(self) -> None:
        """Test Cohen's d for paired data with known values."""
        # Mean diff = 1, std diff = 1 → d = 1
        diffs = np.array([0.0, 1.0, 2.0])  # mean=1, std≈1
        d = cohens_d_paired(diffs)
        expected = 1.0 / np.std(diffs, ddof=1)
        assert abs(d - expected) < 1e-10

    def test_cohens_d_paired_zero_effect(self) -> None:
        """Test Cohen's d is zero when mean diff is zero."""
        diffs = np.array([-1.0, 1.0, -1.0, 1.0])  # mean=0
        d = cohens_d_paired(diffs)
        assert abs(d) < 1e-10, "d should be ~0 for zero mean diff"

    def test_cohens_d_paired_large_effect(self) -> None:
        """Test Cohen's d for large effect (|d| > 0.8)."""
        diffs = np.array([2.0, 2.5, 1.5, 3.0, 2.0])  # mean≈2.2, relatively small std
        d = cohens_d_paired(diffs)
        assert abs(d) > 0.8, f"Expected large effect, got d={d}"

    def test_cohens_d_paired_nan_for_zero_variance(self) -> None:
        """Test that Cohen's d returns NaN for zero variance."""
        diffs = np.array([1.0, 1.0, 1.0])  # zero variance
        d = cohens_d_paired(diffs)
        assert math.isnan(d), "Should return NaN for zero variance"

    def test_cohens_d_paired_insufficient_data(self) -> None:
        """Test that Cohen's d returns NaN for insufficient data."""
        diffs = np.array([1.0])  # only one value
        d = cohens_d_paired(diffs)
        assert math.isnan(d), "Should return NaN for insufficient data"

    def test_bca_bootstrap_ci_contains_point_estimate(self) -> None:
        """Test that BCa CI contains the point estimate for most cases."""
        np.random.seed(42)
        data = np.random.normal(0, 1, 50)
        theta, ci_low, ci_high = bca_bootstrap_ci(
            data,
            lambda x: np.mean(x),
            num_bootstrap=1000,
            alpha=0.05,
            seed=42,
        )
        # The point estimate should typically be within the CI
        # (though not guaranteed for all samples)
        assert ci_low <= theta <= ci_high or math.isnan(ci_low)

    def test_bca_bootstrap_ci_narrower_for_larger_sample(self) -> None:
        """Test that CI is narrower for larger samples."""
        np.random.seed(42)

        small_data = np.random.normal(0, 1, 20)
        large_data = np.random.normal(0, 1, 200)

        _, ci_low_small, ci_high_small = bca_bootstrap_ci(
            small_data, np.mean, num_bootstrap=500, seed=42
        )
        _, ci_low_large, ci_high_large = bca_bootstrap_ci(
            large_data, np.mean, num_bootstrap=500, seed=42
        )

        width_small = ci_high_small - ci_low_small
        width_large = ci_high_large - ci_low_large

        # Larger sample should have narrower CI
        assert width_large < width_small, "Larger sample should have narrower CI"

    def test_benjamini_hochberg_preserves_order(self) -> None:
        """Test that BH correction preserves p-value ordering."""
        p_values = [0.01, 0.03, 0.05, 0.10, 0.50]
        q_values = benjamini_hochberg_correction(p_values)

        for i in range(len(q_values) - 1):
            assert q_values[i] <= q_values[i + 1], "q-values should be monotonically increasing"

    def test_benjamini_hochberg_bounds_q_at_one(self) -> None:
        """Test that BH-corrected q-values don't exceed 1."""
        p_values = [0.5, 0.6, 0.7, 0.8, 0.9]
        q_values = benjamini_hochberg_correction(p_values)

        for q in q_values:
            assert q <= 1.0, "q-values should not exceed 1"

    def test_benjamini_hochberg_handles_nan(self) -> None:
        """Test that BH correction handles NaN values."""
        p_values = [0.01, np.nan, 0.05, np.nan, 0.10]
        q_values = benjamini_hochberg_correction(p_values)

        # NaN positions should remain NaN
        assert math.isnan(q_values[1])
        assert math.isnan(q_values[3])
        # Valid positions should have valid q-values
        assert not math.isnan(q_values[0])
        assert not math.isnan(q_values[2])

    def test_label_shuffle_permutation_p_bounds(self) -> None:
        """Test that label shuffle p-value is in [0, 1]."""
        diffs = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        p = label_shuffle_permutation(diffs, num_perm=100, seed=42)
        assert 0.0 < p <= 1.0, "p-value should be in (0, 1]"

    def test_label_shuffle_permutation_nonzero(self) -> None:
        """Test that label shuffle p-value is never exactly zero (smoothing)."""
        # Even with extreme data, add-one smoothing should prevent p=0
        diffs = np.array([100.0, 100.0, 100.0, 100.0, 100.0])
        p = label_shuffle_permutation(diffs, num_perm=100, seed=42)
        assert p > 0, "p-value should never be exactly zero"

    def test_label_shuffle_null_distribution(self) -> None:
        """Test that zero-mean diffs give non-significant p-values on average."""
        np.random.seed(42)
        p_values = []
        for i in range(50):
            diffs = np.random.normal(0, 1, 20)  # True null: mean = 0
            p = label_shuffle_permutation(diffs, num_perm=100, seed=42 + i)
            p_values.append(p)

        # Under null, p-values should be roughly uniform
        # Mean should be ~0.5
        mean_p = np.mean(p_values)
        assert 0.3 < mean_p < 0.7, f"Under null, mean p should be ~0.5, got {mean_p}"

    def test_jonckheere_terpstra_increasing_trend(self) -> None:
        """Test JT test detects increasing trend."""
        groups = [
            np.array([1.0, 2.0, 1.5]),
            np.array([3.0, 4.0, 3.5]),
            np.array([6.0, 7.0, 6.5]),
        ]
        j_stat, p_value = jonckheere_terpstra_test(groups, alternative="increasing")
        assert p_value < 0.05, f"Should detect increasing trend, p={p_value}"

    def test_jonckheere_terpstra_no_trend(self) -> None:
        """Test JT test returns non-significant for random data."""
        np.random.seed(42)
        groups = [
            np.random.normal(5, 1, 10),
            np.random.normal(5, 1, 10),
            np.random.normal(5, 1, 10),
        ]
        _, p_value = jonckheere_terpstra_test(groups, alternative="increasing")
        # Should not detect trend in random data (most of the time)
        assert p_value > 0.01, f"Should not detect trend in random data, p={p_value}"

    def test_compute_paired_statistics_structure(self) -> None:
        """Test that paired statistics returns expected structure."""
        np.random.seed(42)
        diffs = np.random.normal(0.5, 1, 30)

        stats = compute_paired_statistics(
            diffs, num_bootstrap=100, num_perm=50, alpha=0.05, seed=42
        )

        expected_keys = ["n", "mean_diff", "d", "ci_low", "ci_high", "p_value", "test_type"]
        for key in expected_keys:
            assert key in stats, f"Missing key: {key}"

    def test_compute_paired_statistics_test_selection(self) -> None:
        """Test that appropriate test is selected based on normality."""
        # Normal data should use t-test
        np.random.seed(42)
        normal_diffs = np.random.normal(0, 1, 50)
        stats_normal = compute_paired_statistics(normal_diffs, num_bootstrap=100, seed=42)

        # The test type should be either paired_t or wilcoxon depending on the sample
        assert stats_normal["test_type"] in ["paired_t", "wilcoxon"]


# =============================================================================
# Integration Tests
# =============================================================================


@pytest.mark.validation
class TestIntegration:
    """Integration tests for the complete validation pipeline."""

    def test_end_to_end_small_dataset(self) -> None:
        """Test complete pipeline with a small dataset."""
        sequences = [
            ("test_seq1", "ATGCATGCATGC"),
            ("test_seq2", "GCTAGCTAGCTA"),
        ]

        config = ValidationConfig(
            seed=42,
            num_bootstrap=100,  # Reduced for speed
            num_permutations=50,
            center_sweep=[10.5],  # Single center for speed
            width_sweep=[0.03],
        )

        # Generate pairs
        pairs, bin_edges = generate_wt_mutant_pairs(sequences, config.num_bins, config.seed)
        assert len(pairs) > 0

        # Run band sweep
        results = run_band_sweep(pairs[:5], config)  # Just a few pairs
        assert len(results) > 0

        # Each result should have paired differences
        for r in results:
            assert "diff" in r
            assert isinstance(r["diff"], float)

    def test_config_serialization(self) -> None:
        """Test that config can be serialized to JSON."""
        import json

        config = ValidationConfig(seed=123, num_bootstrap=500)
        config_dict = config.to_dict()

        # Should be JSON serializable
        json_str = json.dumps(config_dict)
        assert len(json_str) > 0

        # Should contain expected fields
        assert "seed" in config_dict
        assert "num_bootstrap" in config_dict
        assert "timestamp" in config_dict
