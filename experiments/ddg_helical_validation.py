#!/usr/bin/env python3
"""
ΔΔG-Pairs + Helical Band Sweep Validation Experiment.

This module implements a comprehensive, pre-registered validation framework for DNA breathing
dynamics encoding, testing sensitivity to thermodynamic perturbations (ΔΔG) through helical
band sweep analysis.

Research Question:
-----------------
Can the DNA breathing dynamics encoding detect thermodynamic perturbations (ΔΔG) from 
single-point mutations across helical frequency bands centered at ~10.5 bp/turn?

Experimental Design:
-------------------
1. ΔΔG Calculator and Pair Generation:
   - Uses SantaLucia 1998 nearest-neighbor thermodynamic parameters
   - Generates all single-point mutants for each wildtype sequence
   - Bins WT-mutant pairs by |ΔΔG| tertiles (low/mid/high)

2. Helical Band Sweep Analysis:
   - Spectral sweep around helical frequency f₀ ≈ 1/10.5 bp/turn
   - Centers: 10.3, 10.4, 10.5, 10.6, 10.7 bp/turn
   - Widths: 1%, 2%, 3%, 5%, 6% of f₀
   - Features: peak_magnitude, band_power, phase_coherence, SNR

3. Statistical Framework:
   - Paired t-test or Wilcoxon signed-rank (based on normality)
   - Cohen's d for paired designs with 95% BCa bootstrap CI (≥1000 reps)
   - Benjamini-Hochberg FDR correction across all tests
   - Jonckheere-Terpstra trend test for dose-response across ΔΔG bins

4. Control Analyses:
   - Off-band windows (±20% from f₀) → expect null effects
   - Label-shuffle permutations (500 reps) → empirical p_null

5. Artifacts Generation:
   - pairs.csv: All 60,000 WT-mutant pairs with ΔΔG values and bin assignments
   - stats.csv: Statistics for all 300 (bin × center × width × feature) conditions
   - trend.csv: Jonckheere-Terpstra trend test results
   - config.json: Full parameter configuration with seeds
   - report.md: Automated summary with acceptance criteria evaluation
   - DETAILED_FINDINGS.md: Comprehensive scientific analysis

Pre-registered Acceptance Criteria:
----------------------------------
- Primary: At least one (center, width, feature) in high-ΔΔG bin achieves
           |d| ≥ 0.5, FDR-q < 0.05, 95% CI excludes 0
- Robustness: Effect replicates across ≥2 adjacent centers or widths
- Specificity: Off-band and shuffle controls remain non-significant (|d| < 0.2, q ≥ 0.1)

Validation Results:
------------------
Executed on 1,000 real Brunello CRISPR sequences (60,000 WT-mutant pairs):
- Primary: FAIL (max |d| = 0.061, far below 0.5 threshold)
- Robustness: FAIL (no primary hits to evaluate)
- Specificity: FAIL (89 shuffle violations)

Key Finding: DNA breathing dynamics encoding shows statistically detectable but practically
negligible sensitivity to single-mutation thermodynamic perturbations. Larger perturbations
or methodological refinements may be needed.

Usage:
------
    # Create real data subsample
    head -2000 data/raw/brunello_parsed.fasta > data/raw/brunello_1k_subsample.fasta
    
    # Run validation
    python experiments/ddg_helical_validation.py \\
        --input data/raw/brunello_1k_subsample.fasta \\
        --output artifacts/ddg_helical_validation \\
        --seed 42
"""

import argparse
import csv
import json
import sys
import warnings
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Add the gist source to path
gist_src = Path(__file__).parent.parent / "gists/breathing/czt_feature_extractor/src"
sys.path.insert(0, str(gist_src))

from dna_breathing_gist import (
    NEAREST_NEIGHBOR_DG,
    czt_analysis,
    encode_sequence,
    extract_features,
)

try:
    from scipy import stats as scipy_stats
except ImportError:
    raise ImportError("scipy is required for statistical tests")


# =============================================================================
# Constants and Configuration
# =============================================================================

HELICAL_PERIOD_DEFAULT = 10.5  # bp/turn for B-form DNA
SEED_DEFAULT = 42
NUM_BOOTSTRAP_DEFAULT = 10000
NUM_PERMUTATIONS_DEFAULT = 1000
NUM_BINS_DEFAULT = 3  # low/mid/high tertiles

# Center sweep values (bp/turn)
CENTER_SWEEP_DEFAULT = [10.3, 10.4, 10.5, 10.6, 10.7]
# Width sweep as fraction of f₀
WIDTH_SWEEP_DEFAULT = [0.01, 0.02, 0.03, 0.05, 0.06]

# Feature names for extraction
FEATURE_NAMES = ["peak_mag", "band_power", "phase_coherence", "snr"]

# Fallback ΔG for unknown dinucleotides (neutral stability)
DG_FALLBACK = -1.5

# Adjacency thresholds for robustness criterion (in bp/turn and fraction)
CENTER_ADJACENCY_THRESHOLD = 0.15  # bp/turn difference for adjacent centers
WIDTH_ADJACENCY_THRESHOLD = 0.02   # Fraction difference for adjacent widths

# Seed hash modulo for generating sequence-specific deterministic seeds
# Limits the offset range while ensuring different sequences get different seeds
SEED_HASH_MODULO = 1000


@dataclass
class ValidationConfig:
    """Configuration for ΔΔG helical validation experiment."""

    seed: int = SEED_DEFAULT
    num_bootstrap: int = NUM_BOOTSTRAP_DEFAULT
    num_permutations: int = NUM_PERMUTATIONS_DEFAULT
    num_bins: int = NUM_BINS_DEFAULT
    center_sweep: List[float] = field(default_factory=lambda: CENTER_SWEEP_DEFAULT.copy())
    width_sweep: List[float] = field(default_factory=lambda: WIDTH_SWEEP_DEFAULT.copy())
    off_band_offset: float = 0.20  # ±20% away from f₀ for control
    alpha: float = 0.05  # Significance level
    effect_size_threshold: float = 0.5  # |d| ≥ 0.5 for primary criterion
    control_effect_threshold: float = 0.2  # |d| < 0.2 for specificity
    control_q_threshold: float = 0.1  # q ≥ 0.1 for specificity
    temperature_k: float = 310.15  # Physiological temperature (37°C)
    at_lifetime: float = 5.0  # AT breathing lifetime (ms)
    gc_lifetime: float = 25.0  # GC breathing lifetime (ms)

    def to_dict(self) -> Dict:
        """Convert config to dictionary for JSON serialization."""
        return {
            "seed": self.seed,
            "num_bootstrap": self.num_bootstrap,
            "num_permutations": self.num_permutations,
            "num_bins": self.num_bins,
            "center_sweep": self.center_sweep,
            "width_sweep": self.width_sweep,
            "off_band_offset": self.off_band_offset,
            "alpha": self.alpha,
            "effect_size_threshold": self.effect_size_threshold,
            "control_effect_threshold": self.control_effect_threshold,
            "control_q_threshold": self.control_q_threshold,
            "temperature_k": self.temperature_k,
            "at_lifetime": self.at_lifetime,
            "gc_lifetime": self.gc_lifetime,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }


# =============================================================================
# ΔΔG Calculator Functions
# =============================================================================


def compute_total_dg(seq: str) -> float:
    """
    Compute total ΔG° for a sequence using nearest-neighbor model.

    Uses SantaLucia 1998 unified parameters from NEAREST_NEIGHBOR_DG.

    Args:
        seq: DNA sequence string (uppercase ATGC only)

    Returns:
        Total ΔG° in kcal/mol
    """
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0

    total_dg = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]
        if dinuc in NEAREST_NEIGHBOR_DG:
            total_dg += NEAREST_NEIGHBOR_DG[dinuc]
        else:
            # Fallback for unknown dinucleotides
            total_dg += DG_FALLBACK

    return total_dg


def compute_ddg(wt_seq: str, mut_seq: str) -> float:
    """
    Compute ΔΔG (mutation-induced free-energy change).

    ΔΔG = ΔG(mutant) - ΔG(wildtype)
    Positive ΔΔG indicates destabilization, negative indicates stabilization.

    Args:
        wt_seq: Wildtype DNA sequence
        mut_seq: Mutant DNA sequence

    Returns:
        ΔΔG in kcal/mol
    """
    dg_wt = compute_total_dg(wt_seq)
    dg_mut = compute_total_dg(mut_seq)
    return dg_mut - dg_wt


def generate_point_mutant(seq: str, position: int, target_base: str) -> str:
    """
    Generate a point mutant at the specified position.

    Args:
        seq: Original DNA sequence
        position: 0-indexed position to mutate
        target_base: Target base (A, T, G, or C)

    Returns:
        Mutated sequence
    """
    seq = seq.upper()
    target_base = target_base.upper()
    if position < 0 or position >= len(seq):
        raise ValueError(f"Position {position} out of range for sequence length {len(seq)}")
    if target_base not in "ATGC":
        raise ValueError(f"Invalid target base: {target_base}")

    return seq[:position] + target_base + seq[position + 1 :]


def generate_mutants_for_sequence(
    seq: str, seed: int = 42
) -> List[Tuple[str, int, str, float]]:
    """
    Generate all possible single-point mutants for a sequence.

    Args:
        seq: Wildtype DNA sequence
        seed: Random seed for reproducibility

    Returns:
        List of (mutant_seq, position, original_base, ddg) tuples
    """
    seq = seq.upper()
    bases = ["A", "T", "G", "C"]
    mutants = []

    for pos in range(len(seq)):
        original = seq[pos]
        if original not in bases:
            continue
        for target in bases:
            if target == original:
                continue
            mut_seq = generate_point_mutant(seq, pos, target)
            ddg = compute_ddg(seq, mut_seq)
            mutants.append((mut_seq, pos, original, ddg))

    return mutants


def select_mutants_by_ddg_magnitude(
    seq: str,
    target_ddg_magnitudes: List[float],
    tolerance: float = 0.1,
    seed: int = 42,
) -> List[Tuple[str, float]]:
    """
    Select mutants that achieve targeted |ΔΔG| magnitudes.

    Args:
        seq: Wildtype DNA sequence
        target_ddg_magnitudes: List of target |ΔΔG| values
        tolerance: Acceptable deviation from target
        seed: Random seed

    Returns:
        List of (mutant_seq, actual_ddg) tuples, one per target magnitude
    """
    np.random.seed(seed)
    all_mutants = generate_mutants_for_sequence(seq, seed)

    selected = []
    for target_mag in target_ddg_magnitudes:
        # Find mutants within tolerance of target magnitude
        candidates = [
            (mut, ddg)
            for mut, pos, orig, ddg in all_mutants
            if abs(abs(ddg) - target_mag) <= tolerance
        ]
        if candidates:
            # Select randomly among candidates using numpy
            idx = np.random.randint(len(candidates))
            choice = candidates[idx]
            selected.append(choice)
        else:
            # If no exact match, find closest
            sorted_by_distance = sorted(
                all_mutants, key=lambda x: abs(abs(x[3]) - target_mag)
            )
            if sorted_by_distance:
                best = sorted_by_distance[0]
                selected.append((best[0], best[3]))

    return selected


def generate_wt_mutant_pairs(
    sequences: List[Tuple[str, str]], num_bins: int = 3, seed: int = 42
) -> Tuple[List[Dict], List[float]]:
    """
    Generate WT-mutant pairs and bin them by |ΔΔG| tertiles.

    For each WT sequence, generates mutants with varying ΔΔG magnitudes
    and assigns pairs to bins based on |ΔΔG| tertiles.

    Args:
        sequences: List of (seq_id, sequence) tuples
        num_bins: Number of bins (default: 3 for low/mid/high)
        seed: Random seed

    Returns:
        Tuple of (pairs_list, bin_edges) where pairs_list contains dicts
        with keys: wt_id, mut_id, wt_seq, mut_seq, delta_delta_g, bin
    """
    np.random.seed(seed)

    all_pairs = []
    pair_counter = 0

    for seq_id, seq in sequences:
        # Generate all possible mutants for this sequence
        mutants = generate_mutants_for_sequence(seq, seed + hash(seq_id) % SEED_HASH_MODULO)

        for mut_seq, pos, orig, ddg in mutants:
            pair_counter += 1
            all_pairs.append(
                {
                    "wt_id": seq_id,
                    "mut_id": f"{seq_id}_mut_{pair_counter}",
                    "wt_seq": seq,
                    "mut_seq": mut_seq,
                    "delta_delta_g": ddg,
                    "abs_ddg": abs(ddg),
                    "position": pos,
                    "original_base": orig,
                }
            )

    if not all_pairs:
        return [], []

    # Compute tertile edges for binning
    abs_ddg_values = [p["abs_ddg"] for p in all_pairs]
    bin_edges = np.percentile(abs_ddg_values, np.linspace(0, 100, num_bins + 1))

    # Assign bins (0=low, 1=mid, 2=high)
    bin_labels = ["low", "mid", "high"][:num_bins]
    for pair in all_pairs:
        abs_ddg = pair["abs_ddg"]
        bin_idx = np.searchsorted(bin_edges[1:-1], abs_ddg, side="right")
        pair["bin"] = bin_labels[bin_idx]
        pair["bin_idx"] = bin_idx

    return all_pairs, bin_edges.tolist()


# =============================================================================
# Band Sweep Analysis Functions
# =============================================================================


def extract_features_at_band(
    seq: str,
    center_period: float,
    width_fraction: float,
    config: ValidationConfig,
) -> Dict[str, float]:
    """
    Extract spectral features at a specific band configuration.

    Args:
        seq: DNA sequence
        center_period: Center of band in bp/turn (e.g., 10.5)
        width_fraction: Width as fraction of center frequency (e.g., 0.02)
        config: Validation configuration

    Returns:
        Dictionary of feature values
    """
    f0 = 1.0 / center_period
    band_width = f0 * width_fraction

    # Encode sequence
    signal = encode_sequence(
        seq,
        at_lifetime=config.at_lifetime,
        gc_lifetime=config.gc_lifetime,
        helical_period=center_period,
        apply_helical=True,
    )

    # CZT analysis
    freqs, spectrum = czt_analysis(signal, f0=f0, band_width=band_width)

    # Extract features
    features = extract_features(freqs, spectrum, f0, band_width)

    # Map to our feature names
    return {
        "peak_mag": features["peak_mag"],
        "band_power": features["band_energy"],
        "phase_coherence": features["phase_coherence"],
        "snr": features["snr"],
    }


def run_band_sweep(
    pairs: List[Dict],
    config: ValidationConfig,
) -> List[Dict]:
    """
    Run band sweep analysis across all pairs and parameter combinations.

    Args:
        pairs: List of WT-mutant pair dictionaries
        config: Validation configuration

    Returns:
        List of results for each (pair, center, width, feature) combination
    """
    results = []

    for pair in pairs:
        wt_seq = pair["wt_seq"]
        mut_seq = pair["mut_seq"]

        for center in config.center_sweep:
            for width in config.width_sweep:
                try:
                    wt_features = extract_features_at_band(wt_seq, center, width, config)
                    mut_features = extract_features_at_band(mut_seq, center, width, config)

                    for feature_name in FEATURE_NAMES:
                        wt_val = wt_features[feature_name]
                        mut_val = mut_features[feature_name]
                        diff = mut_val - wt_val

                        results.append(
                            {
                                "pair_id": pair["mut_id"],
                                "wt_id": pair["wt_id"],
                                "delta_delta_g": pair["delta_delta_g"],
                                "bin": pair["bin"],
                                "bin_idx": pair["bin_idx"],
                                "center": center,
                                "width": width,
                                "feature": feature_name,
                                "wt_value": wt_val,
                                "mut_value": mut_val,
                                "diff": diff,
                            }
                        )
                except Exception as e:
                    warnings.warn(f"Error processing pair {pair['mut_id']}: {e}")
                    continue

    return results


def run_off_band_control(
    pairs: List[Dict],
    config: ValidationConfig,
) -> List[Dict]:
    """
    Run off-band control analysis at ±20% away from f₀.

    Args:
        pairs: List of WT-mutant pair dictionaries
        config: Validation configuration

    Returns:
        List of off-band control results
    """
    results = []
    f0_center = HELICAL_PERIOD_DEFAULT

    # Off-band centers: ±20% away
    off_band_low = f0_center * (1 - config.off_band_offset)
    off_band_high = f0_center * (1 + config.off_band_offset)

    for pair in pairs:
        wt_seq = pair["wt_seq"]
        mut_seq = pair["mut_seq"]

        for off_center in [off_band_low, off_band_high]:
            for width in [0.03]:  # Use median width for control
                try:
                    wt_features = extract_features_at_band(wt_seq, off_center, width, config)
                    mut_features = extract_features_at_band(mut_seq, off_center, width, config)

                    for feature_name in FEATURE_NAMES:
                        wt_val = wt_features[feature_name]
                        mut_val = mut_features[feature_name]
                        diff = mut_val - wt_val

                        results.append(
                            {
                                "pair_id": pair["mut_id"],
                                "wt_id": pair["wt_id"],
                                "delta_delta_g": pair["delta_delta_g"],
                                "bin": pair["bin"],
                                "center": off_center,
                                "width": width,
                                "feature": feature_name,
                                "wt_value": wt_val,
                                "mut_value": mut_val,
                                "diff": diff,
                                "control_type": "off_band",
                            }
                        )
                except Exception:
                    continue

    return results


# =============================================================================
# Statistical Framework
# =============================================================================


def bca_bootstrap_ci(
    data: np.ndarray,
    stat_func,
    num_bootstrap: int = 10000,
    alpha: float = 0.05,
    seed: int = 42,
) -> Tuple[float, float, float]:
    """
    Compute BCa (Bias-Corrected and Accelerated) bootstrap confidence interval.

    Args:
        data: 1D array of data points
        stat_func: Function to compute statistic (takes array, returns scalar)
        num_bootstrap: Number of bootstrap iterations
        alpha: Significance level (e.g., 0.05 for 95% CI)
        seed: Random seed

    Returns:
        Tuple of (point_estimate, ci_low, ci_high)
    """
    np.random.seed(seed)
    n = len(data)

    if n < 2:
        return np.nan, np.nan, np.nan

    # Point estimate
    theta_hat = stat_func(data)

    # Bootstrap distribution
    boot_thetas = []
    for _ in range(num_bootstrap):
        boot_sample = np.random.choice(data, size=n, replace=True)
        boot_thetas.append(stat_func(boot_sample))
    boot_thetas = np.array(boot_thetas)

    # Handle NaN values
    boot_thetas = boot_thetas[~np.isnan(boot_thetas)]
    if len(boot_thetas) < 10:
        return theta_hat, np.nan, np.nan

    # Bias correction (z0)
    # Handle edge cases where all bootstrap values are on one side of theta_hat
    prop_below = np.mean(boot_thetas < theta_hat)
    # Clamp proportion to avoid -inf/+inf from norm.ppf
    prop_below = np.clip(prop_below, 0.001, 0.999)
    z0 = scipy_stats.norm.ppf(prop_below)

    # Acceleration (a) using jackknife
    jack_thetas = []
    for i in range(n):
        jack_sample = np.concatenate([data[:i], data[i + 1 :]])
        if len(jack_sample) > 0:
            jack_thetas.append(stat_func(jack_sample))
    jack_thetas = np.array(jack_thetas)
    jack_mean = np.mean(jack_thetas)
    jack_diffs = jack_mean - jack_thetas

    numerator = np.sum(jack_diffs**3)
    denominator = 6 * (np.sum(jack_diffs**2)) ** 1.5

    if denominator == 0:
        a = 0.0
    else:
        a = numerator / denominator

    # BCa percentiles
    z_alpha_low = scipy_stats.norm.ppf(alpha / 2)
    z_alpha_high = scipy_stats.norm.ppf(1 - alpha / 2)

    def bca_percentile(z_alpha: float) -> float:
        num = z0 + z_alpha
        denom = 1 - a * (z0 + z_alpha)
        if denom == 0:
            return 0.5
        adjusted = z0 + num / denom
        return scipy_stats.norm.cdf(adjusted)

    p_low = bca_percentile(z_alpha_low)
    p_high = bca_percentile(z_alpha_high)

    # Clamp percentiles
    p_low = max(0.001, min(0.999, p_low))
    p_high = max(0.001, min(0.999, p_high))

    ci_low = np.percentile(boot_thetas, p_low * 100)
    ci_high = np.percentile(boot_thetas, p_high * 100)

    return theta_hat, ci_low, ci_high


def cohens_d_paired(diffs: np.ndarray) -> float:
    """
    Compute Cohen's d for paired designs (dz).

    d = mean(diff) / std(diff)

    Args:
        diffs: Array of paired differences

    Returns:
        Cohen's d effect size
    """
    if len(diffs) < 2:
        return np.nan
    mean_diff = np.mean(diffs)
    std_diff = np.std(diffs, ddof=1)
    if std_diff == 0:
        return np.nan
    return mean_diff / std_diff


def compute_paired_statistics(
    diffs: np.ndarray,
    num_bootstrap: int = 10000,
    num_perm: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> Dict:
    """
    Compute comprehensive paired statistics for a set of differences.

    Args:
        diffs: Array of paired differences
        num_bootstrap: Number of bootstrap samples for CI
        num_perm: Number of permutations for null distribution
        alpha: Significance level
        seed: Random seed

    Returns:
        Dictionary with n, mean_diff, d, ci_low, ci_high, p_value, test_type
    """
    np.random.seed(seed)
    n = len(diffs)

    if n < 2:
        return {
            "n": n,
            "mean_diff": np.nan,
            "d": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p_value": np.nan,
            "test_type": "insufficient_data",
        }

    mean_diff = np.mean(diffs)

    # Cohen's d for paired data
    d = cohens_d_paired(diffs)

    # BCa bootstrap CI for Cohen's d
    _, ci_low, ci_high = bca_bootstrap_ci(
        diffs, cohens_d_paired, num_bootstrap=num_bootstrap, alpha=alpha, seed=seed
    )

    # Test for normality (Shapiro-Wilk if n <= 50, else D'Agostino-Pearson)
    if n <= 50:
        _, normality_p = scipy_stats.shapiro(diffs)
    else:
        _, normality_p = scipy_stats.normaltest(diffs)

    # Choose test based on normality
    if normality_p > 0.05:
        # Use paired t-test
        t_stat, p_value = scipy_stats.ttest_1samp(diffs, 0)
        test_type = "paired_t"
    else:
        # Use Wilcoxon signed-rank test
        try:
            stat, p_value = scipy_stats.wilcoxon(diffs, alternative="two-sided")
            test_type = "wilcoxon"
        except ValueError:
            # Fallback if Wilcoxon fails (e.g., all zeros)
            p_value = 1.0
            test_type = "wilcoxon_failed"

    return {
        "n": n,
        "mean_diff": mean_diff,
        "d": d,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_value": p_value,
        "test_type": test_type,
    }


def benjamini_hochberg_correction(p_values: List[float]) -> List[float]:
    """
    Apply Benjamini-Hochberg FDR correction to p-values.

    Args:
        p_values: List of raw p-values

    Returns:
        List of FDR-corrected q-values
    """
    n = len(p_values)
    if n == 0:
        return []

    # Handle NaN values
    p_array = np.array(p_values)
    valid_mask = ~np.isnan(p_array)
    valid_p = p_array[valid_mask]

    if len(valid_p) == 0:
        return [np.nan] * n

    # Sort and compute q-values
    sorted_indices = np.argsort(valid_p)
    sorted_p = valid_p[sorted_indices]

    q_values = np.zeros_like(sorted_p)
    n_valid = len(sorted_p)

    # BH procedure
    for i in range(n_valid - 1, -1, -1):
        rank = i + 1
        if i == n_valid - 1:
            q_values[i] = sorted_p[i]
        else:
            q_values[i] = min(q_values[i + 1], sorted_p[i] * n_valid / rank)

    # Unsort
    unsorted_q = np.zeros(n_valid)
    unsorted_q[sorted_indices] = q_values

    # Reconstruct with NaNs
    result = np.full(n, np.nan)
    result[valid_mask] = unsorted_q

    return list(result)


def label_shuffle_permutation(
    diffs: np.ndarray, num_perm: int = 1000, seed: int = 42
) -> float:
    """
    Compute empirical p-value via label shuffling within pairs.

    Args:
        diffs: Array of paired differences
        num_perm: Number of permutations
        seed: Random seed

    Returns:
        Empirical p-value
    """
    np.random.seed(seed)
    observed_stat = abs(np.mean(diffs))

    count = 0
    for _ in range(num_perm):
        # Randomly flip signs (equivalent to swapping labels within pairs)
        signs = np.random.choice([-1, 1], size=len(diffs))
        perm_diffs = diffs * signs
        perm_stat = abs(np.mean(perm_diffs))
        if perm_stat >= observed_stat:
            count += 1

    # Add-one smoothing to avoid p=0
    return (count + 1) / (num_perm + 1)


def jonckheere_terpstra_test(
    groups: List[np.ndarray], alternative: str = "increasing"
) -> Tuple[float, float]:
    """
    Perform Jonckheere-Terpstra trend test for ordered groups.

    Tests if there is a monotonic trend across ordered groups.

    Args:
        groups: List of arrays, one per ordered group
        alternative: "increasing" or "decreasing"

    Returns:
        Tuple of (J statistic, p-value)
    """
    k = len(groups)
    if k < 2:
        return np.nan, np.nan

    # Compute J statistic (number of concordant pairs across groups)
    J = 0.0
    for i in range(k - 1):
        for j in range(i + 1, k):
            for xi in groups[i]:
                for xj in groups[j]:
                    if xj > xi:
                        J += 1
                    elif xj == xi:
                        J += 0.5

    # Compute expected value and variance under null
    n = [len(g) for g in groups]
    N = sum(n)

    E_J = (N**2 - sum(ni**2 for ni in n)) / 4

    # Variance computation (assuming no ties for simplicity)
    term1 = N**2 * (2 * N + 3) - sum(ni**2 * (2 * ni + 3) for ni in n)
    var_J = term1 / 72

    if var_J <= 0:
        return J, np.nan

    # Z statistic
    z = (J - E_J) / np.sqrt(var_J)

    # Two-sided p-value
    if alternative == "increasing":
        p_value = 1 - scipy_stats.norm.cdf(z)
    elif alternative == "decreasing":
        p_value = scipy_stats.norm.cdf(z)
    else:
        p_value = 2 * (1 - scipy_stats.norm.cdf(abs(z)))

    return J, p_value


# =============================================================================
# Results Aggregation
# =============================================================================


def aggregate_by_condition(
    sweep_results: List[Dict],
    config: ValidationConfig,
) -> List[Dict]:
    """
    Aggregate sweep results by (bin, center, width, feature) and compute statistics.

    Args:
        sweep_results: Raw results from band sweep
        config: Validation configuration

    Returns:
        List of aggregated statistics per condition
    """
    # Group by condition
    grouped = defaultdict(list)

    for result in sweep_results:
        key = (result["bin"], result["center"], result["width"], result["feature"])
        grouped[key].append(result["diff"])

    stats_results = []
    all_p_values = []
    stats_with_keys = []

    for key, diffs in grouped.items():
        bin_name, center, width, feature = key
        diffs_array = np.array(diffs)

        stats = compute_paired_statistics(
            diffs_array,
            num_bootstrap=config.num_bootstrap,
            num_perm=config.num_permutations,
            alpha=config.alpha,
            seed=config.seed,
        )

        result = {
            "bin": bin_name,
            "center": center,
            "width": width,
            "feature": feature,
            "n": stats["n"],
            "mean_diff": stats["mean_diff"],
            "d": stats["d"],
            "ci_low": stats["ci_low"],
            "ci_high": stats["ci_high"],
            "p": stats["p_value"],
            "test_type": stats["test_type"],
        }

        all_p_values.append(stats["p_value"])
        stats_with_keys.append(result)

    # Apply BH-FDR correction
    q_values = benjamini_hochberg_correction(all_p_values)

    for i, result in enumerate(stats_with_keys):
        result["q"] = q_values[i]
        stats_results.append(result)

    return stats_results


def compute_trend_statistics(
    stats_results: List[Dict],
    config: ValidationConfig,
) -> List[Dict]:
    """
    Compute bin-level summaries and trend test results.

    Tests if |d| increases with |ΔΔG| bin using Jonckheere-Terpstra.

    Args:
        stats_results: Aggregated statistics from aggregate_by_condition
        config: Validation configuration

    Returns:
        List of trend test results per (center, width, feature) combination
    """
    # Group by (center, width, feature) across bins
    grouped = defaultdict(lambda: defaultdict(list))

    for result in stats_results:
        key = (result["center"], result["width"], result["feature"])
        bin_name = result["bin"]
        d = result["d"]
        if not np.isnan(d):
            grouped[key][bin_name].append(abs(d))

    trend_results = []
    bin_order = ["low", "mid", "high"]

    for key, bin_data in grouped.items():
        center, width, feature = key

        # Prepare groups in order
        groups = []
        bin_means = {}
        for bin_name in bin_order:
            if bin_name in bin_data:
                groups.append(np.array(bin_data[bin_name]))
                bin_means[bin_name] = np.mean(bin_data[bin_name])
            else:
                groups.append(np.array([]))
                bin_means[bin_name] = np.nan

        # Skip if not enough data
        valid_groups = [g for g in groups if len(g) > 0]
        if len(valid_groups) < 2:
            continue

        # Jonckheere-Terpstra test for increasing trend
        j_stat, jt_p = jonckheere_terpstra_test(valid_groups, alternative="increasing")

        trend_results.append(
            {
                "center": center,
                "width": width,
                "feature": feature,
                "low_mean_abs_d": bin_means.get("low", np.nan),
                "mid_mean_abs_d": bin_means.get("mid", np.nan),
                "high_mean_abs_d": bin_means.get("high", np.nan),
                "jt_statistic": j_stat,
                "jt_p_value": jt_p,
            }
        )

    return trend_results


# =============================================================================
# Acceptance Criteria Evaluation
# =============================================================================


def evaluate_acceptance_criteria(
    stats_results: List[Dict],
    control_results: List[Dict],
    shuffle_p_values: Dict[str, float],
    config: ValidationConfig,
) -> Dict:
    """
    Evaluate pre-registered acceptance criteria.

    Primary: At least one (center, width, feature) in ΔΔG-high bin achieves
             |d| ≥ 0.5, FDR-q < 0.05, 95% CI excludes 0
    Robustness: Effect replicates across ≥2 adjacent centers or widths
    Specificity: Off-band and shuffle controls remain non-significant

    Args:
        stats_results: Main analysis statistics
        control_results: Off-band control statistics
        shuffle_p_values: Label-shuffle permutation p-values
        config: Validation configuration

    Returns:
        Dictionary with pass/fail for each criterion
    """
    # Primary criterion: high bin significant effects
    high_bin_results = [r for r in stats_results if r["bin"] == "high"]
    primary_hits = []

    for r in high_bin_results:
        d = r.get("d", np.nan)
        q = r.get("q", np.nan)
        ci_low = r.get("ci_low", np.nan)
        ci_high = r.get("ci_high", np.nan)

        if np.isnan(d) or np.isnan(q) or np.isnan(ci_low) or np.isnan(ci_high):
            continue

        # Check criteria
        if abs(d) >= config.effect_size_threshold and q < config.alpha:
            # CI excludes 0
            if ci_low > 0 or ci_high < 0:
                primary_hits.append(r)

    primary_pass = len(primary_hits) >= 1

    # Robustness: Check for replication across adjacent parameters
    robustness_pass = False
    if primary_hits:
        # Group hits by feature
        hits_by_feature = defaultdict(list)
        for h in primary_hits:
            hits_by_feature[h["feature"]].append((h["center"], h["width"]))

        for feature, params in hits_by_feature.items():
            if len(params) >= 2:
                centers = sorted(set(p[0] for p in params))
                widths = sorted(set(p[1] for p in params))

                # Check for adjacent centers
                for i in range(len(centers) - 1):
                    if abs(centers[i + 1] - centers[i]) <= CENTER_ADJACENCY_THRESHOLD:
                        robustness_pass = True
                        break

                # Check for adjacent widths
                if not robustness_pass:
                    for i in range(len(widths) - 1):
                        if abs(widths[i + 1] - widths[i]) <= WIDTH_ADJACENCY_THRESHOLD:
                            robustness_pass = True
                            break

    # Specificity: Off-band controls should be non-significant
    off_band_violations = 0
    for r in control_results:
        d = r.get("d", np.nan)
        q = r.get("q", np.nan)
        if not np.isnan(d) and not np.isnan(q):
            if abs(d) >= config.control_effect_threshold and q < config.control_q_threshold:
                off_band_violations += 1

    # Shuffle controls should be non-significant
    shuffle_violations = 0
    for key, p in shuffle_p_values.items():
        if p < config.alpha:
            shuffle_violations += 1

    specificity_pass = off_band_violations == 0 and shuffle_violations == 0

    # Overall pass
    overall_pass = primary_pass and robustness_pass and specificity_pass

    return {
        "primary_criterion_pass": primary_pass,
        "primary_hits_count": len(primary_hits),
        "robustness_criterion_pass": robustness_pass,
        "specificity_criterion_pass": specificity_pass,
        "off_band_violations": off_band_violations,
        "shuffle_violations": shuffle_violations,
        "overall_pass": overall_pass,
    }


# =============================================================================
# Artifacts Generation
# =============================================================================


def save_artifacts(
    output_dir: Path,
    config: ValidationConfig,
    pairs: List[Dict],
    bin_edges: List[float],
    stats_results: List[Dict],
    control_stats: List[Dict],
    trend_results: List[Dict],
    shuffle_p_values: Dict[str, float],
    acceptance: Dict,
) -> None:
    """
    Save all artifacts to the output directory.

    Creates:
    - config.json: Full sweep and analysis params + seeds
    - pairs.csv: WT_id, mut_id, ΔΔG, bin
    - stats.csv: Per (bin, center, width, feature) statistics
    - trend.csv: Bin-level summaries and trend test results
    - report.md: One-page summary with pass/fail

    Args:
        output_dir: Output directory path
        config: Validation configuration
        pairs: WT-mutant pairs
        bin_edges: Bin edge values
        stats_results: Aggregated statistics
        control_stats: Control analysis statistics
        trend_results: Trend test results
        shuffle_p_values: Label-shuffle p-values
        acceptance: Acceptance criteria evaluation
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. config.json
    config_data = config.to_dict()
    config_data["bin_edges"] = bin_edges
    config_data["n_pairs"] = len(pairs)

    with open(output_dir / "config.json", "w") as f:
        json.dump(config_data, f, indent=2)

    # 2. pairs.csv
    with open(output_dir / "pairs.csv", "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["wt_id", "mut_id", "delta_delta_g", "bin"]
        )
        writer.writeheader()
        for pair in pairs:
            writer.writerow(
                {
                    "wt_id": pair["wt_id"],
                    "mut_id": pair["mut_id"],
                    "delta_delta_g": pair["delta_delta_g"],
                    "bin": pair["bin"],
                }
            )

    # 3. stats.csv
    stats_fieldnames = ["bin", "center", "width", "feature", "n", "mean_diff", "d", "ci_low", "ci_high", "p", "q", "test_type"]
    with open(output_dir / "stats.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=stats_fieldnames)
        writer.writeheader()
        for stat in stats_results:
            writer.writerow({k: stat.get(k, "") for k in stats_fieldnames})

    # 4. trend.csv
    trend_fieldnames = ["center", "width", "feature", "low_mean_abs_d", "mid_mean_abs_d", "high_mean_abs_d", "jt_statistic", "jt_p_value"]
    with open(output_dir / "trend.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=trend_fieldnames)
        writer.writeheader()
        for trend in trend_results:
            writer.writerow({k: trend.get(k, "") for k in trend_fieldnames})

    # 5. report.md
    report = generate_report(
        config, pairs, bin_edges, stats_results, control_stats,
        trend_results, shuffle_p_values, acceptance
    )
    with open(output_dir / "report.md", "w") as f:
        f.write(report)


def generate_report(
    config: ValidationConfig,
    pairs: List[Dict],
    bin_edges: List[float],
    stats_results: List[Dict],
    control_stats: List[Dict],
    trend_results: List[Dict],
    shuffle_p_values: Dict[str, float],
    acceptance: Dict,
) -> str:
    """Generate markdown report summarizing validation results."""
    status = "PASS ✓" if acceptance["overall_pass"] else "FAIL ✗"

    report = f"""# ΔΔG-Pairs + Helical Band Sweep Validation Report

**Status: {status}**

Generated: {datetime.now(timezone.utc).isoformat()}

---

## Summary

This report evaluates the DNA breathing dynamics encoding's sensitivity to
thermodynamic perturbations (ΔΔG) across helical frequency bands.

### Configuration

- **Seed**: {config.seed}
- **Bootstrap samples**: {config.num_bootstrap}
- **Permutations**: {config.num_permutations}
- **Center sweep**: {config.center_sweep}
- **Width sweep**: {config.width_sweep}
- **Temperature**: {config.temperature_k} K

### Data

- **Total pairs**: {len(pairs)}
- **Bin edges (|ΔΔG| kcal/mol)**: {[round(e, 3) for e in bin_edges]}

---

## Acceptance Criteria Evaluation

| Criterion | Result |
|-----------|--------|
| Primary (|d| ≥ 0.5, q < 0.05, CI excl. 0 in high bin) | {"PASS ✓" if acceptance["primary_criterion_pass"] else "FAIL ✗"} ({acceptance["primary_hits_count"]} hits) |
| Robustness (replicates across ≥2 adjacent params) | {"PASS ✓" if acceptance["robustness_criterion_pass"] else "FAIL ✗"} |
| Specificity (off-band/shuffle non-significant) | {"PASS ✓" if acceptance["specificity_criterion_pass"] else "FAIL ✗"} |

### Control Violations

- Off-band violations: {acceptance["off_band_violations"]}
- Shuffle violations: {acceptance["shuffle_violations"]}

---

## Top Results (High Bin, Sorted by |d|)

| Center | Width | Feature | n | d | 95% CI | q |
|--------|-------|---------|---|---|--------|---|
"""

    # Add top 10 high-bin results sorted by |d|
    high_results = [r for r in stats_results if r["bin"] == "high" and not np.isnan(r.get("d", np.nan))]
    high_results_sorted = sorted(high_results, key=lambda x: abs(x.get("d", 0)), reverse=True)[:10]

    for r in high_results_sorted:
        ci_str = f"[{r.get('ci_low', np.nan):.3f}, {r.get('ci_high', np.nan):.3f}]"
        report += f"| {r['center']} | {r['width']} | {r['feature']} | {r['n']} | {r.get('d', np.nan):.3f} | {ci_str} | {r.get('q', np.nan):.4f} |\n"

    report += """
---

## Trend Test Results (Jonckheere-Terpstra)

| Center | Width | Feature | Low |d̄| | Mid |d̄| | High |d̄| | JT p-value |
|--------|-------|---------|---------|---------|----------|------------|
"""

    for t in trend_results[:15]:
        report += f"| {t['center']} | {t['width']} | {t['feature']} | {t.get('low_mean_abs_d', np.nan):.3f} | {t.get('mid_mean_abs_d', np.nan):.3f} | {t.get('high_mean_abs_d', np.nan):.3f} | {t.get('jt_p_value', np.nan):.4f} |\n"

    report += """
---

## Pre-registration Notes

This validation follows pre-registered acceptance criteria:

1. **Primary**: At least one (center, width, feature) combination in the
   high-ΔΔG bin must achieve |d| ≥ 0.5, FDR-corrected q < 0.05, and
   95% BCa confidence interval excluding zero.

2. **Robustness**: The effect must replicate across at least 2 adjacent
   center frequencies or band widths.

3. **Specificity**: Off-band windows (±20% from f₀) and label-shuffle
   permutation controls must remain non-significant (|d| < 0.2; q ≥ 0.1).

---

*End of report*
"""

    return report


# =============================================================================
# Sequence Loading
# =============================================================================


def load_sequences_from_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    """
    Load sequences from FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        List of (seq_id, sequence) tuples
    """
    sequences = []
    current_id = None
    current_seq = ""

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    sequences.append((current_id, current_seq.upper()))
                current_id = line[1:].split()[0]  # Take first word as ID
                current_seq = ""
            else:
                current_seq += line

    if current_id and current_seq:
        sequences.append((current_id, current_seq.upper()))

    return sequences


# =============================================================================
# Main Entry Point
# =============================================================================


def run_validation(
    input_path: Path,
    output_dir: Path,
    config: Optional[ValidationConfig] = None,
) -> Dict:
    """
    Run the complete ΔΔG helical validation experiment.

    Args:
        input_path: Path to input FASTA file
        output_dir: Path to output directory for artifacts
        config: Validation configuration (uses defaults if None)

    Returns:
        Dictionary with acceptance criteria evaluation results
    """
    if config is None:
        config = ValidationConfig()

    # Set seed (using numpy for consistency)
    np.random.seed(config.seed)

    print(f"Loading sequences from {input_path}...")
    sequences = load_sequences_from_fasta(input_path)
    print(f"Loaded {len(sequences)} sequences")

    if len(sequences) < 2:
        raise ValueError("Need at least 2 sequences for validation")

    # Generate WT-mutant pairs
    print("Generating WT-mutant pairs...")
    pairs, bin_edges = generate_wt_mutant_pairs(
        sequences, num_bins=config.num_bins, seed=config.seed
    )
    print(f"Generated {len(pairs)} pairs")
    print(f"Bin edges: {bin_edges}")

    # Run main band sweep
    print("Running band sweep analysis...")
    sweep_results = run_band_sweep(pairs, config)
    print(f"Collected {len(sweep_results)} measurements")

    # Run off-band control
    print("Running off-band control analysis...")
    off_band_results = run_off_band_control(pairs, config)

    # Aggregate main statistics
    print("Aggregating statistics...")
    stats_results = aggregate_by_condition(sweep_results, config)
    control_stats = aggregate_by_condition(off_band_results, config)

    # Compute trend statistics
    print("Computing trend statistics...")
    trend_results = compute_trend_statistics(stats_results, config)

    # Run label-shuffle permutation test on high-bin results
    print("Running label-shuffle permutation tests...")
    shuffle_p_values = {}

    high_bin_diffs = defaultdict(list)
    for result in sweep_results:
        if result["bin"] == "high":
            key = (result["center"], result["width"], result["feature"])
            high_bin_diffs[key].append(result["diff"])

    for key, diffs in high_bin_diffs.items():
        diffs_array = np.array(diffs)
        p_shuffle = label_shuffle_permutation(
            diffs_array, num_perm=config.num_permutations, seed=config.seed
        )
        shuffle_p_values[str(key)] = p_shuffle

    # Evaluate acceptance criteria
    print("Evaluating acceptance criteria...")
    acceptance = evaluate_acceptance_criteria(
        stats_results, control_stats, shuffle_p_values, config
    )

    # Save artifacts
    print(f"Saving artifacts to {output_dir}...")
    save_artifacts(
        output_dir,
        config,
        pairs,
        bin_edges,
        stats_results,
        control_stats,
        trend_results,
        shuffle_p_values,
        acceptance,
    )

    # Print summary
    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)
    print(f"Overall: {'PASS ✓' if acceptance['overall_pass'] else 'FAIL ✗'}")
    print(f"  Primary criterion: {'PASS' if acceptance['primary_criterion_pass'] else 'FAIL'}")
    print(f"  Robustness: {'PASS' if acceptance['robustness_criterion_pass'] else 'FAIL'}")
    print(f"  Specificity: {'PASS' if acceptance['specificity_criterion_pass'] else 'FAIL'}")
    print(f"\nArtifacts saved to: {output_dir}")

    return acceptance


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="ΔΔG-Pairs + Helical Band Sweep Validation Experiment"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).parent.parent / "data/human/depmap_subsample/sequences.fasta",
        help="Input FASTA file with DNA sequences",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).parent.parent / "artifacts/ddg_helical_validation",
        help="Output directory for artifacts",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=SEED_DEFAULT,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--num-bootstrap",
        type=int,
        default=NUM_BOOTSTRAP_DEFAULT,
        help="Number of bootstrap samples",
    )
    parser.add_argument(
        "--num-permutations",
        type=int,
        default=NUM_PERMUTATIONS_DEFAULT,
        help="Number of permutation samples",
    )
    parser.add_argument(
        "--num-bins",
        type=int,
        default=NUM_BINS_DEFAULT,
        help="Number of ΔΔG bins",
    )

    args = parser.parse_args()

    config = ValidationConfig(
        seed=args.seed,
        num_bootstrap=args.num_bootstrap,
        num_permutations=args.num_permutations,
        num_bins=args.num_bins,
    )

    run_validation(args.input, args.output, config)


if __name__ == "__main__":
    main()
