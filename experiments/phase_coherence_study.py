#!/usr/bin/env python3
"""
Phase-Coherent vs Random-Phase CZT Spectra Comparison Study.

This script implements a rigorous statistical comparison between phase-coherent
and random-phase CZT spectral features for DNA breathing dynamics analysis.

Study Design:
- Condition A (phase-coherent): complex CZT features as computed (magnitude + phase)
- Condition B (random-phase control): identical magnitudes with uniform random phases

The study uses stratified K-fold cross-validation with repeated splits, logistic
regression classification, and comprehensive statistical analysis including
paired t-tests, Wilcoxon signed-rank tests, Cohen's d effect sizes with bootstrap
confidence intervals, and Benjamini-Hochberg FDR correction.

Usage:
    python experiments/phase_coherence_study.py --input data.fasta --output-dir results/
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from scipy.signal import czt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler

# Optional matplotlib import for plotting
try:
    import matplotlib.pyplot as plt

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# =============================================================================
# Configuration: Default Seeds (Fixed for Reproducibility)
# =============================================================================
DEFAULT_GLOBAL_SEED = 137
DEFAULT_PHASE_SEED = 271828
DEFAULT_CV_SEED = 161803

# =============================================================================
# Biophysics Parameters (SantaLucia 1998 Unified)
# =============================================================================
NEAREST_NEIGHBOR_DG: Dict[str, float] = {
    "AA": -1.00,
    "TT": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "TG": -1.45,
    "GT": -1.44,
    "AC": -1.44,
    "CT": -1.28,
    "AG": -1.28,
    "GA": -1.30,
    "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
    "CC": -1.84,
}

# Breathing lifetimes [ms]: AT pairs more flexible, GC pairs more stable
BREATH_LIFETIMES_MS: Dict[str, float] = {
    "A": 5.0,
    "T": 5.0,
    "G": 25.0,
    "C": 25.0,
}

# Dynamically derive ΔG bounds for normalization
_DG_MIN = min(NEAREST_NEIGHBOR_DG.values())
_DG_MAX = max(NEAREST_NEIGHBOR_DG.values())

# Helical period of B-form DNA
DEFAULT_HELICAL_PERIOD = 10.5


# =============================================================================
# Data Classes
# =============================================================================
@dataclass
class StudyConfig:
    """Configuration for the phase coherence study."""

    # Seeds
    global_seed: int = DEFAULT_GLOBAL_SEED
    phase_seed: int = DEFAULT_PHASE_SEED
    cv_seed: int = DEFAULT_CV_SEED

    # Cross-validation
    n_folds: int = 5
    n_repeats: int = 5

    # CZT band parameters
    band_center: float = 1 / DEFAULT_HELICAL_PERIOD  # ~0.0952 cycles/bp
    band_width: float = 0.01  # Half-width in cycles/bp
    czt_points: int = 256

    # Bootstrap/permutation parameters
    n_bootstrap: int = 10000
    n_permutations: int = 1000

    # FDR correction
    fdr_alpha: float = 0.05

    # Power analysis
    expected_delta: float = 0.03
    power_alpha: float = 0.05
    desired_power: float = 0.8

    # Pilot study fraction
    pilot_fraction: Optional[float] = None

    # Plotting
    skip_plots: bool = False


@dataclass
class CZTFeatures:
    """
    Features extracted from CZT spectrum.

    Note: rayleigh_p is stored for diagnostic purposes but excluded from
    the ML feature array (to_array()) since p-values are not suitable as
    features for classification.
    """

    peak_magnitude: float
    peak_freq_idx: int
    phase_at_peak: float
    band_power: float
    spectral_centroid: float
    phase_kappa: float  # von Mises concentration
    rayleigh_z: float
    rayleigh_p: float  # Diagnostic only, not used as ML feature

    def to_array(self) -> np.ndarray:
        """Convert features to numpy array for ML models."""
        return np.array(
            [
                self.peak_magnitude,
                self.peak_freq_idx,
                self.phase_at_peak,
                self.band_power,
                self.spectral_centroid,
                self.phase_kappa,
                self.rayleigh_z,
            ]
        )

    @staticmethod
    def feature_names() -> List[str]:
        """Return feature names."""
        return [
            "peak_magnitude",
            "peak_freq_idx",
            "phase_at_peak",
            "band_power",
            "spectral_centroid",
            "phase_kappa",
            "rayleigh_z",
        ]


@dataclass
class FoldMetrics:
    """Metrics for a single CV fold."""

    fold_idx: int
    repeat_idx: int
    auroc_a: float  # Phase-coherent
    auroc_b: float  # Random-phase
    auprc_a: float
    auprc_b: float
    brier_a: float
    brier_b: float

    @property
    def delta_auroc(self) -> float:
        """AUROC difference (A - B)."""
        return self.auroc_a - self.auroc_b

    @property
    def delta_auprc(self) -> float:
        """AUPRC difference (A - B)."""
        return self.auprc_a - self.auprc_b

    @property
    def delta_brier(self) -> float:
        """Brier score difference (B - A, since lower is better)."""
        return self.brier_b - self.brier_a


# =============================================================================
# Sequence Encoding Functions
# =============================================================================
def encode_sequence(
    seq: str,
    at_lifetime: float = 5.0,
    gc_lifetime: float = 25.0,
    helical_period: float = DEFAULT_HELICAL_PERIOD,
    apply_helical: bool = True,
) -> np.ndarray:
    """
    Encode DNA sequence into complex-valued signal.

    The real component captures breathing kinetics (lifetime), while the
    imaginary component encodes thermodynamic stability (normalized ΔG).

    Args:
        seq: DNA sequence string (ATGC characters)
        at_lifetime: Breathing lifetime for AT pairs [ms]
        gc_lifetime: Breathing lifetime for GC pairs [ms]
        helical_period: DNA helical period [bp/turn]
        apply_helical: Whether to apply helical phase modulation

    Returns:
        Complex numpy array representing the encoded sequence
    """
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        raise ValueError("Empty sequence")

    valid_bases = set("ATGC")
    complex_signal = np.zeros(n, dtype=complex)

    for i in range(n):
        base = seq[i]
        if base not in valid_bases:
            raise ValueError(f"Invalid base at position {i}: {base}")

        # Real part: breathing lifetime
        real_part = at_lifetime if base in "AT" else gc_lifetime

        # Imaginary part: thermodynamic stability (dinucleotide ΔG)
        if i < n - 1:
            dinuc = seq[i : i + 2]
            dg = NEAREST_NEIGHBOR_DG.get(dinuc, -1.5)
        else:
            dg = -1.5  # Neutral fallback for last base

        # Normalize ΔG to [0, 1] range
        norm_imag = (dg - _DG_MIN) / (_DG_MAX - _DG_MIN)
        complex_signal[i] = complex(real_part, norm_imag)

    if apply_helical:
        # Apply helical phase modulation (2π/10.5 per base)
        phases = np.exp(1j * 2 * np.pi * np.arange(n) / helical_period)
        complex_signal *= phases

    return complex_signal


# =============================================================================
# CZT Analysis Functions
# =============================================================================
def compute_czt_spectrum(
    signal: np.ndarray,
    f0: float = 1 / DEFAULT_HELICAL_PERIOD,
    band_width: float = 0.01,
    m: int = 256,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute CZT spectrum focused on helical frequency band.

    Args:
        signal: Complex-valued input signal
        f0: Center frequency (cycles/bp)
        band_width: Half-bandwidth (cycles/bp)
        m: Number of frequency points

    Returns:
        Tuple of (frequencies, complex spectrum)
    """
    f_low = f0 - band_width
    f_high = f0 + band_width
    fs = 1.0  # Spatial sampling at 1 bp

    # CZT parameters: contour from f_low to f_high
    a = np.exp(-2j * np.pi * f_low / fs)
    w = np.exp(2j * np.pi * (f_high - f_low) / (m * fs))

    spectrum = czt(signal, a=a, w=w, m=m)
    freqs = np.linspace(f_low, f_high, m)

    return freqs, spectrum


def extract_czt_features(
    freqs: np.ndarray,
    spectrum: np.ndarray,
) -> CZTFeatures:
    """
    Extract features from CZT spectrum.

    Features include:
    - Peak magnitude and frequency index
    - Phase at peak
    - Band-integrated power
    - Spectral centroid
    - Circular phase concentration (von Mises κ)
    - Rayleigh z-statistic for phase uniformity

    Args:
        freqs: Frequency array
        spectrum: Complex CZT spectrum

    Returns:
        CZTFeatures dataclass with extracted values
    """
    mags = np.abs(spectrum)
    phases = np.angle(spectrum)

    # Peak detection
    peak_idx = np.argmax(mags)
    peak_magnitude = mags[peak_idx]
    phase_at_peak = phases[peak_idx]

    # Band-integrated power
    band_power = np.sum(mags**2)

    # Spectral centroid (weighted mean frequency)
    if np.sum(mags) > 0:
        spectral_centroid = np.sum(freqs * mags) / np.sum(mags)
    else:
        spectral_centroid = freqs[len(freqs) // 2]

    # Circular statistics on phases
    phase_kappa, rayleigh_z, rayleigh_p = _compute_circular_stats(phases)

    return CZTFeatures(
        peak_magnitude=peak_magnitude,
        peak_freq_idx=peak_idx,
        phase_at_peak=phase_at_peak,
        band_power=band_power,
        spectral_centroid=spectral_centroid,
        phase_kappa=phase_kappa,
        rayleigh_z=rayleigh_z,
        rayleigh_p=rayleigh_p,
    )


def _compute_circular_stats(phases: np.ndarray) -> Tuple[float, float, float]:
    """
    Compute circular statistics on phase array.

    Returns:
        Tuple of (von Mises κ, Rayleigh z, Rayleigh p-value)
    """
    n = len(phases)
    if n == 0:
        return 0.0, 0.0, 1.0

    # Resultant vector (mean direction)
    c = np.mean(np.cos(phases))
    s = np.mean(np.sin(phases))
    r_bar = np.sqrt(c**2 + s**2)  # Mean resultant length

    # Von Mises concentration parameter κ estimation
    # Using Mardia & Jupp approximation for κ
    if r_bar < 0.53:
        kappa = 2 * r_bar + r_bar**3 + 5 * r_bar**5 / 6
    elif r_bar < 0.85:
        kappa = -0.4 + 1.39 * r_bar + 0.43 / (1 - r_bar)
    else:
        kappa = 1 / (r_bar**3 - 4 * r_bar**2 + 3 * r_bar)

    # Rayleigh test for uniformity
    rayleigh_z = n * r_bar**2
    # Approximate p-value (valid for large n)
    rayleigh_p = np.exp(-rayleigh_z) * (
        1
        + (2 * rayleigh_z - rayleigh_z**2) / (4 * n)
        - (
            24 * rayleigh_z
            - 132 * rayleigh_z**2
            + 76 * rayleigh_z**3
            - 9 * rayleigh_z**4
        )
        / (288 * n**2)
    )
    rayleigh_p = max(0.0, min(1.0, rayleigh_p))

    return kappa, rayleigh_z, rayleigh_p


def scramble_phases(
    spectrum: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Create random-phase control spectrum.

    Preserves magnitudes but replaces phases with uniform random values
    in [-π, π).

    Args:
        spectrum: Complex CZT spectrum
        rng: NumPy random generator for reproducibility

    Returns:
        Phase-scrambled complex spectrum
    """
    mags = np.abs(spectrum)
    random_phases = rng.uniform(-np.pi, np.pi, len(spectrum))
    return mags * np.exp(1j * random_phases)


# =============================================================================
# Data Loading Functions
# =============================================================================
def load_sequences_and_labels(
    input_path: Path,
    labels_path: Optional[Path] = None,
) -> Tuple[List[str], List[str], List[int]]:
    """
    Load sequences and labels from FASTA file or CSV.

    Supports two formats:
    1. FASTA with embedded labels: >seq_id|label
    2. Separate labels CSV with columns: id, label

    Args:
        input_path: Path to input FASTA/CSV file
        labels_path: Optional path to labels CSV

    Returns:
        Tuple of (sequence_ids, sequences, binary_labels)
    """
    seq_ids: List[str] = []
    sequences: List[str] = []
    labels_str: List[str] = []

    suffix = input_path.suffix.lower()

    if suffix in [".fasta", ".fa"]:
        with open(input_path, "r") as f:
            current_id = ""
            current_seq = ""
            current_label = "unknown"

            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        seq_ids.append(current_id)
                        sequences.append(current_seq)
                        labels_str.append(current_label)
                    current_seq = ""

                    header = line[1:].strip()
                    if "|" in header:
                        parts = header.split("|", 1)
                        current_id = parts[0]
                        current_label = parts[1] if len(parts) > 1 else "unknown"
                    else:
                        current_id = header
                        current_label = "unknown"
                else:
                    current_seq += line

            if current_seq:
                seq_ids.append(current_id)
                sequences.append(current_seq)
                labels_str.append(current_label)

    else:
        # Assume CSV format: id,sequence,label or sequence,label
        with open(input_path, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if not row:
                    continue
                if len(row) >= 3:
                    # Format: id, sequence, label
                    seq_ids.append(row[0])
                    sequences.append(row[1])
                    labels_str.append(row[2])
                elif len(row) == 2:
                    # Format: sequence, label
                    seq_ids.append(str(len(sequences)))
                    sequences.append(row[0])
                    labels_str.append(row[1])

    # If external labels provided, override
    if labels_path is not None:
        label_map: Dict[str, str] = {}
        with open(labels_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                label_map[row["id"]] = row["label"]

        for i, sid in enumerate(seq_ids):
            if sid in label_map:
                labels_str[i] = label_map[sid]

    # Convert string labels to binary
    unique_labels = sorted(set(labels_str))
    if len(unique_labels) < 2:
        raise ValueError(
            f"Need at least 2 distinct labels for binary classification, "
            f"found: {unique_labels}"
        )

    # Map to 0/1 (alphabetically sorted)
    label_to_int = {lbl: i for i, lbl in enumerate(unique_labels[:2])}
    binary_labels = [label_to_int.get(lbl, 0) for lbl in labels_str]

    return seq_ids, sequences, binary_labels


# =============================================================================
# Feature Extraction Pipeline
# =============================================================================
def extract_features_for_sequences(
    sequences: List[str],
    config: StudyConfig,
    phase_rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract phase-coherent and random-phase features for all sequences.

    Args:
        sequences: List of DNA sequences
        config: Study configuration
        phase_rng: Random generator for phase scrambling

    Returns:
        Tuple of (X_coherent, X_random) feature arrays
    """
    features_coherent: List[np.ndarray] = []
    features_random: List[np.ndarray] = []

    for seq in sequences:
        try:
            signal = encode_sequence(seq)
            freqs, spectrum = compute_czt_spectrum(
                signal,
                f0=config.band_center,
                band_width=config.band_width,
                m=config.czt_points,
            )

            # Phase-coherent features
            feat_coherent = extract_czt_features(freqs, spectrum)
            features_coherent.append(feat_coherent.to_array())

            # Random-phase control
            scrambled_spectrum = scramble_phases(spectrum, phase_rng)
            feat_random = extract_czt_features(freqs, scrambled_spectrum)
            features_random.append(feat_random.to_array())

        except Exception as e:
            warnings.warn(f"Error processing sequence: {e}")
            # Use zeros for failed sequences
            n_features = len(CZTFeatures.feature_names())
            features_coherent.append(np.zeros(n_features))
            features_random.append(np.zeros(n_features))

    X_coherent = np.array(features_coherent)
    X_random = np.array(features_random)

    return X_coherent, X_random


# =============================================================================
# Cross-Validation and Evaluation
# =============================================================================
def run_cv_evaluation(
    X_coherent: np.ndarray,
    X_random: np.ndarray,
    y: np.ndarray,
    config: StudyConfig,
) -> List[FoldMetrics]:
    """
    Run stratified repeated K-fold cross-validation.

    Uses same splits for both conditions (phase-coherent and random-phase).

    Args:
        X_coherent: Feature matrix (phase-coherent)
        X_random: Feature matrix (random-phase)
        y: Binary labels
        config: Study configuration

    Returns:
        List of FoldMetrics for each fold
    """
    cv = RepeatedStratifiedKFold(
        n_splits=config.n_folds,
        n_repeats=config.n_repeats,
        random_state=config.cv_seed,
    )

    fold_metrics: List[FoldMetrics] = []
    fold_idx = 0

    for train_idx, test_idx in cv.split(X_coherent, y):
        repeat_idx = fold_idx // config.n_folds
        fold_in_repeat = fold_idx % config.n_folds

        X_train_a, X_test_a = X_coherent[train_idx], X_coherent[test_idx]
        X_train_b, X_test_b = X_random[train_idx], X_random[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Skip if not enough samples or all same class
        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            fold_idx += 1
            continue

        # Standardize features (fit on train, transform both)
        scaler_a = StandardScaler()
        scaler_b = StandardScaler()
        X_train_a_scaled = scaler_a.fit_transform(X_train_a)
        X_test_a_scaled = scaler_a.transform(X_test_a)
        X_train_b_scaled = scaler_b.fit_transform(X_train_b)
        X_test_b_scaled = scaler_b.transform(X_test_b)

        # Train logistic regression (no hyperparameter tuning)
        clf_a = LogisticRegression(max_iter=1000, random_state=config.global_seed)
        clf_b = LogisticRegression(max_iter=1000, random_state=config.global_seed)

        try:
            clf_a.fit(X_train_a_scaled, y_train)
            clf_b.fit(X_train_b_scaled, y_train)
        except Exception:
            fold_idx += 1
            continue

        y_prob_a = clf_a.predict_proba(X_test_a_scaled)[:, 1]
        y_prob_b = clf_b.predict_proba(X_test_b_scaled)[:, 1]

        # Compute metrics
        auroc_a = roc_auc_score(y_test, y_prob_a)
        auroc_b = roc_auc_score(y_test, y_prob_b)
        auprc_a = average_precision_score(y_test, y_prob_a)
        auprc_b = average_precision_score(y_test, y_prob_b)
        brier_a = brier_score_loss(y_test, y_prob_a)
        brier_b = brier_score_loss(y_test, y_prob_b)

        metrics = FoldMetrics(
            fold_idx=fold_in_repeat,
            repeat_idx=repeat_idx,
            auroc_a=auroc_a,
            auroc_b=auroc_b,
            auprc_a=auprc_a,
            auprc_b=auprc_b,
            brier_a=brier_a,
            brier_b=brier_b,
        )
        fold_metrics.append(metrics)
        fold_idx += 1

    return fold_metrics


# =============================================================================
# Statistical Analysis
# =============================================================================
def compute_paired_statistics(
    fold_metrics: List[FoldMetrics],
    config: StudyConfig,
) -> Dict[str, Any]:
    """
    Compute paired statistics for condition comparison.

    Includes:
    - Paired t-test
    - Wilcoxon signed-rank test
    - Cohen's d (paired) with bootstrap CI
    - Effect size interpretation

    Args:
        fold_metrics: List of per-fold metrics
        config: Study configuration

    Returns:
        Dictionary with statistical results
    """
    if not fold_metrics:
        return {"error": "No fold metrics to analyze"}

    deltas_auroc = [m.delta_auroc for m in fold_metrics]
    deltas_auprc = [m.delta_auprc for m in fold_metrics]
    deltas_brier = [m.delta_brier for m in fold_metrics]

    results: Dict[str, Any] = {}

    for metric_name, deltas in [
        ("auroc", deltas_auroc),
        ("auprc", deltas_auprc),
        ("brier", deltas_brier),
    ]:
        deltas_arr = np.array(deltas)
        n = len(deltas_arr)

        if n < 2:
            results[metric_name] = {"error": "Insufficient data"}
            continue

        # Mean and std of deltas
        mean_delta = float(np.mean(deltas_arr))
        std_delta = float(np.std(deltas_arr, ddof=1))

        # Paired t-test
        t_stat, t_pval = stats.ttest_1samp(deltas_arr, 0)

        # Wilcoxon signed-rank test
        try:
            w_stat, w_pval = stats.wilcoxon(deltas_arr, zero_method="wilcox")
        except ValueError:
            # All zeros or insufficient variation
            w_stat, w_pval = np.nan, 1.0

        # Cohen's d (paired): d = mean(delta) / std(delta)
        cohens_d = mean_delta / std_delta if std_delta > 0 else 0.0

        # Bootstrap CI for Cohen's d
        rng = np.random.default_rng(config.global_seed)
        bootstrap_d = []
        for _ in range(config.n_bootstrap):
            boot_sample = rng.choice(deltas_arr, size=n, replace=True)
            boot_mean = np.mean(boot_sample)
            boot_std = np.std(boot_sample, ddof=1)
            if boot_std > 0:
                bootstrap_d.append(boot_mean / boot_std)
            else:
                bootstrap_d.append(0.0)

        ci_low = float(np.percentile(bootstrap_d, 2.5))
        ci_high = float(np.percentile(bootstrap_d, 97.5))

        # Effect size interpretation
        abs_d = abs(cohens_d)
        if abs_d < 0.2:
            effect_interpretation = "negligible"
        elif abs_d < 0.5:
            effect_interpretation = "small"
        elif abs_d < 0.8:
            effect_interpretation = "medium"
        else:
            effect_interpretation = "large"

        results[metric_name] = {
            "mean_delta": mean_delta,
            "std_delta": std_delta,
            "n_folds": n,
            "t_statistic": float(t_stat),
            "t_pvalue": float(t_pval),
            "wilcoxon_statistic": float(w_stat) if not np.isnan(w_stat) else None,
            "wilcoxon_pvalue": float(w_pval),
            "cohens_d": cohens_d,
            "cohens_d_ci_low": ci_low,
            "cohens_d_ci_high": ci_high,
            "effect_interpretation": effect_interpretation,
        }

    return results


def apply_fdr_correction(
    p_values: List[float],
    alpha: float = 0.05,
) -> Tuple[List[bool], List[float]]:
    """
    Apply Benjamini-Hochberg FDR correction.

    Args:
        p_values: List of p-values
        alpha: FDR threshold

    Returns:
        Tuple of (significant indicators, adjusted p-values)
    """
    n = len(p_values)
    if n == 0:
        return [], []

    # Sort p-values with indices
    sorted_pairs = sorted(enumerate(p_values), key=lambda x: x[1])

    # Compute adjusted p-values
    adjusted = [0.0] * n
    prev_adj = 1.0
    for i in range(n - 1, -1, -1):
        orig_idx, pval = sorted_pairs[i]
        adj = min(prev_adj, pval * n / (i + 1))
        adjusted[orig_idx] = adj
        prev_adj = adj

    significant = [adj <= alpha for adj in adjusted]
    return significant, adjusted


def run_permutation_test(
    X_coherent: np.ndarray,
    X_random: np.ndarray,
    y: np.ndarray,
    config: StudyConfig,
) -> Dict[str, Any]:
    """
    Run label permutation test.

    Permutes labels and reruns CV to establish null distribution.
    The phase-coherence advantage should vanish under permutation.

    Args:
        X_coherent: Feature matrix (phase-coherent)
        X_random: Feature matrix (random-phase)
        y: Binary labels
        config: Study configuration

    Returns:
        Dictionary with permutation test results
    """
    rng = np.random.default_rng(config.global_seed + 1000)

    # Observed delta
    observed_metrics = run_cv_evaluation(X_coherent, X_random, y, config)
    if not observed_metrics:
        return {"error": "No observed metrics"}

    observed_delta = np.mean([m.delta_auroc for m in observed_metrics])

    # Permutation distribution
    perm_deltas: List[float] = []
    for perm_i in range(config.n_permutations):
        y_perm = rng.permutation(y)
        perm_metrics = run_cv_evaluation(X_coherent, X_random, y_perm, config)
        if perm_metrics:
            perm_delta = np.mean([m.delta_auroc for m in perm_metrics])
            perm_deltas.append(perm_delta)

    if not perm_deltas:
        return {"error": "No permutation results"}

    # Compute p-value
    perm_deltas_arr = np.array(perm_deltas)
    n_extreme = np.sum(perm_deltas_arr >= observed_delta)
    p_value = (n_extreme + 1) / (len(perm_deltas_arr) + 1)

    return {
        "observed_delta_auroc": observed_delta,
        "permutation_mean": float(np.mean(perm_deltas_arr)),
        "permutation_std": float(np.std(perm_deltas_arr)),
        "permutation_pvalue": p_value,
        "n_permutations": len(perm_deltas_arr),
    }


def run_ablation_study(
    sequences: List[str],
    y: np.ndarray,
    config: StudyConfig,
) -> Dict[str, Any]:
    """
    Run phase-only vs magnitude-only ablation study.

    Args:
        sequences: List of DNA sequences
        y: Binary labels
        config: Study configuration

    Returns:
        Dictionary with ablation results
    """
    phase_features: List[np.ndarray] = []
    mag_features: List[np.ndarray] = []

    for seq in sequences:
        try:
            signal = encode_sequence(seq)
            freqs, spectrum = compute_czt_spectrum(
                signal,
                f0=config.band_center,
                band_width=config.band_width,
                m=config.czt_points,
            )

            phases = np.angle(spectrum)
            mags = np.abs(spectrum)

            # Phase-only: use phase values directly
            phase_feat = np.array(
                [
                    np.mean(phases),
                    np.std(phases),
                    np.max(phases) - np.min(phases),
                ]
            )

            # Magnitude-only: use magnitude values directly
            mag_feat = np.array(
                [
                    np.max(mags),
                    np.mean(mags),
                    np.sum(mags**2),
                ]
            )

            phase_features.append(phase_feat)
            mag_features.append(mag_feat)

        except Exception:
            phase_features.append(np.zeros(3))
            mag_features.append(np.zeros(3))

    X_phase = np.array(phase_features)
    X_mag = np.array(mag_features)

    cv = RepeatedStratifiedKFold(
        n_splits=config.n_folds,
        n_repeats=1,
        random_state=config.cv_seed,
    )

    aurocs_phase: List[float] = []
    aurocs_mag: List[float] = []

    for train_idx, test_idx in cv.split(X_phase, y):
        y_train, y_test = y[train_idx], y[test_idx]

        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            continue

        try:
            scaler_phase = StandardScaler()
            X_phase_train = scaler_phase.fit_transform(X_phase[train_idx])
            X_phase_test = scaler_phase.transform(X_phase[test_idx])
            clf_phase = LogisticRegression(
                max_iter=1000, random_state=config.global_seed
            )
            clf_phase.fit(X_phase_train, y_train)
            y_prob_phase = clf_phase.predict_proba(X_phase_test)[:, 1]
            aurocs_phase.append(roc_auc_score(y_test, y_prob_phase))
        except Exception:
            pass

        try:
            scaler_mag = StandardScaler()
            X_mag_train = scaler_mag.fit_transform(X_mag[train_idx])
            X_mag_test = scaler_mag.transform(X_mag[test_idx])
            clf_mag = LogisticRegression(max_iter=1000, random_state=config.global_seed)
            clf_mag.fit(X_mag_train, y_train)
            y_prob_mag = clf_mag.predict_proba(X_mag_test)[:, 1]
            aurocs_mag.append(roc_auc_score(y_test, y_prob_mag))
        except Exception:
            pass

    return {
        "phase_only_auroc_mean": float(np.mean(aurocs_phase)) if aurocs_phase else None,
        "phase_only_auroc_std": float(np.std(aurocs_phase)) if aurocs_phase else None,
        "magnitude_only_auroc_mean": float(np.mean(aurocs_mag)) if aurocs_mag else None,
        "magnitude_only_auroc_std": float(np.std(aurocs_mag)) if aurocs_mag else None,
    }


def compute_power_analysis(
    fold_metrics: List[FoldMetrics],
    config: StudyConfig,
) -> Dict[str, Any]:
    """
    Compute post-hoc power analysis.

    Args:
        fold_metrics: List of per-fold metrics
        config: Study configuration

    Returns:
        Dictionary with power analysis results
    """
    if len(fold_metrics) < 2:
        return {"error": "Insufficient data for power analysis"}

    deltas = [m.delta_auroc for m in fold_metrics]
    n = len(deltas)
    observed_delta = np.mean(deltas)
    observed_std = np.std(deltas, ddof=1)

    if observed_std == 0:
        return {"error": "Zero variance in deltas"}

    # Effect size (Cohen's d for one-sample)
    observed_d = abs(observed_delta) / observed_std

    # Non-centrality parameter
    ncp = observed_d * np.sqrt(n)

    # Critical t-value
    t_crit = stats.t.ppf(1 - config.power_alpha / 2, df=n - 1)

    # Power (probability of rejecting H0 when H1 is true)
    power = (
        1
        - stats.nct.cdf(t_crit, df=n - 1, nc=ncp)
        + stats.nct.cdf(-t_crit, df=n - 1, nc=ncp)
    )

    # Sample size needed for desired power
    # Using approximation: n ≈ (z_α + z_β)² / d²
    z_alpha = stats.norm.ppf(1 - config.power_alpha / 2)
    z_beta = stats.norm.ppf(config.desired_power)
    if observed_d > 0:
        n_needed = ((z_alpha + z_beta) / observed_d) ** 2
    else:
        n_needed = float("inf")

    return {
        "n_observations": n,
        "observed_delta": observed_delta,
        "observed_std": observed_std,
        "observed_effect_size": observed_d,
        "achieved_power": power,
        "desired_power": config.desired_power,
        "alpha": config.power_alpha,
        "n_needed_for_desired_power": (
            int(np.ceil(n_needed)) if np.isfinite(n_needed) else None
        ),
    }


# =============================================================================
# Output Generation
# =============================================================================
def save_fold_metrics_csv(
    fold_metrics: List[FoldMetrics],
    config: StudyConfig,
    output_path: Path,
) -> None:
    """Save per-fold metrics to CSV."""
    fieldnames = [
        "fold_idx",
        "repeat_idx",
        "auroc_a",
        "auroc_b",
        "auprc_a",
        "auprc_b",
        "brier_a",
        "brier_b",
        "delta_auroc",
        "delta_auprc",
        "delta_brier",
        "global_seed",
        "phase_seed",
        "cv_seed",
        "band_center",
        "band_width",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for m in fold_metrics:
            writer.writerow(
                {
                    "fold_idx": m.fold_idx,
                    "repeat_idx": m.repeat_idx,
                    "auroc_a": f"{m.auroc_a:.6f}",
                    "auroc_b": f"{m.auroc_b:.6f}",
                    "auprc_a": f"{m.auprc_a:.6f}",
                    "auprc_b": f"{m.auprc_b:.6f}",
                    "brier_a": f"{m.brier_a:.6f}",
                    "brier_b": f"{m.brier_b:.6f}",
                    "delta_auroc": f"{m.delta_auroc:.6f}",
                    "delta_auprc": f"{m.delta_auprc:.6f}",
                    "delta_brier": f"{m.delta_brier:.6f}",
                    "global_seed": config.global_seed,
                    "phase_seed": config.phase_seed,
                    "cv_seed": config.cv_seed,
                    "band_center": f"{config.band_center:.6f}",
                    "band_width": f"{config.band_width:.6f}",
                }
            )


def save_summary_json(
    results: Dict[str, Any],
    output_path: Path,
) -> None:
    """Save summary results to JSON."""

    def convert_numpy(obj: Any) -> Any:
        """Convert numpy types to Python natives for JSON serialization."""
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj) if isinstance(obj, np.floating) else int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_numpy(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_numpy(v) for v in obj]
        return obj

    with open(output_path, "w") as f:
        json.dump(convert_numpy(results), f, indent=2)


# =============================================================================
# Plotting Functions
# =============================================================================
def generate_plots(
    fold_metrics: List[FoldMetrics],
    stats_results: Dict[str, Any],
    output_dir: Path,
) -> None:
    """
    Generate diagnostic plots.

    Produces:
    - Paired spaghetti plot of per-fold AUROC
    - Violin/box plot of deltas
    - Phase circular histogram (if phase data available)
    """
    if not MATPLOTLIB_AVAILABLE:
        warnings.warn("matplotlib not available; skipping plots")
        return

    # Spaghetti plot: paired AUROC
    fig, ax = plt.subplots(figsize=(8, 6))
    for m in fold_metrics:
        ax.plot([0, 1], [m.auroc_a, m.auroc_b], "o-", color="gray", alpha=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Phase-Coherent (A)", "Random-Phase (B)"])
    ax.set_ylabel("AUROC")
    ax.set_title("Paired AUROC: Phase-Coherent vs Random-Phase")
    ax.axhline(y=0.5, color="red", linestyle="--", alpha=0.5, label="Chance")
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "spaghetti_auroc.png", dpi=150)
    plt.close()

    # Violin/box plot of deltas
    deltas_auroc = [m.delta_auroc for m in fold_metrics]
    deltas_auprc = [m.delta_auprc for m in fold_metrics]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].boxplot(deltas_auroc, vert=True)
    axes[0].axhline(y=0, color="red", linestyle="--", alpha=0.5)
    axes[0].set_title("ΔAUROC (A - B)")
    axes[0].set_ylabel("Delta")

    axes[1].boxplot(deltas_auprc, vert=True)
    axes[1].axhline(y=0, color="red", linestyle="--", alpha=0.5)
    axes[1].set_title("ΔAUPRC (A - B)")
    axes[1].set_ylabel("Delta")

    plt.tight_layout()
    plt.savefig(output_dir / "delta_boxplot.png", dpi=150)
    plt.close()

    print(f"Plots saved to {output_dir}")


# =============================================================================
# Main Entry Point
# =============================================================================
def run_study(
    input_path: Path,
    labels_path: Optional[Path],
    output_dir: Path,
    config: StudyConfig,
) -> Dict[str, Any]:
    """
    Run the complete phase coherence study.

    Args:
        input_path: Path to input FASTA/CSV
        labels_path: Optional path to labels CSV
        output_dir: Output directory
        config: Study configuration

    Returns:
        Dictionary with all results
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase-Coherent vs Random-Phase CZT Spectra Comparison Study")
    print("=" * 60)
    print()
    print("Configuration:")
    print(f"  Global Seed: {config.global_seed}")
    print(f"  Phase Seed: {config.phase_seed}")
    print(f"  CV Seed: {config.cv_seed}")
    print(f"  K-Folds: {config.n_folds}")
    print(f"  Repeats: {config.n_repeats}")
    print(f"  Band Center: {config.band_center:.6f} cycles/bp")
    print(f"  Band Width: {config.band_width:.6f} cycles/bp")
    print()

    # Load data
    print("Loading sequences...")
    seq_ids, sequences, labels = load_sequences_and_labels(input_path, labels_path)
    y = np.array(labels)
    n_samples = len(sequences)
    print(f"  Loaded {n_samples} sequences")
    print(f"  Label distribution: {dict(zip(*np.unique(y, return_counts=True)))}")
    print()

    # Apply pilot fraction if specified
    if config.pilot_fraction is not None and 0 < config.pilot_fraction < 1:
        n_pilot = int(n_samples * config.pilot_fraction)
        rng = np.random.default_rng(config.global_seed)
        pilot_idx = rng.choice(n_samples, size=n_pilot, replace=False)
        sequences = [sequences[i] for i in pilot_idx]
        y = y[pilot_idx]
        print(
            f"  Pilot study: using {n_pilot} sequences ({config.pilot_fraction * 100:.0f}%)"
        )
        print()

    # Extract features
    print("Extracting CZT features...")
    phase_rng = np.random.default_rng(config.phase_seed)
    X_coherent, X_random = extract_features_for_sequences(sequences, config, phase_rng)
    print(f"  Features shape: {X_coherent.shape}")
    print()

    # Run cross-validation
    print("Running cross-validation...")
    fold_metrics = run_cv_evaluation(X_coherent, X_random, y, config)
    print(f"  Completed {len(fold_metrics)} folds")
    print()

    # Compute statistics
    print("Computing statistics...")
    stats_results = compute_paired_statistics(fold_metrics, config)

    # Power analysis
    print("Running power analysis...")
    power_results = compute_power_analysis(fold_metrics, config)

    # Ablation study
    print("Running ablation study (phase-only vs magnitude-only)...")
    ablation_results = run_ablation_study(sequences, y, config)

    # Label permutation test (can be slow)
    print(f"Running label permutation test ({config.n_permutations} permutations)...")
    permutation_config = StudyConfig(
        global_seed=config.global_seed,
        phase_seed=config.phase_seed,
        cv_seed=config.cv_seed,
        n_folds=config.n_folds,
        n_repeats=1,  # Reduce for permutation test
        band_center=config.band_center,
        band_width=config.band_width,
        n_permutations=config.n_permutations,
    )
    perm_results = run_permutation_test(X_coherent, X_random, y, permutation_config)

    # Compile results
    results: Dict[str, Any] = {
        "study_design": {
            "n_samples": len(sequences),
            "n_folds": config.n_folds,
            "n_repeats": config.n_repeats,
            "band_center": config.band_center,
            "band_width": config.band_width,
            "seeds": {
                "global": config.global_seed,
                "phase": config.phase_seed,
                "cv": config.cv_seed,
            },
        },
        "statistics": stats_results,
        "power_analysis": power_results,
        "ablation": ablation_results,
        "permutation_test": perm_results,
    }

    # Print summary
    print()
    print("=" * 60)
    print("Results Summary")
    print("=" * 60)

    if "auroc" in stats_results:
        auroc_stats = stats_results["auroc"]
        print()
        print("AUROC Analysis (Primary Metric):")
        print(f"  Mean Δ (A-B): {auroc_stats.get('mean_delta', 0):.4f}")
        print(
            f"  Cohen's d: {auroc_stats.get('cohens_d', 0):.4f} "
            f"[{auroc_stats.get('cohens_d_ci_low', 0):.4f}, "
            f"{auroc_stats.get('cohens_d_ci_high', 0):.4f}]"
        )
        print(f"  Effect: {auroc_stats.get('effect_interpretation', 'N/A')}")
        print(f"  Paired t-test p: {auroc_stats.get('t_pvalue', 1):.6f}")
        print(f"  Wilcoxon p: {auroc_stats.get('wilcoxon_pvalue', 1):.6f}")

    if power_results and "achieved_power" in power_results:
        print()
        print("Power Analysis:")
        print(f"  Achieved power: {power_results['achieved_power']:.4f}")
        if power_results.get("n_needed_for_desired_power"):
            print(
                f"  N needed for 80% power: {power_results['n_needed_for_desired_power']}"
            )

    if perm_results and "permutation_pvalue" in perm_results:
        print()
        print("Permutation Test:")
        print(f"  Observed Δ AUROC: {perm_results['observed_delta_auroc']:.4f}")
        print(f"  Permutation p-value: {perm_results['permutation_pvalue']:.4f}")

    print()

    # Save outputs
    csv_path = output_dir / "fold_metrics.csv"
    save_fold_metrics_csv(fold_metrics, config, csv_path)
    print(f"Saved fold metrics to: {csv_path}")

    json_path = output_dir / "summary.json"
    save_summary_json(results, json_path)
    print(f"Saved summary to: {json_path}")

    # Generate plots
    if not config.skip_plots:
        generate_plots(fold_metrics, stats_results, output_dir)

    return results


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Phase-Coherent vs Random-Phase CZT Spectra Comparison Study",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python experiments/phase_coherence_study.py \\
        --input data/sequences.fasta \\
        --output-dir results/phase_study/

The study compares phase-coherent CZT features (Condition A) against
random-phase controls (Condition B) using stratified K-fold cross-validation
with logistic regression classification.
        """,
    )

    # Input/output arguments
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to FASTA file with sequences (labels in headers: >id|label)",
    )
    parser.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="Optional path to CSV with columns: id, label",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/phase_coherence_study"),
        help="Directory for output artifacts",
    )

    # Seed arguments
    parser.add_argument(
        "--global-seed",
        type=int,
        default=DEFAULT_GLOBAL_SEED,
        help=f"Global random seed (default: {DEFAULT_GLOBAL_SEED})",
    )
    parser.add_argument(
        "--phase-seed",
        type=int,
        default=DEFAULT_PHASE_SEED,
        help=f"Phase scrambling seed (default: {DEFAULT_PHASE_SEED})",
    )
    parser.add_argument(
        "--cv-seed",
        type=int,
        default=DEFAULT_CV_SEED,
        help=f"Cross-validation seed (default: {DEFAULT_CV_SEED})",
    )

    # CV parameters
    parser.add_argument(
        "--n-folds",
        type=int,
        default=5,
        help="Number of CV folds (default: 5)",
    )
    parser.add_argument(
        "--n-repeats",
        type=int,
        default=5,
        help="Number of CV repeats (default: 5)",
    )

    # CZT band parameters
    parser.add_argument(
        "--band-center",
        type=float,
        default=1 / DEFAULT_HELICAL_PERIOD,
        help=f"CZT band center frequency in cycles/bp (default: 1/{DEFAULT_HELICAL_PERIOD})",
    )
    parser.add_argument(
        "--band-width",
        type=float,
        default=0.01,
        help="CZT band half-width in cycles/bp (default: 0.01)",
    )

    # Study options
    parser.add_argument(
        "--pilot-fraction",
        type=float,
        default=None,
        help="Fraction of data to use for pilot study (0-1)",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=10000,
        help="Number of bootstrap iterations for CI (default: 10000)",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=1000,
        help="Number of label permutations (default: 1000)",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip plot generation",
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if args.labels and not args.labels.exists():
        print(f"Error: Labels file not found: {args.labels}", file=sys.stderr)
        sys.exit(1)

    # Build configuration
    config = StudyConfig(
        global_seed=args.global_seed,
        phase_seed=args.phase_seed,
        cv_seed=args.cv_seed,
        n_folds=args.n_folds,
        n_repeats=args.n_repeats,
        band_center=args.band_center,
        band_width=args.band_width,
        pilot_fraction=args.pilot_fraction,
        n_bootstrap=args.n_bootstrap,
        n_permutations=args.n_permutations,
        skip_plots=args.skip_plots,
    )

    # Run study
    try:
        run_study(args.input, args.labels, args.output_dir, config)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
