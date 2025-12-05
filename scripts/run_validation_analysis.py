#!/usr/bin/env python3
"""
DNA Breathing Dynamics Validation Analysis Script.

This script computes scientifically valid statistics (effect sizes,
confidence intervals) on spectral features derived from DNA breathing
dynamics in real CRISPR guide RNA sequences.

Usage:
    python scripts/run_validation_analysis.py

Outputs:
    - params.yaml: Analysis parameters
    - features.csv: Per-sequence spectral features
    - reports/validation_stats.md: Summary statistics report

Reference: Issue #44 - Revised Test Plan for DNA Breathing Dynamics Analysis
"""

# Standard library imports
import csv
import random
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Third-party imports
import numpy as np
from Bio import SeqIO
from scipy import stats
from scipy.signal import czt

# =============================================================================
# Constants and Parameters
# =============================================================================

# CRISPR guide RNA sequence length (standard for Cas9)
CRISPR_GUIDE_LENGTH = 20

# SantaLucia 1998 Unified nearest-neighbor parameters (kcal/mol)
NEAREST_NEIGHBOR_DG = {
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

# Thermodynamic constants derived from nearest-neighbor table
DG_MIN = min(NEAREST_NEIGHBOR_DG.values())  # Most stable (most negative)
DG_MAX = max(NEAREST_NEIGHBOR_DG.values())  # Least stable

# Default ΔG value for unpaired/terminal bases (neutral, near mean of NN values)
DEFAULT_DG_VALUE = -1.5

# Default noise floor for SNR calculations when no off-band data available
DEFAULT_NOISE_FLOOR = 1.0

# Minimum bootstrap samples for valid confidence interval estimation
MIN_BOOTSTRAP_SAMPLES = 50


# =============================================================================
# Encoding Functions
# =============================================================================


def compute_gc_content(seq: str) -> float:
    """Calculate GC content fraction for a sequence."""
    seq = seq.upper()
    gc_count = sum(1 for b in seq if b in "GC")
    return gc_count / len(seq) if len(seq) > 0 else 0.0


def encode_sequence(
    seq: str,
    at_lifetime: float = 1.0,
    gc_lifetime: float = 50.0,
    helical_period: float = 10.5,
    apply_helical: bool = True,
) -> np.ndarray:
    """
    Encode DNA sequence as complex-valued signal for spectral analysis.

    Real part: Breathing kinetics (base pair lifetimes)
        AT pairs: ~1 ms lifetime (faster opening, 2 H-bonds)
        GC pairs: ~50 ms lifetime (slower opening, 3 H-bonds)

    Imaginary part: Thermodynamic stability (ΔG from nearest-neighbor model)

    Args:
        seq: DNA sequence (uppercase ACGT)
        at_lifetime: AT pair opening lifetime (ms)
        gc_lifetime: GC pair opening lifetime (ms)
        helical_period: B-form DNA period (10.5 bp)
        apply_helical: Apply helical phase modulation

    Returns:
        Complex-valued signal array
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

        # Real part: breathing kinetics
        real_part = at_lifetime if base in "AT" else gc_lifetime

        # Imaginary part: thermodynamic stability from dinucleotide
        if i < n - 1:
            dinuc = seq[i : i + 2]
            dg = NEAREST_NEIGHBOR_DG.get(dinuc, DEFAULT_DG_VALUE)
        else:
            dg = DEFAULT_DG_VALUE  # Neutral for last base

        # Normalize DG to 0..1 range
        norm_imag = (dg - DG_MIN) / (DG_MAX - DG_MIN)
        complex_signal[i] = complex(real_part, norm_imag)

    if apply_helical:
        phases = np.exp(1j * 2 * np.pi * np.arange(n) / helical_period)
        complex_signal *= phases

    return complex_signal


# =============================================================================
# Spectral Analysis Functions
# =============================================================================


def czt_analysis(
    signal: np.ndarray,
    f0: float = 1 / 10.5,
    band_width: float = 0.01,
    noise_band_width: float = 0.03,
    M: int = 256,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute Chirp Z-Transform focused on helical frequency.

    Computes a wider frequency band than the analysis band to allow
    proper noise floor estimation from off-band frequencies.

    Args:
        signal: Complex-valued input signal
        f0: Center frequency (1/10.5 bp⁻¹ for helical period)
        band_width: Frequency band width around f0 for signal analysis
        noise_band_width: Extended band for noise floor estimation
        M: Number of frequency points

    Returns:
        Tuple of (frequencies, spectrum)
    """
    fs = 1.0  # Sampling frequency (1 sample per base)
    # Compute wider band to capture noise floor
    f_low = f0 - noise_band_width
    f_high = f0 + noise_band_width

    # CZT spiral parameters (matching scipy.signal.czt conventions)
    a = np.exp(-2j * np.pi * f_low / fs)
    w = np.exp(2j * np.pi * (f_high - f_low) / (M * fs))

    spectrum = czt(signal, a=a, w=w, m=M)
    freqs = np.linspace(f_low, f_high, M)

    return freqs, spectrum


def extract_features(
    freqs: np.ndarray,
    spectrum: np.ndarray,
    f0: float = 1 / 10.5,
    band_width: float = 0.01,
    guard_band: float = 0.005,
) -> Dict[str, Any]:
    """
    Extract spectral features from CZT spectrum.

    Args:
        freqs: Frequency array
        spectrum: Complex spectrum array
        f0: Target frequency (helical period)
        band_width: Width of analysis band
        guard_band: Guard band for skirt estimation

    Returns:
        Dictionary of spectral features
    """
    mags = np.abs(spectrum)
    phases = np.angle(spectrum)

    # In-band selection
    in_band = (freqs >= f0 - band_width) & (freqs <= f0 + band_width)
    band_mags = mags[in_band]
    band_phases = phases[in_band]

    if len(band_mags) == 0:
        raise ValueError("No points in frequency band")

    # Peak detection
    peak_local_idx = np.argmax(band_mags)
    peak_idx = np.where(in_band)[0][peak_local_idx]
    peak_mag = float(band_mags[peak_local_idx])
    peak_freq = float(freqs[peak_idx])
    peak_phase = float(band_phases[peak_local_idx])

    # Skirt estimation (off-band noise floor)
    # Use frequencies outside the signal band but within the computed spectrum
    edge_low = f0 - band_width - guard_band
    edge_high = f0 + band_width + guard_band
    off_mask = (freqs < edge_low) | (freqs > edge_high)
    skirt_mags = mags[off_mask]

    # Compute noise floor from off-band frequencies
    if len(skirt_mags) > 0:
        skirt_mean = float(np.mean(skirt_mags))
    else:
        # Fallback if no off-band data (should not happen with wider CZT)
        skirt_mean = DEFAULT_NOISE_FLOOR

    # Feature calculations with finite value guarantees
    band_mean = float(np.mean(band_mags))
    band_energy = float(np.sum(band_mags**2))
    coherence = float(np.abs(np.mean(np.exp(1j * band_phases))))

    # Use large but finite values instead of infinity for ratios
    max_ratio = 1e6  # Large but finite upper bound
    if skirt_mean > 0:
        peak_to_skirt = min(peak_mag / skirt_mean, max_ratio)
        snr = min(band_mean / skirt_mean, max_ratio)
    else:
        peak_to_skirt = max_ratio
        snr = max_ratio

    return {
        "peak_mag": peak_mag,
        "peak_freq": peak_freq,
        "peak_phase": peak_phase,
        "peak_to_skirt": peak_to_skirt,
        "phase_coherence": coherence,
        "band_energy": band_energy,
        "snr": snr,
    }


# =============================================================================
# Statistical Analysis Functions
# =============================================================================


def compute_hedges_g(group1: np.ndarray, group2: np.ndarray) -> float:
    """
    Compute Hedges' g (bias-corrected Cohen's d).

    Args:
        group1: First group values
        group2: Second group values

    Returns:
        Hedges' g effect size
    """
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return np.nan

    mean1, mean2 = np.mean(group1), np.mean(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)

    # Pooled standard deviation
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return np.nan

    # Cohen's d
    cohens_d = (mean1 - mean2) / pooled_std

    # Bias correction factor (Hedges' g)
    correction = 1 - (3 / (4 * (n1 + n2) - 9))
    hedges_g = cohens_d * correction

    return float(hedges_g)


def bootstrap_hedges_g_ci(
    group1: np.ndarray,
    group2: np.ndarray,
    n_bootstrap: int = 500,
    confidence: float = 0.95,
    seed: int = 42,
) -> Tuple[float, float, float]:
    """
    Compute Hedges' g with bootstrap confidence interval.

    Args:
        group1: First group values
        group2: Second group values
        n_bootstrap: Number of bootstrap iterations
        confidence: Confidence level (e.g., 0.95 for 95% CI)
        seed: Random seed for reproducibility

    Returns:
        Tuple of (hedges_g, ci_low, ci_high)
    """
    rng = np.random.RandomState(seed)
    hedges_g = compute_hedges_g(group1, group2)

    if np.isnan(hedges_g):
        return np.nan, np.nan, np.nan

    bootstrap_gs = []
    for _ in range(n_bootstrap):
        bs_g1 = rng.choice(group1, len(group1), replace=True)
        bs_g2 = rng.choice(group2, len(group2), replace=True)
        bs_g = compute_hedges_g(bs_g1, bs_g2)
        if not np.isnan(bs_g):
            bootstrap_gs.append(bs_g)

    if len(bootstrap_gs) < MIN_BOOTSTRAP_SAMPLES:
        return hedges_g, np.nan, np.nan

    alpha = 1 - confidence
    ci_low = float(np.percentile(bootstrap_gs, 100 * alpha / 2))
    ci_high = float(np.percentile(bootstrap_gs, 100 * (1 - alpha / 2)))

    return hedges_g, ci_low, ci_high


def wilcoxon_ranksum_test(group1: np.ndarray, group2: np.ndarray) -> float:
    """
    Perform Wilcoxon rank-sum (Mann-Whitney U) test.

    Args:
        group1: First group values
        group2: Second group values

    Returns:
        Two-sided p-value
    """
    if len(group1) < 2 or len(group2) < 2:
        return 1.0

    _, p_value = stats.mannwhitneyu(group1, group2, alternative="two-sided")
    return float(p_value)


def benjamini_hochberg_correction(p_values: List[float]) -> List[float]:
    """
    Apply Benjamini-Hochberg FDR correction.

    Args:
        p_values: List of p-values

    Returns:
        FDR-adjusted p-values
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort by p-value with original indices
    sorted_pairs = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n

    # Calculate adjusted p-values
    min_so_far = 1.0
    for i in range(n - 1, -1, -1):
        orig_idx, p = sorted_pairs[i]
        adjusted_p = min(min_so_far, p * n / (i + 1))
        adjusted_p = min(adjusted_p, 1.0)
        adjusted[orig_idx] = adjusted_p
        min_so_far = adjusted_p

    return adjusted


# =============================================================================
# Data Loading Functions
# =============================================================================


def load_sequences_from_fasta(fasta_path: str, max_sequences: int = 2000) -> List[str]:
    """
    Load DNA sequences from FASTA file.

    Args:
        fasta_path: Path to FASTA file
        max_sequences: Maximum number of sequences to load

    Returns:
        List of sequences
    """
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        # Validate sequence (CRISPR guide length, ACGT only)
        if len(seq) == CRISPR_GUIDE_LENGTH and all(b in "ACGT" for b in seq):
            sequences.append(seq)
            if len(sequences) >= max_sequences:
                break
    return sequences


def group_by_gc_content(
    sequences: List[str],
    threshold: float = 0.5,
    seed: int = 42,
) -> Tuple[List[str], List[str]]:
    """
    Group sequences by GC content.

    Args:
        sequences: List of sequences
        threshold: GC content threshold (default 0.5 = 50%)
        seed: Random seed for subsampling

    Returns:
        Tuple of (high_gc_sequences, low_gc_sequences)
    """
    rng = random.Random(seed)

    high_gc = [s for s in sequences if compute_gc_content(s) > threshold]
    low_gc = [s for s in sequences if compute_gc_content(s) <= threshold]

    # Ensure balanced groups
    min_size = min(len(high_gc), len(low_gc))
    if min_size < 100:
        warnings.warn(f"Small group size: high_gc={len(high_gc)}, low_gc={len(low_gc)}")

    # Subsample to equal sizes (max 500 per group as per issue spec)
    target_size = min(min_size, 500)
    if len(high_gc) > target_size:
        high_gc = rng.sample(high_gc, target_size)
    if len(low_gc) > target_size:
        low_gc = rng.sample(low_gc, target_size)

    return high_gc, low_gc


# =============================================================================
# Main Analysis Pipeline
# =============================================================================


def process_sequence(seq: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a single sequence and extract features.

    Args:
        seq: DNA sequence
        params: Analysis parameters

    Returns:
        Dictionary of features
    """
    try:
        signal = encode_sequence(
            seq,
            at_lifetime=params["at_lifetime"],
            gc_lifetime=params["gc_lifetime"],
            helical_period=params["helical_period"],
            apply_helical=params["apply_helical"],
        )
        freqs, spectrum = czt_analysis(
            signal,
            f0=params["target_frequency"],
            band_width=params["band_width"],
            noise_band_width=params.get("noise_band_width", 0.03),
            M=params["czt_points"],
        )
        features = extract_features(
            freqs,
            spectrum,
            f0=params["target_frequency"],
            band_width=params["band_width"],
        )
        features["sequence"] = seq
        features["gc_content"] = compute_gc_content(seq)
        return features
    except (ValueError, IndexError, KeyError) as e:
        return {"error": str(e), "sequence": seq}


def run_validation_analysis(
    data_path: Optional[str] = None,
    output_dir: Optional[str] = None,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Run the full validation analysis pipeline.

    Args:
        data_path: Path to input FASTA file
        output_dir: Output directory for results
        seed: Random seed for reproducibility

    Returns:
        Dictionary of analysis results
    """
    # Set random seeds
    random.seed(seed)
    np.random.seed(seed)

    # Set default paths
    repo_root = Path(__file__).parent.parent
    if data_path is None:
        data_path = str(repo_root / "data" / "raw" / "brunello_parsed.fasta")
    if output_dir is None:
        output_dir = str(repo_root)

    output_path = Path(output_dir)
    reports_dir = output_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Analysis parameters
    params = {
        "dataset": Path(data_path).name,
        "seed": seed,
        "temperature_k": 310.15,  # 37°C physiological
        "na_concentration_m": 1.0,  # 1M Na⁺
        "at_lifetime": 1.0,  # ms (faster opening, 2 H-bonds)
        "gc_lifetime": 50.0,  # ms (slower opening, 3 H-bonds)
        "helical_period": 10.5,  # B-form DNA
        "apply_helical": True,
        "target_frequency": 1.0 / 10.5,  # bp⁻¹
        "band_width": 0.01,
        "noise_band_width": 0.03,  # Extended band for noise estimation
        "czt_points": 256,
        "gc_threshold": 0.5,  # 50%
        "max_sequences": 2000,
        "max_per_group": 500,
        "n_bootstrap": 500,
        "confidence_level": 0.95,
        "analysis_timestamp": datetime.now(timezone.utc).isoformat(),
    }

    print(f"Loading sequences from: {data_path}")

    # Load and group sequences
    max_seqs: int = params["max_sequences"]  # type: ignore[assignment]
    sequences = load_sequences_from_fasta(data_path, max_seqs)
    if len(sequences) < 100:
        raise ValueError(
            f"Insufficient sequences: {len(sequences)} (need at least 100)"
        )

    print(f"Loaded {len(sequences)} valid sequences")

    gc_threshold: float = params["gc_threshold"]  # type: ignore[assignment]
    high_gc, low_gc = group_by_gc_content(
        sequences,
        threshold=gc_threshold,
        seed=seed,
    )

    print(f"Groups: high_gc={len(high_gc)}, low_gc={len(low_gc)}")

    # Process all sequences
    all_features = []

    for seq in high_gc:
        features = process_sequence(seq, params)
        if "error" not in features:
            features["group"] = "high_gc"
            all_features.append(features)

    for seq in low_gc:
        features = process_sequence(seq, params)
        if "error" not in features:
            features["group"] = "low_gc"
            all_features.append(features)

    print(f"Processed {len(all_features)} sequences successfully")

    # Extract feature arrays by group
    high_gc_features = [f for f in all_features if f["group"] == "high_gc"]
    low_gc_features = [f for f in all_features if f["group"] == "low_gc"]

    # Key metrics to analyze
    metrics = ["peak_mag", "snr", "phase_coherence", "band_energy"]
    results: Dict[str, Any] = {"metrics": {}}

    raw_p_values: List[float] = []
    metric_results: List[Dict[str, Any]] = []

    n_bootstrap: int = params["n_bootstrap"]  # type: ignore[assignment]
    confidence_level: float = params["confidence_level"]  # type: ignore[assignment]

    for metric in metrics:
        high_vals = np.array([f[metric] for f in high_gc_features])
        low_vals = np.array([f[metric] for f in low_gc_features])

        # Hedges' g with bootstrap CI
        g, ci_low, ci_high = bootstrap_hedges_g_ci(
            high_vals,
            low_vals,
            n_bootstrap=n_bootstrap,
            confidence=confidence_level,
            seed=seed,
        )

        # Wilcoxon test
        p_value = wilcoxon_ranksum_test(high_vals, low_vals)
        raw_p_values.append(p_value)

        # Normality test (Kolmogorov-Smirnov)
        combined = np.concatenate([high_vals, low_vals])
        ks_stat, ks_p = stats.kstest(
            combined, "norm", args=(np.mean(combined), np.std(combined))
        )

        metric_results.append(
            {
                "metric": metric,
                "hedges_g": g,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "wilcoxon_p": p_value,
                "high_gc_mean": float(np.mean(high_vals)),
                "high_gc_std": float(np.std(high_vals)),
                "low_gc_mean": float(np.mean(low_vals)),
                "low_gc_std": float(np.std(low_vals)),
                "ks_stat": float(ks_stat),
                "ks_p": float(ks_p),
                "normal_distribution": ks_p > 0.05,
            }
        )

    # Apply FDR correction
    adjusted_p_values = benjamini_hochberg_correction(raw_p_values)
    for i, res in enumerate(metric_results):
        res["wilcoxon_p_fdr"] = adjusted_p_values[i]
        results["metrics"][res["metric"]] = res

    # Summary statistics
    results["summary"] = {
        "total_sequences": len(all_features),
        "high_gc_count": len(high_gc_features),
        "low_gc_count": len(low_gc_features),
        "primary_metric": "peak_mag",
        "primary_hedges_g": results["metrics"]["peak_mag"]["hedges_g"],
        "primary_ci": [
            results["metrics"]["peak_mag"]["ci_low"],
            results["metrics"]["peak_mag"]["ci_high"],
        ],
        "primary_p_fdr": results["metrics"]["peak_mag"]["wilcoxon_p_fdr"],
    }

    # Save params.yaml
    params_path = output_path / "params.yaml"
    with open(params_path, "w") as f:
        f.write("# DNA Breathing Dynamics Validation Analysis Parameters\n")
        f.write(f"# Generated: {params['analysis_timestamp']}\n\n")
        for key, value in params.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.6f}\n")
            else:
                f.write(f"{key}: {value}\n")
    print(f"Saved parameters to: {params_path}")

    # Save features.csv
    features_path = output_path / "features.csv"
    fieldnames = [
        "sequence",
        "group",
        "gc_content",
        "peak_mag",
        "peak_freq",
        "peak_phase",
        "peak_to_skirt",
        "phase_coherence",
        "band_energy",
        "snr",
    ]
    with open(features_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for feat in all_features:
            row = {k: feat.get(k, "") for k in fieldnames}
            writer.writerow(row)
    print(f"Saved features to: {features_path}")

    # Generate validation_stats.md report
    report_path = reports_dir / "validation_stats.md"
    _generate_report(report_path, params, results, metric_results)
    print(f"Saved report to: {report_path}")

    return results


def _generate_report(
    report_path: Path,
    params: Dict[str, Any],
    results: Dict[str, Any],
    metric_results: List[Dict[str, Any]],
) -> None:
    """Generate the validation statistics markdown report."""
    with open(report_path, "w") as f:
        f.write("# DNA Breathing Dynamics Validation Statistics\n\n")
        f.write(f"**Generated:** {params['analysis_timestamp']}\n\n")

        f.write("## Executive Summary\n\n")
        summary = results["summary"]
        f.write(f"- **Dataset:** {params['dataset']}\n")
        f.write(f"- **Total Sequences:** {summary['total_sequences']}\n")
        f.write(f"- **High GC Group (>50%):** {summary['high_gc_count']} sequences\n")
        f.write(f"- **Low GC Group (≤50%):** {summary['low_gc_count']} sequences\n\n")

        f.write("### Primary Finding\n\n")
        g = summary["primary_hedges_g"]
        ci_low, ci_high = summary["primary_ci"]
        p_fdr = summary["primary_p_fdr"]
        f.write("**Spectral Power at 1/10.5 bp⁻¹ (helical frequency):**\n\n")
        f.write(f"- Hedges' g = {g:.4f} [95% CI: {ci_low:.4f}, {ci_high:.4f}]\n")
        f.write(f"- Wilcoxon p-value (FDR-corrected) = {p_fdr:.4e}\n\n")

        # Interpretation
        if abs(g) >= 0.8:
            effect_size_interp = "large"
        elif abs(g) >= 0.5:
            effect_size_interp = "medium"
        elif abs(g) >= 0.2:
            effect_size_interp = "small"
        else:
            effect_size_interp = "negligible"

        f.write(f"**Interpretation:** A {effect_size_interp} effect size ")
        if g > 0:
            f.write(
                "indicates that high-GC sequences show stronger spectral resonance "
            )
        else:
            f.write("indicates that low-GC sequences show stronger spectral resonance ")
        f.write("at the B-form DNA helical frequency.\n\n")

        f.write("## Detailed Metrics\n\n")
        header = "| Metric | Hedges' g | 95% CI | p-value (FDR) "
        header += "| High GC Mean ± SD | Low GC Mean ± SD |\n"
        f.write(header)
        separator = "|--------|-----------|--------|---------------"
        separator += "|-------------------|------------------|\n"
        f.write(separator)

        for res in metric_results:
            metric = res["metric"]
            g = res["hedges_g"]
            ci_low = res["ci_low"]
            ci_high = res["ci_high"]
            p_fdr = res["wilcoxon_p_fdr"]
            high_mean = res["high_gc_mean"]
            high_std = res["high_gc_std"]
            low_mean = res["low_gc_mean"]
            low_std = res["low_gc_std"]

            f.write(
                f"| {metric} | {g:.4f} | [{ci_low:.4f}, {ci_high:.4f}] | "
                f"{p_fdr:.4e} | {high_mean:.4f} ± {high_std:.4f} | "
                f"{low_mean:.4f} ± {low_std:.4f} |\n"
            )

        f.write("\n## Normality Assessment\n\n")
        f.write("| Metric | K-S Statistic | K-S p-value | Normal? |\n")
        f.write("|--------|---------------|-------------|--------|\n")
        for res in metric_results:
            normal = "Yes" if res["normal_distribution"] else "No"
            f.write(
                f"| {res['metric']} | {res['ks_stat']:.4f} | "
                f"{res['ks_p']:.4e} | {normal} |\n"
            )

        f.write(
            "\n*Note: Non-parametric Wilcoxon test used regardless of normality.*\n\n"
        )

        f.write("## Analysis Parameters\n\n")
        f.write("```yaml\n")
        for key, value in params.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.6f}\n")
            else:
                f.write(f"{key}: {value}\n")
        f.write("```\n\n")

        f.write("## Reproducibility\n\n")
        f.write(f"- **Random Seed:** {params['seed']}\n")
        f.write(f"- **Bootstrap Iterations:** {params['n_bootstrap']}\n")
        f.write(f"- **Confidence Level:** {params['confidence_level'] * 100:.0f}%\n")
        f.write(
            f"- **Temperature:** {params['temperature_k']} K (37°C physiological)\n\n"
        )

        f.write("## Methodology\n\n")
        f.write("1. **Sequence Encoding:** Complex-valued signal with real part = ")
        f.write("breathing lifetimes (AT=1.0 ms, GC=50.0 ms) and imaginary part = ")
        f.write("SantaLucia nearest-neighbor thermodynamic stability.\n\n")
        f.write("2. **Spectral Analysis:** Chirp Z-Transform (CZT) focused on ")
        f.write("1/10.5 bp⁻¹ helical frequency with helical phase modulation.\n\n")
        f.write("3. **Effect Size:** Hedges' g (bias-corrected Cohen's d) with ")
        f.write("500-iteration bootstrap for 95% confidence interval.\n\n")
        f.write("4. **Significance Testing:** Two-sided Wilcoxon rank-sum test ")
        f.write("with Benjamini-Hochberg FDR correction for multiple comparisons.\n")


if __name__ == "__main__":
    try:
        results = run_validation_analysis()
        print("\n" + "=" * 60)
        print("VALIDATION ANALYSIS COMPLETE")
        print("=" * 60)
        summary = results["summary"]
        g = summary["primary_hedges_g"]
        ci_low, ci_high = summary["primary_ci"]
        p = summary["primary_p_fdr"]
        print(f"Peak Magnitude: g = {g:.4f} [95% CI: {ci_low:.4f}, {ci_high:.4f}]")
        print(f"Wilcoxon p-value (FDR): {p:.4e}")
        print(f"Total sequences analyzed: {summary['total_sequences']}")
        print("=" * 60)
    except (ValueError, FileNotFoundError, IOError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
