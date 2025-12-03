#!/usr/bin/env python3
"""
Prototype Alpha: CZT-based DNA Breathing Dynamics Feature Extractor.

This script implements a minimal, self-contained pipeline for extracting
CZT-based spectral features from DNA sequences using breathing lifetime
and thermodynamic stability encoding.

Features:
- Biophysics encoding: breathing lifetimes (τ) + thermodynamic stability (ΔG)
- CZT spectral analysis focused on 10.5 bp/turn helical frequency
- Feature extraction: mag_peak, band_energy, phase_peak, freq_eff
- GC-content baseline comparison
- Random encoding control (shuffled τ/ΔG)
- Minimal stats: ROC AUC, Cohen's d, t-test/Mann-Whitney U

Usage:
    python prototype_alpha.py --fasta sequences.fasta --labels labels.csv
"""

import argparse
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.signal import czt
from scipy.stats import mannwhitneyu, ttest_ind
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None  # type: ignore[assignment]

# =============================================================================
# Configuration: Deterministic seed for reproducibility
# =============================================================================
RNG_SEED = 137
np.random.seed(RNG_SEED)

# =============================================================================
# Biophysics Parameters (editable tables)
# =============================================================================

# Breathing lifetimes [ms]: AT pairs are more flexible, GC pairs are more stable
BREATH_LIFETIMES_MS: Dict[str, float] = {
    "A": 5.0,
    "T": 5.0,
    "G": 25.0,
    "C": 25.0,
}

# Thermodynamic stability: per-bp ΔG in kcal/mol (placeholder values)
# Can be replaced by SantaLucia nearest-neighbor parameters or EPBD outputs
DELTA_G: Dict[str, float] = {
    "A": -1.0,
    "T": -1.0,
    "G": -2.0,
    "C": -2.0,
}


def encode_sequence(seq: str) -> Optional[np.ndarray]:
    """
    Encode DNA sequence into complex signal: x_i = τ_i + 1j * ΔG_i.

    Args:
        seq: DNA sequence string (ATGC characters only)

    Returns:
        Complex numpy array or None if sequence has no valid bases
    """
    taus: List[float] = []
    dgs: List[float] = []
    for b in seq.upper():
        if b not in BREATH_LIFETIMES_MS:
            continue
        taus.append(BREATH_LIFETIMES_MS[b])
        dgs.append(DELTA_G[b])
    if not taus:
        return None
    x: np.ndarray = np.array(taus) + 1j * np.array(dgs)
    return x


def gc_content(seq: str) -> float:
    """
    Calculate GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC fraction (0.0 to 1.0)
    """
    s = seq.upper()
    gc = sum(b in "GC" for b in s)
    atgc = sum(b in "ATGC" for b in s)
    return gc / atgc if atgc > 0 else 0.0


def czt_band_features(
    x: np.ndarray,
    f0: float = 1 / 10.5,
    rel_bw: float = 0.1,
    m: int = 64,
) -> np.ndarray:
    """
    Apply CZT around target helical frequency and extract spectral features.

    Args:
        x: Complex signal array
        f0: Target spatial frequency (cycles/bp), default 1/10.5 for DNA helix
        rel_bw: Relative bandwidth (±10% of f0 by default)
        m: Number of frequency bins in CZT

    Returns:
        Feature array: [mag_peak, band_energy, phase_peak, freq_eff]
    """
    # Spatial frequency in cycles/bp mapped to digital radian frequency
    w0 = 2 * np.pi * f0
    # Center step and band
    w_start = w0 * (1 - rel_bw)
    w_end = w0 * (1 + rel_bw)
    # CZT: chirp z-transform mapping from w_start to w_end
    A = np.exp(1j * w_start)
    W = np.exp(1j * (w_end - w_start) / (m - 1))
    X = czt(x, m=m, w=W, a=A)

    mags = np.abs(X)
    idx = np.argmax(mags)
    mag_peak = mags[idx]
    band_energy = np.sum(mags**2)
    phase_peak = np.angle(X[idx])
    freq_eff = f0 * (1 - rel_bw + 2 * rel_bw * idx / (m - 1))

    return np.array([mag_peak, band_energy, phase_peak, freq_eff])


def random_encoding_control(x: np.ndarray) -> np.ndarray:
    """
    Create random encoding control by shuffling real/imag parts separately.

    This preserves marginal distributions but destroys position-dependent
    biophysical structure.

    Args:
        x: Complex signal array

    Returns:
        Shuffled complex signal array
    """
    r = x.real.copy()
    i = x.imag.copy()
    np.random.shuffle(r)
    np.random.shuffle(i)
    result: np.ndarray = r + 1j * i
    return result


def fasta_to_df(fasta_path: str, labels_csv: str) -> pd.DataFrame:
    """
    Load sequences from FASTA and labels from CSV.

    Args:
        fasta_path: Path to FASTA file with sequences
        labels_csv: Path to CSV with columns: id, label

    Returns:
        DataFrame with columns: id, seq, label
    """
    if SeqIO is None:
        raise ImportError(
            "BioPython is required for FASTA parsing. "
            "Install with: pip install biopython"
        )

    labels = pd.read_csv(labels_csv)
    label_map = dict(zip(labels["id"], labels["label"]))
    rows: List[Dict[str, object]] = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq_id = rec.id
        if seq_id not in label_map:
            continue
        rows.append({"id": seq_id, "seq": str(rec.seq), "label": label_map[seq_id]})
    return pd.DataFrame(rows)


def build_feature_table(
    df: pd.DataFrame,
) -> List[Dict[str, object]]:
    """
    Build feature vectors for all sequences.

    Features include:
    - CZT features: mag_peak, band_energy, phase_peak, freq_eff
    - Simple stats: mean(τ), var(τ), mean(ΔG), var(ΔG), gc
    - Random control CZT features

    Args:
        df: DataFrame with columns: id, seq, label

    Returns:
        List of feature dictionaries
    """
    records: List[Dict[str, object]] = []
    for _, row in df.iterrows():
        seq = row["seq"]
        x = encode_sequence(seq)
        if x is None:
            continue
        czt_feat = czt_band_features(x)
        gc = gc_content(seq)
        tau = x.real
        dg = x.imag
        stats = np.array([tau.mean(), tau.var(), dg.mean(), dg.var(), gc])
        # Random-encoding CZT control
        xr = random_encoding_control(x)
        czt_rand = czt_band_features(xr)
        feat = np.concatenate([czt_feat, stats, czt_rand])
        records.append(
            {
                "id": row["id"],
                "label": row["label"],
                "feat": feat,
                "mag_peak": czt_feat[0],
                "mag_peak_rand": czt_rand[0],
                "gc": gc,
            }
        )
    return records


def stack_features(records: List[Dict[str, object]]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Stack feature vectors into arrays for sklearn.

    Args:
        records: List of feature dictionaries

    Returns:
        Tuple of (X, y) arrays
    """
    feat_list = [r["feat"] for r in records]
    X: np.ndarray = np.stack(feat_list, axis=0)  # type: ignore[arg-type]
    y: np.ndarray = np.array([r["label"] for r in records])
    return X, y


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    """
    Calculate Cohen's d effect size.

    Args:
        a: First group values
        b: Second group values

    Returns:
        Cohen's d value
    """
    a = np.asarray(a)
    b = np.asarray(b)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 0.0
    va, vb = float(a.var(ddof=1)), float(b.var(ddof=1))
    sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    return float((a.mean() - b.mean()) / sp) if sp > 0 else 0.0


def run_experiment(fasta_path: str, labels_csv: str, test_size: float = 0.3) -> None:
    """
    Run the full CZT breathing dynamics experiment.

    This function:
    1. Loads sequences and labels
    2. Extracts CZT and baseline features
    3. Trains logistic regression classifiers
    4. Reports ROC AUC, Cohen's d, and statistical tests

    Args:
        fasta_path: Path to FASTA file
        labels_csv: Path to labels CSV
        test_size: Test set proportion (default 0.3)
    """
    print("=" * 60)
    print("Prototype Alpha: CZT Breathing Dynamics Experiment")
    print("=" * 60)
    print("\nConfiguration:")
    print(f"  RNG Seed: {RNG_SEED}")
    print(f"  FASTA: {fasta_path}")
    print(f"  Labels: {labels_csv}")
    print(f"  Test Size: {test_size}")
    print(f"  Target Frequency: 1/10.5 = {1/10.5:.6f} cycles/bp")
    print()

    # Load data
    df = fasta_to_df(fasta_path, labels_csv)
    records = build_feature_table(df)
    X, y = stack_features(records)

    if len(records) == 0:
        print("ERROR: No valid sequences found.")
        return

    print("Data Summary:")
    print(f"  N = {len(records)}")
    print(f"  Label 0: {sum(y == 0)}")
    print(f"  Label 1: {sum(y == 1)}")
    print()

    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=RNG_SEED, stratify=y
    )

    # (a) CZT breathing features model
    clf = LogisticRegression(max_iter=1000)
    clf.fit(X_train, y_train)
    y_prob = clf.predict_proba(X_test)[:, 1]
    auc_czt = roc_auc_score(y_test, y_prob)

    # (b) GC-only baseline model
    gc_vals = np.array([r["gc"] for r in records]).reshape(-1, 1)
    Xg_tr, Xg_te, yg_tr, yg_te = train_test_split(
        gc_vals, y, test_size=test_size, random_state=RNG_SEED, stratify=y
    )
    clf_gc = LogisticRegression(max_iter=1000)
    clf_gc.fit(Xg_tr, yg_tr)
    y_prob_gc = clf_gc.predict_proba(Xg_te)[:, 1]
    auc_gc = roc_auc_score(yg_te, y_prob_gc)

    # Scalar feature statistics
    mag_peak = np.array([r["mag_peak"] for r in records])
    mag_peak_rand = np.array([r["mag_peak_rand"] for r in records])
    lbl = np.array([r["label"] for r in records])

    m1 = mag_peak[lbl == 1]
    m0 = mag_peak[lbl == 0]
    d = cohen_d(m1, m0)
    tstat, pval = ttest_ind(m1, m0, equal_var=False)
    ustat, pval_mw = mannwhitneyu(m1, m0, alternative="two-sided")

    # Control: random encoding
    mr1 = mag_peak_rand[lbl == 1]
    mr0 = mag_peak_rand[lbl == 0]
    d_rand = cohen_d(mr1, mr0)
    _, pval_rand = ttest_ind(mr1, mr0, equal_var=False)

    # Report results
    print("=" * 60)
    print("Results")
    print("=" * 60)
    print("\nROC AUC Comparison:")
    print(f"  CZT Breathing Features: {auc_czt:.4f}")
    print(f"  GC-Only Baseline:       {auc_gc:.4f}")
    print()
    print("Effect Size (mag_peak, true encoding):")
    print(f"  Cohen's d = {d:.4f}")
    print(f"  t-test p  = {pval:.6f}")
    print(f"  Mann-Whitney U p = {pval_mw:.6f}")
    print()
    print("Control (mag_peak, random encoding):")
    print(f"  Cohen's d = {d_rand:.4f}")
    print(f"  t-test p  = {pval_rand:.6f}")
    print()
    print("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prototype Alpha: CZT Breathing Dynamics Feature Extractor"
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to FASTA file with DNA sequences",
    )
    parser.add_argument(
        "--labels",
        required=True,
        help="Path to CSV file with columns: id, label",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.3,
        help="Test set proportion (default: 0.3)",
    )
    args = parser.parse_args()
    run_experiment(args.fasta, args.labels, args.test_size)
