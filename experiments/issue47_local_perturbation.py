#!/usr/bin/env python3
"""
Issue #47 local perturbation validation for DNA breathing dynamics encoding.

Generates matched wild‑type vs single GC‑altering mutant guide pairs (20 nt),
computes helical-band features, and runs paired statistical tests per metric.

Outputs:
  - experiments/issue47_pairs.fasta : WT/mutant FASTA (paired order)
  - experiments/issue47_stats.csv   : per‑metric statistics with BH-corrected p-values
  - experiments/issue47_diffplots.png : difference distributions per metric
"""

import csv
import itertools
import math
import random
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Import gist utilities
GIST_SRC = Path(__file__).parent.parent / "gists/breathing/czt_feature_extractor/src"
sys.path.insert(0, str(GIST_SRC))
from dna_breathing_gist import encode_sequence, czt_analysis, extract_features  # type: ignore

DATA_PATH = (
    Path(__file__).parent.parent
    / "gists/breathing/czt_feature_extractor/data/processed/brunello_gc1000.fasta"
)
OUT_FASTA = Path(__file__).parent / "issue47_pairs.fasta"
OUT_STATS = Path(__file__).parent / "issue47_stats.csv"
OUT_PLOT = Path(__file__).parent / "issue47_diffplots.png"

RNG_SEED = 47
NUM_PAIRS = 50  # 25 GC-increase, 25 GC-decrease
HELICAL_FREQ = 1 / 10.5
METRICS = ["peak_mag", "snr", "phase_coherence", "band_energy"]


def load_sequences(path: Path) -> List[str]:
    """Load FASTA sequences as uppercase strings."""
    seqs = []
    with path.open() as f:
        current = ""
        for line in f:
            if line.startswith(">"):
                if current:
                    seqs.append(current.upper())
                current = ""
            else:
                current += line.strip()
        if current:
            seqs.append(current.upper())
    return seqs


def mutate_one_base(seq: str, to_gc: bool) -> Tuple[str, int, str]:
    """
    Mutate a single position to increase (to_gc=True) or decrease GC content.
    Returns (mutated_seq, pos, new_base). Raises if no suitable position.
    """
    candidates = []
    if to_gc:
        for i, b in enumerate(seq):
            if b in ("A", "T"):
                candidates.append(i)
    else:
        for i, b in enumerate(seq):
            if b in ("G", "C"):
                candidates.append(i)
    if not candidates:
        raise ValueError("No eligible position for requested GC change")

    pos = candidates[0]  # deterministic choice
    orig = seq[pos]
    if to_gc:
        new_b = "G" if orig == "A" else "C"
    else:
        new_b = "A" if orig == "G" else "T"

    mutated = seq[:pos] + new_b + seq[pos + 1 :]
    return mutated, pos, new_b


def build_pairs(seqs: List[str]) -> List[Tuple[str, str, str]]:
    """
    Build paired WT/mutant sequences.
    Returns list of (wt, mut, direction) with direction in {"increase_gc","decrease_gc"}.
    """
    random.seed(RNG_SEED)
    increase_pairs = []
    decrease_pairs = []

    for seq in seqs:
        if len(increase_pairs) < NUM_PAIRS // 2:
            try:
                mut, pos, nb = mutate_one_base(seq, to_gc=True)
                if mut != seq:
                    increase_pairs.append((seq, mut, "increase_gc"))
            except ValueError:
                pass
        if len(decrease_pairs) < NUM_PAIRS // 2:
            try:
                mut, pos, nb = mutate_one_base(seq, to_gc=False)
                if mut != seq:
                    decrease_pairs.append((seq, mut, "decrease_gc"))
            except ValueError:
                pass
        if len(increase_pairs) + len(decrease_pairs) >= NUM_PAIRS:
            break

    if len(increase_pairs) + len(decrease_pairs) < NUM_PAIRS:
        raise RuntimeError("Insufficient sequences to construct required pairs")

    pairs = increase_pairs[: NUM_PAIRS // 2] + decrease_pairs[: NUM_PAIRS // 2]
    assert len(pairs) == NUM_PAIRS
    return pairs


def compute_metric_vector(seq: str) -> dict:
    signal = encode_sequence(seq)
    freqs, spectrum = czt_analysis(signal)
    feats = extract_features(freqs, spectrum, HELICAL_FREQ, 0.01)
    return feats


def paired_stats(diffs: np.ndarray, alpha=0.05, boot_iters=5000) -> dict:
    """Compute paired t, Wilcoxon, Cohen's d, Hedges' g, and bootstrap CI."""
    n = len(diffs)
    mean_diff = float(np.mean(diffs))
    sd_diff = float(np.std(diffs, ddof=1)) if n > 1 else float("nan")

    t_stat, t_p = stats.ttest_rel(diffs, np.zeros_like(diffs))
    # Wilcoxon signed-rank; use Pratt to handle zeros without dropping them
    try:
        w_stat, w_p = stats.wilcoxon(diffs, zero_method="pratt", alternative="two-sided")
    except ValueError:
        w_stat, w_p = math.nan, math.nan

    d = mean_diff / sd_diff if sd_diff and not math.isnan(sd_diff) else math.nan
    j = 1 - 3 / (4 * n - 1) if n > 1 else 1.0
    g = d * j if not math.isnan(d) else math.nan

    rng = np.random.default_rng(RNG_SEED)
    boot_ds = []
    for _ in range(boot_iters):
        sample = rng.choice(diffs, size=n, replace=True)
        sd = sample.std(ddof=1)
        boot_ds.append(sample.mean() / sd if sd > 0 else math.nan)
    boot_ds = np.array(boot_ds)
    ci_low, ci_high = np.nanpercentile(boot_ds, [2.5, 97.5])

    return {
        "n": n,
        "mean_diff": mean_diff,
        "sd_diff": sd_diff,
        "t_stat": t_stat,
        "t_p_raw": t_p,
        "w_stat": w_stat,
        "w_p_raw": w_p,
        "cohens_d": d,
        "hedges_g": g,
        "ci_low": ci_low,
        "ci_high": ci_high,
    }


def main():
    seqs = load_sequences(DATA_PATH)
    pairs = build_pairs(seqs)

    # Persist FASTA for auditability
    with OUT_FASTA.open("w") as f:
        for idx, (wt, mut, direction) in enumerate(pairs, 1):
            f.write(f">{idx}|wt|{direction}\n{wt}\n")
            f.write(f">{idx}|mut|{direction}\n{mut}\n")

    records = []
    for idx, (wt, mut, direction) in enumerate(pairs, 1):
        wt_feats = compute_metric_vector(wt)
        mut_feats = compute_metric_vector(mut)
        rec = {"pair_id": idx, "direction": direction}
        for m in METRICS:
            rec[f"wt_{m}"] = wt_feats[m]
            rec[f"mut_{m}"] = mut_feats[m]
            rec[f"diff_{m}"] = mut_feats[m] - wt_feats[m]
        records.append(rec)

    df = pd.DataFrame(records)

    stats_rows = []
    # Collect raw p-values for BH correction
    t_pvals = []
    w_pvals = []
    for m in METRICS:
        diffs = df[f"diff_{m}"].to_numpy()
        s = paired_stats(diffs)
        s["metric"] = m
        stats_rows.append(s)
        t_pvals.append(s["t_p_raw"])
        w_pvals.append(s["w_p_raw"])

    # Benjamini–Hochberg across metrics (separately for t and Wilcoxon)
    t_adj = multipletests(t_pvals, method="fdr_bh")[1]
    w_adj = multipletests(w_pvals, method="fdr_bh")[1]
    for i, s in enumerate(stats_rows):
        s["t_p_adj"] = t_adj[i]
        s["w_p_adj"] = w_adj[i]
        s["success_flag"] = (
            abs(s["cohens_d"]) >= 4
            and s["t_p_adj"] < 0.001
            and not (s["ci_low"] <= 0 <= s["ci_high"])
        )

    # Save stats table
    fieldnames = [
        "metric",
        "n",
        "mean_diff",
        "sd_diff",
        "t_stat",
        "t_p_raw",
        "t_p_adj",
        "w_stat",
        "w_p_raw",
        "w_p_adj",
        "cohens_d",
        "hedges_g",
        "ci_low",
        "ci_high",
        "success_flag",
    ]
    with OUT_STATS.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in stats_rows:
            writer.writerow({k: s.get(k) for k in fieldnames})

    # Plot difference distributions
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.ravel()
    for ax, m in zip(axes, METRICS):
        diffs = df[f"diff_{m}"].to_numpy()
        ax.hist(diffs, bins=15, color="#3b5b92", alpha=0.8)
        ax.axvline(0, color="k", linestyle="--", linewidth=1)
        ax.set_title(f"{m} (mut - wt)")
        ax.set_xlabel("Difference")
        ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=200)

    successes = [s["success_flag"] for s in stats_rows]
    any_success = any(successes)
    print("Completed Issue #47 local perturbation validation.")
    for s in stats_rows:
        print(
            f"{s['metric']}: d={s['cohens_d']:.3f} "
            f"[{s['ci_low']:.3f}, {s['ci_high']:.3f}], "
            f"t_p_adj={s['t_p_adj']:.3e}, success={s['success_flag']}"
        )
    if not any_success:
        print("Success criterion not met (|d|>=4 with FDR-adjusted p<0.001 and CI excluding 0).")


if __name__ == "__main__":
    main()
