#!/usr/bin/env python3
# DNA Breathing Dynamics Encoding Gist: CZT-based Feature Extractor. See README.md for details.

import argparse
import sys
import csv
import numpy as np
from typing import List, Dict, Tuple
from pathlib import Path
import random
import collections

try:
    from scipy import stats
    from scipy.signal import czt
except ImportError as e:
    print(f"Error: SciPy >=1.12 required: {e}", file=sys.stderr)
    sys.exit(1)


NEAREST_NEIGHBOR_DG = {
    # SantaLucia 1998 Unified parameters (kcal/mol)
    "AA": -1.00, "TT": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45, "TG": -1.45,
    "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28,
    "GA": -1.30, "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84, "CC": -1.84,
}


def encode_sequence(
    seq,
    at_lifetime=5.0,
    gc_lifetime=25.0,
    helical_period=10.5,
    apply_helical=True,
    iupac_policy="average",
):
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        raise ValueError("Empty sequence")

    valid_bases = set("ATGC")
    # Simplified IUPAC handling for brevity in gist
    iupac_codes = {
        "N": "ATGC", "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
        "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    }

    complex_signal = np.zeros(n, dtype=complex)

    for i in range(n):
        base = seq[i]
        
        # 1. Real Part: Lifetime
        if base in iupac_codes:
            if iupac_policy == "mask":
                real_part = 0.0
            else:
                bases = list(iupac_codes[base])
                real_vals = [at_lifetime if b in "AT" else gc_lifetime for b in bases]
                real_part = np.mean(real_vals)
        elif base in valid_bases:
            real_part = at_lifetime if base in "AT" else gc_lifetime
        else:
             raise ValueError(f"Invalid base at pos {i}: {base}")

        # 2. Imaginary Part: Stability (Dinucleotide)
        # Assign DG(seq[i], seq[i+1]) to index i. Last base gets previous val or 0.
        if i < n - 1:
            dinuc = seq[i : i + 2]
            # Handle IUPAC in dinuc
            if any(b not in valid_bases for b in dinuc):
                if iupac_policy == "mask":
                    dg = 0.0
                else:
                    # Average over all expansions
                    b1s = list(iupac_codes.get(dinuc[0], dinuc[0]))
                    b2s = list(iupac_codes.get(dinuc[1], dinuc[1]))
                    dgs = []
                    for b1 in b1s:
                        for b2 in b2s:
                            if b1+b2 in NEAREST_NEIGHBOR_DG:
                                dgs.append(NEAREST_NEIGHBOR_DG[b1+b2])
                    dg = np.mean(dgs) if dgs else -1.5 # fallback
            else:
                dg = NEAREST_NEIGHBOR_DG.get(dinuc, -1.5)
        else:
            # Boundary condition: repeat last valid DG or neutral
            dg = -1.5 

        # Normalize DG to -1..1 or similar range for unit circle
        # Range is approx -0.58 to -2.24. 
        min_dg, max_dg = -2.5, -0.5 
        norm_imag = (dg - min_dg) / (max_dg - min_dg) # 0..1 approx
        
        complex_signal[i] = complex(real_part, norm_imag)

    if apply_helical:
        phases = np.exp(1j * 2 * np.pi * np.arange(n) / helical_period)
        complex_signal *= phases

    return complex_signal


def czt_analysis(
    signal, f0=1 / 10.5, band_width=0.01, M=256, apply_taper=False, taper_k=0.3
):
    n = len(signal)
    fs = 1.0

    if apply_taper:
        phi = (1 + np.sqrt(5)) / 2
        mod21 = np.arange(n) % 21 / 21.0
        taper = phi * (mod21**taper_k)
        signal = signal * taper

    f_low = f0 - band_width
    f_high = f0 + band_width
    a = np.exp(-2j * np.pi * f_low / fs)
    w = np.exp(2j * np.pi * (f_high - f_low) / (M * fs))

    spectrum = czt(signal, a=a, w=w, m=M)
    freqs = np.linspace(f_low, f_high, M)

    return freqs, spectrum


def extract_features(freqs, spectrum, f0, band_width, guard_band=0.005):
    mags = np.abs(spectrum)
    phases = np.angle(spectrum)

    in_band = (freqs >= f0 - band_width) & (freqs <= f0 + band_width)
    band_mags = mags[in_band]
    band_phases = phases[in_band]

    if len(band_mags) == 0:
        raise ValueError("No points in band")

    peak_local_idx = np.argmax(band_mags)
    peak_idx = np.where(in_band)[0][peak_local_idx]
    peak_mag = band_mags[peak_local_idx]
    peak_freq = freqs[peak_idx]
    peak_phase = band_phases[peak_local_idx]

    edge_low = f0 - band_width - guard_band
    edge_high = f0 + band_width + guard_band
    off_mask = (freqs < edge_low) | (freqs > edge_high)
    skirt_mags = mags[off_mask]
    skirt_mean = np.mean(skirt_mags) if len(skirt_mags) > 0 else 1.0
    peak_to_skirt = peak_mag / skirt_mean if skirt_mean > 0 else float("inf")

    coherence = np.abs(np.mean(np.exp(1j * band_phases)))

    band_energy = np.sum(band_mags**2)

    band_mean = np.mean(band_mags)
    snr = band_mean / skirt_mean if skirt_mean > 0 else float("inf")

    edge_threshold = 0.9 * band_width
    qc_flag = 1.0 if abs(peak_freq - f0) > edge_threshold else 0.0

    return {
        "peak_mag": peak_mag,
        "peak_idx": peak_idx,
        "peak_freq": peak_freq,
        "peak_phase": peak_phase,
        "peak_to_skirt": peak_to_skirt,
        "phase_coherence": coherence,
        "band_energy": band_energy,
        "snr": snr,
        "qc_flag": qc_flag,
    }


def generate_dinuc_shuffles(seq, num_shuffles=10, seed=42):
    """
    Generate dinucleotide-preserving shuffles using a randomized Eulerian path 
    approach with rejection sampling. Efficient for short sequences (e.g., CRISPR).
    """
    random.seed(seed)
    n = len(seq)
    if n < 2:
        return [seq] * num_shuffles

    # Build adjacency graph of transitions
    edges = collections.defaultdict(list)
    for i in range(n - 1):
        edges[seq[i]].append(seq[i+1])
    
    start_node = seq[0]
    shuffles = []
    
    max_retries = 100
    
    for _ in range(num_shuffles):
        for attempt in range(max_retries):
            # Deep copy edges to consume
            current_edges = {k: list(v) for k, v in edges.items()}
            # Shuffle the order of outgoing edges
            for k in current_edges:
                random.shuffle(current_edges[k])
            
            # Walk the graph
            curr = start_node
            path = [curr]
            while current_edges[curr]:
                next_node = current_edges[curr].pop()
                path.append(next_node)
                curr = next_node
            
            # Check if we visited all edges (length match)
            if len(path) == n:
                shuffles.append("".join(path))
                break
        else:
            # Fallback if rejection sampling fails (rare for short seqs)
            shuffles.append(seq)
            
    return shuffles


def phase_scramble(spectrum, seed=42):
    np.random.seed(seed)
    mags = np.abs(spectrum)
    random_phases = np.random.uniform(0, 2 * np.pi, len(spectrum))
    return mags * np.exp(1j * random_phases)


def compute_stats(features_list, groups=None, num_bootstrap=100, num_perm=20, seed=42):
    np.random.seed(seed)

    if groups is None or len(set(groups)) < 2:
        return {"cohens_d": 0.0, "ci_low": 0.0, "ci_high": 0.0, "p_perm": 1.0}

    unique_groups = list(set(groups))
    if len(unique_groups) < 2:
        return {"cohens_d": 0.0, "ci_low": 0.0, "ci_high": 0.0, "p_perm": 1.0}

    group1 = [
        f["peak_mag"] for f, g in zip(features_list, groups) if g == unique_groups[0]
    ]
    group2 = [
        f["peak_mag"] for f, g in zip(features_list, groups) if g == unique_groups[1]
    ]

    if len(group1) < 2 or len(group2) < 2:
        return {"cohens_d": 0.0, "ci_low": 0.0, "ci_high": 0.0, "p_perm": 1.0}

    # Cohen's d
    mean1, mean2 = np.mean(group1), np.mean(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt((var1 + var2) / 2)
    cohens_d = (mean1 - mean2) / pooled_std if pooled_std > 0 else 0.0

    # Bootstrap CI
    def bootstrap_d(bs_samples):
        bs_d = []
        for _ in range(bs_samples):
            bs_g1 = np.random.choice(group1, len(group1), replace=True)
            bs_g2 = np.random.choice(group2, len(group2), replace=True)
            m1, m2 = np.mean(bs_g1), np.mean(bs_g2)
            v1, v2 = np.var(bs_g1, ddof=1), np.var(bs_g2, ddof=1)
            p_std = np.sqrt((v1 + v2) / 2)
            d_val = (m1 - m2) / p_std if p_std > 0 else 0.0
            bs_d.append(d_val)
        return np.percentile(bs_d, [2.5, 97.5])

    ci = bootstrap_d(num_bootstrap)

    # Permutation Test
    def perm_p(perm_samples):
        obs_diff = abs(mean1 - mean2)
        perm_diffs = []
        all_data = np.concatenate([group1, group2])
        n1 = len(group1)
        for _ in range(perm_samples):
            np.random.shuffle(all_data)
            p_g1 = all_data[:n1]
            p_g2 = all_data[n1:]
            perm_diff = abs(np.mean(p_g1) - np.mean(p_g2))
            perm_diffs.append(perm_diff)
        p = np.mean([1 if diff >= obs_diff else 0 for diff in perm_diffs])
        return p

    p_perm = perm_p(num_perm)

    return {"cohens_d": cohens_d, "ci_low": ci[0], "ci_high": ci[1], "p_perm": p_perm}


def load_sequences(input_file):
    path = Path(input_file)
    seqs = []
    labels = []

    if path.suffix in [".fasta", ".fa"]:
        with open(input_file, "r") as f:
            current_seq = ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_seq:
                        seqs.append(current_seq)
                        labels.append("unknown")
                    current_seq = ""
                else:
                    current_seq += line
            if current_seq:
                seqs.append(current_seq)
                labels.append("unknown")
    else:
        with open(input_file, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if not row: continue
                if len(row) >= 1:
                    seqs.append(row[0])
                if len(row) >= 2:
                    labels.append(row[1])
                else:
                    labels.append("unknown")

    return seqs, labels


def main():
    parser = argparse.ArgumentParser(
        description="DNA Breathing Dynamics CZT Feature Extractor"
    )
    parser.add_argument(
        "--input",
        default=str(Path(__file__).parent / "sample.fasta"),
        help="Input file (FASTA/TSV/CSV; default: sample.fasta in gist dir)",
    )
    parser.add_argument("--output", default="features.csv", help="Output CSV")
    parser.add_argument(
        "--groups", help="Comma-separated group labels (e.g., groupA,groupB)"
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--at-lifetime", type=float, default=5.0)
    parser.add_argument("--gc-lifetime", type=float, default=25.0)
    parser.add_argument("--helical-period", type=float, default=10.5)
    parser.add_argument("--apply-helical", action="store_true", default=True)
    parser.add_argument(
        "--iupac-policy", default="average", choices=["average", "mask"]
    )
    parser.add_argument("--band-width", type=float, default=0.01)
    parser.add_argument("--apply-taper", action="store_true", default=False)
    parser.add_argument("--taper-k", type=float, default=0.3)
    parser.add_argument("--num-shuffles", type=int, default=10)
    parser.add_argument("--bootstrap-samples", type=int, default=100)
    parser.add_argument("--num-perm", type=int, default=20)

    args = parser.parse_args()

    np.random.seed(args.seed)
    random.seed(args.seed)

    seqs, labels = load_sequences(args.input)
    if not seqs:
        print("No sequences loaded.", file=sys.stderr)
        sys.exit(1)

    features_list = []
    for i, seq in enumerate(seqs):
        try:
            # 1. Original Sequence Analysis
            signal = encode_sequence(
                seq,
                args.at_lifetime,
                args.gc_lifetime,
                args.helical_period,
                args.apply_helical,
                args.iupac_policy,
            )
            freqs, spectrum = czt_analysis(
                signal,
                band_width=args.band_width,
                apply_taper=args.apply_taper,
                taper_k=args.taper_k,
            )
            features = extract_features(freqs, spectrum, 1 / 10.5, args.band_width)
            features["seq_id"] = i
            features["seq_len"] = len(seq)
            features["label"] = labels[i] if i < len(labels) else "unknown"

            # 2. Dinucleotide Shuffle Analysis
            shuffles = generate_dinuc_shuffles(seq, args.num_shuffles, args.seed + i)
            shuffle_features = []
            for shuf in shuffles:
                shuf_signal = encode_sequence(
                    shuf,
                    args.at_lifetime,
                    args.gc_lifetime,
                    args.helical_period,
                    args.apply_helical,
                    args.iupac_policy,
                )
                shuf_freqs, shuf_spectrum = czt_analysis(
                    shuf_signal,
                    band_width=args.band_width,
                    apply_taper=args.apply_taper,
                    taper_k=args.taper_k,
                )
                shuf_feat = extract_features(
                    shuf_freqs, shuf_spectrum, 1 / 10.5, args.band_width
                )
                shuffle_features.append(shuf_feat["peak_mag"])
            features["shuffle_mean_peak"] = np.mean(shuffle_features) if shuffle_features else 0.0

            # 3. Phase Scramble Null Analysis
            # Scramble the spectrum of the original signal
            scramble_peaks = []
            for k in range(args.num_shuffles):
                # Use same count as shuffles for consistency
                scrambled_spec = phase_scramble(spectrum, seed=args.seed + i + k*1000)
                # Re-extract features from scrambled spectrum
                scram_feat = extract_features(
                   freqs, scrambled_spec, 1 / 10.5, args.band_width
                )
                scramble_peaks.append(scram_feat["peak_mag"])
            features["phase_scramble_mean_peak"] = np.mean(scramble_peaks) if scramble_peaks else 0.0

            features_list.append(features)

        except Exception as e:
            print(f"Error processing seq {i}: {e}", file=sys.stderr)
            continue

    group_labels = None
    if args.groups:
        group_labels = args.groups.split(",")
    
    # If labels are provided in CSV, override CLI groups if CLI groups not enough
    # Actually, we prefer using the 'label' column from input if 'groups' arg not set
    if not group_labels and len(set(labels)) > 1:
        group_labels = labels

    stats_dict = compute_stats(
        features_list, group_labels, args.bootstrap_samples, args.num_perm, args.seed
    )

    fieldnames = [
        "seq_id",
        "seq_len",
        "label",
        "peak_mag",
        "peak_idx",
        "peak_freq",
        "peak_phase",
        "peak_to_skirt",
        "phase_coherence",
        "band_energy",
        "snr",
        "qc_flag",
        "shuffle_mean_peak",
        "phase_scramble_mean_peak",
    ] + list(stats_dict.keys())

    if args.output:
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for f in features_list:
                row = {k: v for k, v in f.items()}
                writer.writerow(row)
            stats_row = {"seq_id": "stats", "label": "stats", "seq_len": 0}
            stats_row.update(stats_dict)
            writer.writerow(stats_row)

    print("Features extracted and saved.")
    print(
        f"Cohen's d: {stats_dict['cohens_d']:.4f} [CI: {stats_dict['ci_low']:.4f}, {stats_dict['ci_high']:.4f}], p_perm: {stats_dict['p_perm']:.4f}"
    )


if __name__ == "__main__":
    main()
