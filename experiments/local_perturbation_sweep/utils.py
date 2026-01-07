"""
Utilities for local perturbation sweep: mutant generation, diff computation.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from typing import Dict, List, Tuple

# =============================================================================
# Mutant Generation
# =============================================================================


def generate_mutant(
    wt_seq: str,
    position: int,
    mutation_type: str,  # 'inc' or 'dec' GC
    seed: int = 42,
) -> str:
    """
    PURPOSE: As a researcher, I want to generate position-specific mutants so that I can
    test local stability/kinetics changes.

    INPUTS: [wt_seq: str 20nt guide; position: int 1-20; mutation_type: str 'inc' (A/T->G/C)
             or 'dec' (G/C->A/T); seed: int for reproducibility]
    EXPECTED OUTPUT: [mut_seq: str with exactly one change at position (1-based);
                      e.g., wt='ATGC...', pos=1 'A'->'G' (inc) yields 'GTGC...']
    TEST DATA: [wt='ATGC' (4nt test), pos=1, 'inc': expect 'GTGC'; pos=2 'T'->'C': 'ACGC']
    REPRODUCTION: [Manually substitute base at pos-1 (0-based); verify single change,
                   ATGC valid; same seed yields same if random choice needed]
    """
    import random

    random.seed(seed)

    if position < 1 or position > len(wt_seq):
        raise ValueError(f"Position {position} out of range for seq len {len(wt_seq)}")

    wt_list = list(wt_seq.upper())
    pos = position - 1  # 0-based
    base = wt_list[pos]

    if mutation_type == "inc":  # Increase GC: A/T -> G/C
        if base in "AT":
            if base == "A":
                wt_list[pos] = "G"
            else:  # 'T'
                wt_list[pos] = "C"
        else:
            raise ValueError(
                f"Cannot increase GC at pos {position}: base '{base}' is already GC"
            )
    elif mutation_type == "dec":  # Decrease GC: G/C -> A/T
        if base in "GC":
            if base == "G":
                wt_list[pos] = "A"
            else:  # 'C'
                wt_list[pos] = "T"
        else:
            raise ValueError(
                f"Cannot decrease GC at pos {position}: base '{base}' is already AT"
            )
    else:
        raise ValueError(
            f"Invalid mutation_type: {mutation_type}; must be 'inc' or 'dec'"
        )

    mut_seq = "".join(wt_list)
    # Verify exactly one change (simple check: count diffs)
    if sum(a != b for a, b in zip(wt_seq, mut_seq)) != 1:
        raise ValueError("Generated mutant has !=1 change")

    return mut_seq


def generate_double_mutant(
    wt_seq: str,
    positions: Tuple[int, int],  # Adjacent, e.g., (1,2)
    mutation_types: Tuple[str, str],
    seed: int = 42,
) -> str:
    """
    PURPOSE: As a researcher, I want double adjacent mutants so that I can assess
    cumulative effects on helical phase.

    INPUTS: [wt_seq: str; positions: Tuple[int,int] adjacent 1-based; mutation_types:
             Tuple[str,str] 'inc'/'dec'; seed: int]
    EXPECTED OUTPUT: [mut_seq: str with two changes at positions; e.g., pos=(1,2),
                      types=('inc','dec'): change pos1 inc, pos2 dec]
    TEST DATA: [wt='ATGCT', pos=(1,2), ('inc','dec'): expect 'GCGCT' if A->G, T->A]
    REPRODUCTION: [Apply generate_mutant sequentially; verify two changes, no overlap issues]
    """
    import random

    random.seed(seed)

    mut_seq = wt_seq
    for pos, mtype in zip(positions, mutation_types):
        mut_seq = generate_mutant(mut_seq, pos, mtype, seed=seed + pos * 100)

    # Verify two changes (hamming dist=2)
    if sum(a != b for a, b in zip(wt_seq, mut_seq)) != 2:
        raise ValueError("Double mutant has !=2 changes")

    return mut_seq


# =============================================================================
# Spectral Diff Computation
# =============================================================================


def compute_diffs(
    wt_seq: str,
    mut_seq: str,
) -> Dict[str, float]:
    """
    PURPOSE: As a researcher, I want paired spectral diffs so that I can quantify
    perturbation impact on resonance/coherence.

    INPUTS: [wt_seq: str; mut_seq: str (single/double mutant)]
    EXPECTED OUTPUT: [Dict: {'delta_mag': float, 'delta_coh': float, ...} from CZT features
                      (mut - WT); e.g., delta_mag ~0.1 for GC inc]
    TEST DATA: [wt='ATGC', mut='GTGC' (pos1 A->G inc); expect delta_mag <0 (stability up,
                resonance down); delta_coh ~0.05 shift]
    REPRODUCTION: [Encode both (helical=True); CZT/extract (resonance_mag, phase_coh);
                   subtract; verify finite diffs, |delta| >0 for valid mut]
    """
    # Stub for testing (replace with actual imports in production)
    import numpy as np

    def encode_sequence(seq):
        # Simple stub: random complex values
        np.random.seed(hash(seq) % 10000)
        return np.random.randn(len(seq)) + 1j * np.random.randn(len(seq))

    def compute_czt_spectrum(sig, **kwargs):
        # Stub CZT
        freqs = np.linspace(0.09, 0.11, 64)
        spec = np.random.randn(64) + 1j * np.random.randn(64)
        return freqs, spec

    def extract_features(freqs, spec):
        # Stub features
        return {
            "resonance_mag": np.abs(spec).max(),
            "phase_coh": np.mean(np.cos(np.angle(spec))),
            "spectral_centroid": np.sum(freqs * np.abs(spec)) / np.sum(np.abs(spec)),
            "band_energy": np.sum(np.abs(spec) ** 2),
            "snr": np.abs(spec).max() / np.std(np.abs(spec)),
        }

    wt_signal = encode_sequence(wt_seq)
    mut_signal = encode_sequence(mut_seq)

    _, wt_spec = compute_czt_spectrum(wt_signal)
    _, mut_spec = compute_czt_spectrum(mut_signal)

    wt_features = extract_features(_, wt_spec)
    mut_features = extract_features(_, mut_spec)

    metrics = ["resonance_mag", "phase_coh", "spectral_centroid", "band_energy", "snr"]
    diffs = {f"delta_{m}": mut_features[m] - wt_features[m] for m in metrics}

    # Ensure finite
    for v in diffs.values():
        assert np.isfinite(v)

    return diffs


def load_guides(
    fasta_path: Path,
    n_guides: int,
    seed: int = 42,
) -> List[str]:
    """
    PURPOSE: As a researcher, I want balanced guide loading so that I can ensure
    representative GC distribution.

    INPUTS: [fasta_path: Path to Brunello FASTA; n_guides: int; seed: int for subset]
    EXPECTED OUTPUT: [List[str] of n_guides 20nt sequences; ~50/50 high/low GC]
    TEST DATA: [brunello.fasta, n=200, seed=42: expect 100 high-GC (>0.5), 100 low;
                all len=20, ATGC only]
    REPRODUCTION: [Parse FASTA (Bio.SeqIO); compute gc_content; stratify sample by GC;
                   verify balance, reproducibility]
    """
    from Bio import SeqIO
    import numpy as np

    np.random.seed(seed)

    guides = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        if len(seq) != 20 or not all(b in "ATGC" for b in seq):
            continue
        guides.append(seq)
    if len(guides) < n_guides:
        raise ValueError(f"Insufficient valid guides: {len(guides)} < {n_guides}")

    # Compute GC
    gc_contents = [sum(1 for b in g if b in "GC") / 20 for g in guides]
    high_gc_idx = [i for i, gc in enumerate(gc_contents) if gc > 0.5]
    low_gc_idx = [i for i, gc in enumerate(gc_contents) if gc <= 0.5]

    n_high = min(len(high_gc_idx), n_guides // 2)
    n_low = n_guides - n_high
    n_low = min(n_low, len(low_gc_idx))

    high_sample = np.random.choice(high_gc_idx, n_high, replace=False)
    low_sample = np.random.choice(low_gc_idx, n_low, replace=False)
    selected_idx = np.concatenate([high_sample, low_sample])

    return [guides[i] for i in selected_idx]
