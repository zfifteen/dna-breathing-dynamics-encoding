"""
Utilities for local perturbation sweep: mutant generation, diff computation.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from typing import Dict, List, Tuple

# =============================================================================
# Mutant Generation Stub
# =============================================================================


def generate_mutant(\n    wt_seq: str,\n    position: int,\n    mutation_type: str,  # 'inc' or 'dec' GC\n    seed: int = 42,\n) -> str:\n    \"\"\"\n    PURPOSE: As a researcher, I want to generate position-specific mutants so that I can\n    test local stability/kinetics changes.\n\n    INPUTS: [wt_seq: str 20nt guide; position: int 1-20; mutation_type: str 'inc' (A/T->G/C)\n             or 'dec' (G/C->A/T); seed: int for reproducibility]\n    EXPECTED OUTPUT: [mut_seq: str with exactly one change at position (1-based);\n                      e.g., wt='ATGC...', pos=1 'A'->'G' (inc) yields 'GTGC...']\n    TEST DATA: [wt='ATGC' (4nt test), pos=1, 'inc': expect 'GTGC'; pos=2 'T'->'C': 'ACGC']\n    REPRODUCTION: [Manually substitute base at pos-1 (0-based); verify single change,\n                   ATGC valid; same seed yields same if random choice needed]\n    \"\"\"\n    import random\n    random.seed(seed)\n\n    if position < 1 or position > len(wt_seq):\n        raise ValueError(f\"Position {position} out of range for seq len {len(wt_seq)}\")\n\n    wt_list = list(wt_seq.upper())\n    pos = position - 1  # 0-based\n    base = wt_list[pos]\n\n    if mutation_type == 'inc':  # Increase GC: A/T -> G/C\n        if base in 'AT':\n            if base == 'A':\n                wt_list[pos] = 'G'\n            else:  # 'T'\n                wt_list[pos] = 'C'\n        else:\n            raise ValueError(f\"Cannot increase GC at pos {position}: base '{base}' is already GC\")\n    elif mutation_type == 'dec':  # Decrease GC: G/C -> A/T\n        if base in 'GC':\n            if base == 'G':\n                wt_list[pos] = 'A'\n            else:  # 'C'\n                wt_list[pos] = 'T'\n        else:\n            raise ValueError(f\"Cannot decrease GC at pos {position}: base '{base}' is already AT\")\n    else:\n        raise ValueError(f\"Invalid mutation_type: {mutation_type}; must be 'inc' or 'dec'\")\n\n    mut_seq = ''.join(wt_list)\n    # Verify exactly one change (simple check: count diffs)\n    if sum(a != b for a, b in zip(wt_seq, mut_seq)) != 1:\n        raise ValueError(\"Generated mutant has !=1 change\")\n\n    return mut_seq


def generate_double_mutant(\n    wt_seq: str,\n    positions: Tuple[int, int],  # Adjacent, e.g., (1,2)\n    mutation_types: Tuple[str, str],\n    seed: int = 42,\n) -> str:\n    \"\"\"\n    PURPOSE: As a researcher, I want double adjacent mutants so that I can assess\n    cumulative effects on helical phase.\n\n    INPUTS: [wt_seq: str; positions: Tuple[int,int] adjacent 1-based; mutation_types:\n             Tuple[str,str] 'inc'/'dec'; seed: int]\n    EXPECTED OUTPUT: [mut_seq: str with two changes at positions; e.g., pos=(1,2),\n                      types=('inc','dec'): change pos1 inc, pos2 dec]\n    TEST DATA: [wt='ATGCT', pos=(1,2), ('inc','dec'): expect 'GCGCT' if A->G, T->A]\n    REPRODUCTION: [Apply generate_mutant sequentially; verify two changes, no overlap issues]\n    \"\"\"\n    import random\n    random.seed(seed)\n\n    mut_seq = wt_seq\n    for pos, mtype in zip(positions, mutation_types):\n        mut_seq = generate_mutant(mut_seq, pos, mtype, seed=seed + pos)  # Offset seed per mut\n\n    # Verify two changes (hamming dist=2)\n    if sum(a != b for a, b in zip(wt_seq, mut_seq)) != 2:\n        raise ValueError(\"Double mutant has !=2 changes\")\n\n    return mut_seq


# =============================================================================
# Spectral Diff Computation Stub
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
    # In full impl: import from src/core/params.py and gist
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
            'resonance_mag': np.abs(spec).max(),
            'phase_coh': np.mean(np.cos(np.angle(spec))),
            'spectral_centroid': np.sum(freqs * np.abs(spec)) / np.sum(np.abs(spec)),
            'band_energy': np.sum(np.abs(spec)**2),
            'snr': np.abs(spec).max() / np.std(np.abs(spec)),
        }

    wt_signal = encode_sequence(wt_seq)
    mut_signal = encode_sequence(mut_seq)

    _, wt_spec = compute_czt_spectrum(wt_signal)
    _, mut_spec = compute_czt_spectrum(mut_signal)

    wt_features = extract_features(_, wt_spec)
    mut_features = extract_features(_, mut_spec)

    metrics = ['resonance_mag', 'phase_coh', 'spectral_centroid', 'band_energy', 'snr']
    diffs = {f'delta_{m}': mut_features[m] - wt_features[m] for m in metrics}

    # Ensure finite
    for v in diffs.values():
        assert np.isfinite(v)

    return diffs


def load_guides(\n    fasta_path: Path,\n    n_guides: int,\n    seed: int = 42,\n) -> List[str]:\n    \"\"\"\n    PURPOSE: As a researcher, I want balanced guide loading so that I can ensure\n    representative GC distribution.\n\n    INPUTS: [fasta_path: Path to Brunello FASTA; n_guides: int; seed: int for subset]\n    EXPECTED OUTPUT: [List[str] of n_guides 20nt sequences; ~50/50 high/low GC]\n    TEST DATA: [brunello.fasta, n=200, seed=42: expect 100 high-GC (>0.5), 100 low;\n                all len=20, ATGC only]\n    REPRODUCTION: [Parse FASTA (Bio.SeqIO); compute gc_content; stratify sample by GC;\n                   verify balance, reproducibility]\n    \"\"\"\n    from Bio import SeqIO  # Assumes biopython available\n    import numpy as np\n    np.random.seed(seed)\n\n    guides = []\n    for rec in SeqIO.parse(fasta_path, \"fasta\"):\n        seq = str(rec.seq).upper()\n        if len(seq) != 20 or not all(b in 'ATGC' for b in seq):\n            continue\n        guides.append(seq)\n    if len(guides) < n_guides:\n        raise ValueError(f\"Insufficient valid guides: {len(guides)} < {n_guides}\")\n\n    # Compute GC\n    gc_contents = [sum(1 for b in g if b in 'GC') / 20 for g in guides]\n    high_gc_idx = [i for i, gc in enumerate(gc_contents) if gc > 0.5]\n    low_gc_idx = [i for i, gc in enumerate(gc_contents) if gc <= 0.5]\n\n    n_high = min(len(high_gc_idx), n_guides // 2)\n    n_low = n_guides - n_high\n    n_low = min(n_low, len(low_gc_idx))\n\n    high_sample = np.random.choice(high_gc_idx, n_high, replace=False)\n    low_sample = np.random.choice(low_gc_idx, n_low, replace=False)\n    selected_idx = np.concatenate([high_sample, low_sample])\n\n    return [guides[i] for i in selected_idx]
