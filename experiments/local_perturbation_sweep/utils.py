"""
Utilities for local perturbation sweep: mutant generation, diff computation.

No logic implemented—stubs for biophysical mutations and spectral analysis.
"""

from typing import Dict, List, Tuple
from pathlib import Path

# =============================================================================
# Mutant Generation Stub
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
    # TODO: Implement: 1-based pos; if base allows (e.g., inc: if A/T, A->G/T->C random);
    # else skip/raise; ensure exactly one change; use seed for ties
    pass
    return ""


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
    # TODO: Implement: Sequential single mutants; handle if position invalid
    pass
    return ""


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
    # TODO: Implement: Reuse src/core/params.encode_sequence; CZT from gist;
    # extract: resonance_mag=peak_mag, phase_coh=kappa, centroid, energy, snr;
    # return {f'delta_{k}': mut[k] - wt[k] for k in metrics}
    pass
    return {}


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
    # TODO: Implement: Parse FASTA; filter valid 20nt; compute gc; balanced subsample
    pass
    return []
