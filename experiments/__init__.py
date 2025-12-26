"""
Experiments module for DNA Breathing Dynamics framework.

This module contains experimental scripts and prototypes for
DNA spectral analysis using CZT-based feature extraction.

Experiments:
-----------
- prototype_alpha.py: Initial CZT feature extractor prototype
- phase_coherence_study.py: Phase-coherent vs random-phase CZT spectra comparison
- ddg_helical_validation.py: ΔΔG-Pairs + Helical Band Sweep validation framework

Study Results:
-------------
1. Phase Coherence Study (experiments/results/phase_coherence_study/FINDINGS.md):
   Phase coherence at the helical frequency (~10.5 bp/turn) does NOT provide
   meaningful predictive advantage over random-phase controls for GC-content
   classification (AUROC Δ = +0.0025, p = 0.596, Cohen's d = 0.11 [negligible]).

2. ΔΔG Helical Validation (artifacts/ddg_helical_validation/):
   DNA breathing dynamics encoding shows statistically detectable but practically
   negligible sensitivity to single-mutation thermodynamic perturbations (ΔΔG).
   Max Cohen's d = 0.061, far below pre-registered threshold of 0.5.
   Validation executed on 1,000 real Brunello CRISPR sequences (60,000 WT-mutant pairs).
"""
