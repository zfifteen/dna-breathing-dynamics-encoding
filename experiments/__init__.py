"""
Experiments module for DNA Breathing Dynamics framework.

This module contains experimental scripts and prototypes for
DNA spectral analysis using CZT-based feature extraction.

Experiments:
-----------
- prototype_alpha.py: Initial CZT feature extractor prototype
- phase_coherence_study.py: Phase-coherent vs random-phase CZT spectra comparison

Study Results:
-------------
See experiments/results/phase_coherence_study/FINDINGS.md for the phase
coherence study conclusions. Key finding: Phase coherence at the helical
frequency (~10.5 bp/turn) does NOT provide meaningful predictive advantage
over random-phase controls for GC-content classification (AUROC Δ = +0.0025,
p = 0.596, Cohen's d = 0.11 [negligible]).
"""
