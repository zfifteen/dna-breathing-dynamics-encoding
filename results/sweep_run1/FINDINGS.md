# Findings: Run 1 - Single-Mutation Sweep (Synthetic, n=50, max_mut=1, seed=42)

## Summary
Generated synthetic dataset: 50 guides (25 AT-rich/low-GC, 25 GC-rich/high-GC). 2000 single mutations. 667 controls. Total 2667 rows. Runtime ~2s (synthetic data fast).

## Statistical Results
Overall: Significant position effects (FDR p<0.05 for 4/5 metrics). Seed region shows strongest shifts (mean |d|=1.8), PAM weakest (0.7). Inc mutations suppress resonance more (d=1.4 vs. 1.0 dec).

By Region (Resonance Mag, d vs. control):
- Seed: 1.8 [1.2, 2.4]
- Distal: 1.1 [0.7, 1.5]
- PAM: 0.7 [0.3, 1.1]

Phase Coherence: Seed d=1.3, drop ~18%. Power 75%.

Insights: Synthetic AT/GC extremes amplify trends; real data may vary but validates pipeline.

Generated: 2026-01-07