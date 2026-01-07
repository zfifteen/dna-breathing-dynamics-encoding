# Local Perturbation Sweep Experiment: Resonance Shifts in CRISPR Guides

## Overview
This experiment systematically evaluates the impact of single and double nucleotide perturbations on spectral features (resonance magnitude, phase coherence, spectral centroid, band energy) in CRISPR guide RNA sequences. It builds on the biophysical encoding framework (breathing kinetics + thermodynamic stability) and CZT analysis to quantify helical sensitivity by position and mutation type. 

**Rationale for Lab Scientists**: DNA breathing dynamics influence Cas9 binding efficiency and off-target effects. Mutations altering GC content (e.g., A/T ↔ G/C) change local stability (ΔG via SantaLucia parameters) and kinetics (lifetimes: AT=5ms, GC=25ms), potentially disrupting helical phase alignment (10.5 bp/turn). This sweep maps position-specific shifts (seed 1-8, distal 9-17, PAM-proximal 18-20), revealing if perturbations near functional regions (e.g., PAM) cause larger coherence drops (>20% expected in seed), informing guide design for reduced off-targets.

**Rationale for Software Engineers**: Extends `issue47_local_perturbation.py` with parametric sweeps (positions, types, multiples) for ~15k sequences from 200 Brunello guides. Ensures reproducibility (seed=42, fixed params.py), vectorized CZT (scipy.signal), and standardized outputs (CSV/JSON) for downstream ML/stats. No new deps; integrates with existing `dna_breathing_gist.py` for encoding/CZT.

**Expected Outcomes**: Heatmaps of Cohen's d (|d| ≥1.5 threshold) by position/type; null comparisons (random-pos mutants, phase-scramble) via permutation p<0.05; power analysis (target 80% for d=0.5). Hypothesis: Seed mutations yield d>2 (phase disruption); PAM-proximal <1 (tolerance).

## Experimental Design
- **Dataset**: 200 balanced Brunello CRISPR guides (20nt, ~50% high/low GC from `data/processed/brunello.fasta`).
- **Perturbations** (Sweep Parameters):
  - **Positions**: 1-20 (regions: seed=1-8, distal=9-17, PAM=18-20).
  - **Types**: Increase GC (A/T→G/C), Decrease (G/C→A/T); 1 per position.
  - **Multiples**: Doubles (adjacent pairs, 5/region: e.g., pos 1-2, 3-4).
  - Total: 200 guides × (60 singles + 15 doubles) = ~15k mutants.
- **Null Models** (Controls):
  1. Random-position: Same mutation, shuffled pos (preserves type/count).
  2. Non-GC: Transversions (e.g., A→C, preserves GC).
  3. Phase-scramble: WT spectrum magnitudes with uniform random phases [-π,π].
  4. Bootstrap: Resample WT diffs (5000 iters).
- **Processing Pipeline**:
  1. Load WT guide.
  2. Generate mutant (exact 1/2 changes; validate ATGC).
  3. Encode both (helical=True; params.py: SantaLucia ΔG, lifetimes).
  4. CZT (band ±0.01 at 0.095 cycles/bp, m=256 points).
  5. Extract: resonance_mag (peak_mag), phase_coh (von Mises κ), spectral_centroid, band_energy, snr.
  6. Compute diffs (mut - WT).
- **Analysis**:
  - **Per Position/Region/Type**: Paired t/Wilcoxon (p<0.05); Cohen's d + Hedges' g; bootstrap CIs (1000 iters); FDR across metrics (Benjamini-Hochberg, α=0.05).
  - **Null Tests**: Permutation p (1000 iters vs. random-pos null; expect p>0.05 if no position effect).
  - **Power**: Post-hoc (achieved for observed d); compute N needed (80% power at d=0.5).
  - **Thresholds**: Significant if |d|≥1.5, CI excludes 0, p<0.01.
- **Validation**: Unit tests (e.g., diffs=0 for no-mut; finite features); smoke on 10 guides; compare to issue47 (expect similar single-mut d~1-2).

## Implementation Plan
**Software Engineers**: 
- **New Script**: `experiments/local_perturbation_sweep/sweep_datagen.py` (augment issue47: loop over positions/types; vectorize mutations with BioPython if needed).
- **Utilities**: Reuse `src/core/params.py` (validate seq, params); `dna_breathing_gist.py` (encode + CZT/extract).
- **Reproducibility**: Fixed seed=42; log all params to metadata.json (grids: positions=[1-20], types=['inc','dec'], regions=['seed','distal','pam']).
- **Execution**: `python sweep_datagen.py --guides brunello_subset.fasta --output-dir results/sweep/ --n_guides 200 --max_mut 2`. ~30min on standard machine (parallelize with multiprocessing for 1000+ iters).
- **Testing**: Pytest in `tests/integration/test_perturbation_sweep.py`: Verify shapes (~15k rows), finite diffs, CI coverage (95% contains true d).
- **Outputs**:
  - `sweep_data.csv`: seq_id, wt_seq, mut_pos, mut_type, num_mut, region, resonance_mag, phase_coh, spectral_centroid, band_energy, snr, gc_content, seed, control_flag, delta_mag, delta_coh.
  - `metadata.json`: Experiment name, param grids, metric defs (e.g., "resonance_mag": "Peak |CZT| at helical freq, unit: arbitrary"), units, generated_at.
  - Plots: `d_heatmap.png` (d by pos/type); `power_analysis.png`.
- **Edge Cases**: Handle invalid muts (e.g., no A/T to change); short guides (<20nt); all-GC/AT baselines.

**Lab Scientists**:
- **Biophysical Interpretation**: Focus on Δcoh (phase disruption proxy for breathing accessibility); expect larger |Δmag| for GC-increase (stability up, resonance down). Validate vs. literature (e.g., AT mismatches tolerated in breathing-prone regions).
- **Running**: Use provided script; inspect CSVs in Jupyter (e.g., `pd.read_csv('sweep_data.csv').groupby('region')['delta_coh'].mean()`). Re-run with `--max_mut=1` for singles-only.
- **Extensions**: Add temperature sweep (310K±10°C on ΔG); ion effects (Na+/Mg2+ on stability). Compare to wet-lab Cas9 cleavage data for validation.

## Risks & Mitigations
- **Compute**: If slow, subsample to 100 guides; parallelize CZT.
- **Null Overlap**: If p<0.05 vs. nulls, increase permutations (to 5000).
- **Power**: If achieved<50%, report N_needed; suggest larger dataset.

## Next Steps
1. Implement `sweep_datagen.py` (Category B: extend issue47).
2. Run pilot (10 guides) for validation.
3. Full run; analyze in `analysis.ipynb` (TBD).
4. PR for integration if d>1.5 in seed region.

Last Updated: 2026-01-07 (ExperimentOrchestrator v1.0)