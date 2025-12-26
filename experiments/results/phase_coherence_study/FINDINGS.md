# Phase-Coherent vs Random-Phase CZT Spectra: Study Findings

**Study Date:** December 5, 2024  
**Dataset:** Brunello CRISPR library subset (n=1000 gRNAs, balanced 500 high-GC / 500 low-GC)  
**Primary Outcome:** Binary GC-content classification as efficiency proxy

---

## Executive Summary

**Conclusion: Phase coherence at the helical frequency (~10.5 bp/turn) does NOT provide meaningful predictive advantage over random-phase controls for this classification task.**

The study found:
- **AUROC Δ = +0.0025** (negligible effect, Cohen's d = 0.11)
- **p-value = 0.596** (not statistically significant)
- **Achieved power = 8.1%** (severely underpowered)
- **Effect interpretation: Negligible**

---

## Study Design

| Parameter | Value |
|-----------|-------|
| Samples | 1,000 gRNA sequences |
| Labels | 500 high-GC, 500 low-GC (balanced) |
| CV Strategy | 5-fold × 5 repeats stratified |
| Model | Logistic regression (no tuning) |
| Seeds | global=137, phase=271828, cv=161803 |
| Band Center | 1/10.5 = 0.0952 cycles/bp |
| Band Width | ±0.01 cycles/bp |

### Features Extracted (7 total)
1. Peak magnitude
2. Peak frequency index  
3. Phase at peak
4. Band-integrated power
5. Spectral centroid
6. Von Mises κ (circular phase concentration)
7. Rayleigh z-statistic

---

## Primary Results

### AUROC Analysis

| Metric | Value |
|--------|-------|
| Mean Δ (A−B) | +0.0025 |
| Std Δ | 0.0230 |
| Cohen's d | 0.107 |
| 95% CI | [−0.312, +0.488] |
| Paired t-test p | 0.596 |
| Wilcoxon p | 0.692 |
| **Effect Size** | **Negligible** |

### Secondary Metrics

| Metric | Mean Δ | Cohen's d | p-value | Interpretation |
|--------|--------|-----------|---------|----------------|
| AUPRC | +0.0097 | 0.324 | 0.119 | Small (not significant) |
| Brier Score | +0.0008 | 0.215 | 0.293 | Small (not significant) |

---

## Power Analysis

| Parameter | Value |
|-----------|-------|
| Observed effect size (d) | 0.107 |
| Achieved power | **8.1%** |
| Desired power | 80% |
| N folds needed for 80% power | **680** |

**Interpretation:** The study is severely underpowered to detect the observed effect. However, the effect size itself is negligible (d < 0.2), meaning even with adequate power, the practical significance would be minimal.

---

## Robustness Checks

### Ablation Study: Phase-Only vs Magnitude-Only

| Condition | Mean AUROC | Std |
|-----------|------------|-----|
| Phase-only features | 0.542 | 0.033 |
| Magnitude-only features | **0.601** | 0.028 |

**Finding:** Magnitude features alone outperform phase features alone by ~6 percentage points AUROC. Phase information adds little predictive value for this task.

### Label Permutation Test

| Parameter | Value |
|-----------|-------|
| Observed Δ AUROC | +0.0132 |
| Permutation mean | −0.0012 |
| Permutation std | 0.0234 |
| Permutation p-value | **0.264** |

**Finding:** The observed advantage is not distinguishable from random label shuffling (p > 0.05), confirming the null result.

---

## Interpretation

### Why No Phase Coherence Advantage?

1. **Task mismatch:** GC-content classification may not require helical phase information. GC-content is a sequence composition property captured by magnitude/stability features.

2. **Feature redundancy:** The 7 spectral features may already capture relevant information in their magnitudes, making phase redundant.

3. **Small effect:** Any true phase coherence effect is likely <0.01 AUROC, requiring thousands of folds to detect with adequate power.

### What This Means for DNA Breathing Dynamics

- Phase coherence at the helical frequency does not appear to drive GC-classification performance
- Magnitude and stability features are the primary predictors
- For tasks requiring helical positioning (e.g., CRISPR on-target prediction with positional effects), phase may still be relevant—this study does not rule that out

---

## Recommendations

1. **Do not prioritize phase features** for GC-content or simple composition-based classification

2. **Focus on magnitude/stability encoding** which shows stronger predictive signal (0.60 AUROC)

3. **Test phase coherence on different outcomes:**
   - Actual CRISPR efficiency scores (continuous)
   - On-target vs off-target classification
   - Mismatch tolerance prediction

4. **Consider different frequency bands** or multi-band analysis

---

## Reproducibility

All seeds and parameters are fixed. To replicate:

```bash
python experiments/phase_coherence_study.py \
    --input gists/breathing/czt_feature_extractor/data/processed/brunello_gc1000.fasta \
    --output-dir experiments/results/phase_coherence_study \
    --n-folds 5 --n-repeats 5 \
    --n-bootstrap 1000 --n-permutations 500
```

---

## Files Generated

| File | Description |
|------|-------------|
| `fold_metrics.csv` | Per-fold AUROC, AUPRC, Brier scores for conditions A and B |
| `summary.json` | Complete statistical summary |
| `spaghetti_auroc.png` | Paired fold comparison plot |
| `delta_boxplot.png` | Distribution of Δ AUROC across folds |
| `FINDINGS.md` | This report |

---

## Conclusion

**H₁ (phase coherence advantage) is NOT supported.** The phase-coherent condition shows no statistically or practically significant advantage over random-phase controls for predicting GC-content labels. Effect size is negligible (d=0.11), and the ablation study confirms magnitude features carry the predictive signal.

This null result is scientifically informative: it suggests that for sequence composition tasks, the CZT magnitude spectrum—not phase—captures the relevant biophysical information encoded in DNA breathing dynamics.
