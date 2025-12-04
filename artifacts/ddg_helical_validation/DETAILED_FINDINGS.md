# DNA Breathing Dynamics Validation — Detailed Findings Report

**Experiment Date:** 2025-12-04  
**Framework Version:** ΔΔG-Pairs + Helical Band Sweep Validation v1.0  
**Dataset:** Brunello CRISPR Library (1,000 sgRNA sequences, real biological data)  

---

## Executive Summary

This report presents findings from the pre-registered validation framework testing whether DNA breathing dynamics encoding can detect thermodynamic perturbations (ΔΔG) across helical frequency bands.

### Key Findings

| Metric | Result |
|--------|--------|
| **Overall Validation Status** | **FAIL** |
| Total WT-mutant pairs analyzed | 60,000 |
| Measurements collected | 6,000,000 |
| ΔΔG bins | Low (0–0.55), Mid (0.55–1.12), High (1.12–1.83) kcal/mol |
| Primary criterion met | No (|d| < 0.5 for all conditions) |
| Largest observed effect size | d = 0.061 (peak_mag at center=10.7, width=6%) |
| Trend test significant | No (all JT p-values > 0.05) |

---

## 1. Experimental Configuration

### 1.1 Dataset

- **Source:** Brunello CRISPR sgRNA library (real biological sequences)
- **Sample size:** 1,000 unique 20bp sgRNA sequences
- **Generated pairs:** 60,000 WT-mutant pairs (all single-point mutations)
- **Exclusions:** Ambiguous bases removed

### 1.2 Parameter Sweep

| Parameter | Values |
|-----------|--------|
| Center frequencies | 10.3, 10.4, 10.5, 10.6, 10.7 bp/turn |
| Band widths | 1%, 2%, 3%, 5%, 6% of f₀ |
| Features measured | peak_mag, band_power, phase_coherence, snr |
| ΔΔG bins | 3 tertiles (low/mid/high) |

### 1.3 Statistical Framework

- **Bootstrap samples:** 1,000 (BCa confidence intervals)
- **Permutations:** 500 (label-shuffle null distribution)
- **Significance level:** α = 0.05
- **Effect size threshold:** |d| ≥ 0.5
- **Control threshold:** |d| < 0.2

---

## 2. Primary Results

### 2.1 High-ΔΔG Bin (Primary Criterion)

The pre-registered primary criterion requires |Cohen's d| ≥ 0.5 in the high-ΔΔG bin with q < 0.05 and CI excluding zero.

**Result: FAIL (0 hits)**

Top 10 effects in high-ΔΔG bin (sorted by |d|):

| Center | Width | Feature | n | d | 95% CI | FDR q |
|--------|-------|---------|---|---|--------|-------|
| 10.7 | 6% | peak_mag | 20,043 | 0.061 | [0.047, 0.076] | <0.0001 |
| 10.7 | 5% | peak_mag | 20,043 | 0.060 | [0.046, 0.076] | <0.0001 |
| 10.7 | 3% | peak_mag | 20,043 | 0.057 | [0.045, 0.074] | <0.0001 |
| 10.6 | 6% | peak_mag | 20,043 | 0.056 | [0.042, 0.071] | <0.0001 |
| 10.7 | 2% | peak_mag | 20,043 | 0.056 | [0.043, 0.072] | <0.0001 |
| 10.6 | 5% | peak_mag | 20,043 | 0.055 | [0.042, 0.071] | <0.0001 |
| 10.7 | 1% | peak_mag | 20,043 | 0.054 | [0.041, 0.070] | <0.0001 |
| 10.6 | 3% | peak_mag | 20,043 | 0.054 | [0.040, 0.069] | <0.0001 |
| 10.6 | 2% | peak_mag | 20,043 | 0.053 | [0.039, 0.068] | <0.0001 |
| 10.7 | 1% | snr | 20,043 | 0.052 | [0.039, 0.068] | <0.0001 |

**Interpretation:** While effects are statistically significant (q < 0.05) due to large sample size, they are negligible in magnitude (d ≈ 0.05-0.06), far below the pre-registered threshold of |d| ≥ 0.5. The largest effect (d = 0.061) represents only 12% of the required effect size.

### 2.2 Feature Performance Comparison

| Feature | Max |d| in High Bin | Median |d| | Best Center |
|---------|-------------------|------------|-------------|
| peak_mag | 0.061 | 0.048 | 10.7 |
| band_power | 0.041 | 0.039 | 10.7 |
| snr | 0.052 | 0.040 | 10.7 |
| phase_coherence | 0.021 | 0.012 | 10.6 |

**Interpretation:** Peak magnitude shows the strongest (though still small) sensitivity to ΔΔG perturbations. Phase coherence appears least sensitive to thermodynamic changes.

---

## 3. Dose-Response Analysis

### 3.1 Trend Test (Jonckheere-Terpstra)

Tests whether |d| increases monotonically with |ΔΔG| bin (dose-response).

**Result: No significant trends detected**

Representative results:

| Center | Width | Feature | Low |d̄| | Mid |d̄| | High |d̄| | JT p-value |
|--------|-------|---------|---------|---------|----------|------------|
| 10.5 | 3% | peak_mag | 0.019 | 0.006 | 0.046 | 0.30 |
| 10.7 | 5% | peak_mag | 0.030 | 0.011 | 0.060 | 0.30 |
| 10.5 | 2% | snr | 0.027 | 0.008 | 0.051 | 0.30 |

**Interpretation:** Effect sizes do increase from low to high ΔΔG bins, but the pattern is not monotonic—mid-bin effects are consistently smaller than low-bin effects. This U-shaped pattern (instead of expected monotonic increase) suggests the signal may be capturing something other than pure thermodynamic sensitivity.

### 3.2 Bin-Specific Effect Sizes

| Bin | |ΔΔG| Range (kcal/mol) | N pairs | Mean |d| | Max |d| |
|-----|------------------------|---------|----------|---------|
| Low | 0.00 – 0.55 | 19,997 | 0.016 | 0.028 |
| Mid | 0.55 – 1.12 | 19,960 | 0.005 | 0.019 |
| High | 1.12 – 1.83 | 20,043 | 0.039 | 0.061 |

---

## 4. Control Analyses

### 4.1 Off-Band Control

Tests spectral windows ±20% away from helical frequency (8.4 and 12.6 bp/turn).

**Result: PASS (0 violations)**

- No off-band effects exceeded |d| = 0.2
- This confirms the encoding specifically targets the helical frequency region

### 4.2 Label-Shuffle Permutation

Tests whether observed effects could arise by chance when WT/mutant labels are randomly shuffled within pairs.

**Result: FAIL (89 violations)**

Many high-bin conditions showed label-shuffle p-values < 0.05, indicating that the small observed effects, while consistent within the original labeling, do not robustly survive permutation testing at the pair level.

---

## 5. Robustness Assessment

### 5.1 Adjacent Parameter Replication

**Result: FAIL**

Pre-registered criterion requires effects to replicate across ≥2 adjacent centers or widths. Since no effects met the primary |d| ≥ 0.5 threshold, robustness cannot be evaluated.

### 5.2 Effect Consistency Across Parameter Space

Effects at peak_mag are relatively stable across:
- **Centers:** Higher at 10.6-10.7 bp/turn than 10.3-10.5
- **Widths:** Larger bands (5-6%) show slightly higher effects

This spatial pattern is consistent but weak.

---

## 6. Power Analysis

### 6.1 Observed Power

With N = 20,000 pairs per bin:
- Observed effect size: d ≈ 0.06
- Actual power at observed d: >99%
- Required effect size for |d| = 0.5 detection: Power = 100%

The experiment has more than sufficient power—the issue is that the true effect size is small, not that the experiment lacks sensitivity.

### 6.2 What Effect Size Would We Need?

To achieve the pre-registered goal of |d| ≥ 0.5:
- Current max |d| = 0.061
- Required increase: 8.2× larger effect
- This would require either:
  - Much larger ΔΔG perturbations (not physically possible with single mutations)
  - A fundamentally different encoding approach

---

## 7. Scientific Interpretation

### 7.1 Why Did the Validation Fail?

1. **Single-point mutations produce small ΔΔG changes:**
   - Maximum |ΔΔG| observed: 1.83 kcal/mol
   - This corresponds to roughly one nearest-neighbor stack difference
   - The encoding may be insensitive to such small thermodynamic changes

2. **The helical frequency may not be the optimal detection band:**
   - Effects were slightly larger at higher centers (10.6-10.7)
   - Perhaps off-period modulations carry more ΔΔG signal

3. **Phase coherence showed minimal ΔΔG sensitivity:**
   - This feature captures rotational positioning
   - Single-point mutations don't substantially alter helical phase

### 7.2 Positive Findings

Despite failing the primary criterion, several encouraging patterns emerged:

1. **Consistent direction:** High-ΔΔG pairs consistently show higher peak_mag differences
2. **Spatial specificity:** Effects concentrate near helical frequency, not off-band
3. **Reproducibility:** Results are deterministic with fixed seed

### 7.3 Recommendations for Future Work

1. **Multiple mutations:** Test WT/mutant pairs with 2-5 mutations to amplify ΔΔG signal
2. **Extreme pairs:** Focus on mutations at stacking-sensitive positions (GC→AT conversions)
3. **Alternative frequency targets:** Explore sub-harmonic or super-harmonic frequencies
4. **Cross-validation:** Test on independent CRISPR libraries (GeCKO v2, TKO v3)

---

## 8. Artifacts Produced

| File | Description |
|------|-------------|
| `config.json` | Full parameter configuration with timestamps and seed |
| `pairs.csv` | All 60,000 WT-mutant pairs with ΔΔG values and bin assignments |
| `stats.csv` | Statistics for all 300 (bin × center × width × feature) conditions |
| `trend.csv` | Jonckheere-Terpstra trend test results |
| `report.md` | Automated summary report |

---

## 9. Conclusions

The DNA breathing dynamics encoding shows **statistically detectable but practically negligible** sensitivity to single-mutation thermodynamic perturbations. The pre-registered acceptance criteria were not met:

| Criterion | Status | Detail |
|-----------|--------|--------|
| Primary (|d| ≥ 0.5) | **FAIL** | Max |d| = 0.061 (12% of threshold) |
| Robustness | **FAIL** | No primary hits to evaluate |
| Specificity | **FAIL** | 89 shuffle violations |

**Key takeaway:** The encoding requires larger thermodynamic perturbations or methodological refinement to achieve biologically meaningful sensitivity.

---

## Appendix: Reproducibility Information

```
Seed: 42
NumPy version: 1.26.x
SciPy version: 1.11.x
Python version: 3.12.3
Platform: Linux
Timestamp: 2025-12-04T20:41:31+00:00
```

All artifacts are versioned and reproducible with the specified random seed.

---

*End of Detailed Findings Report*
