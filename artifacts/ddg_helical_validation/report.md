# ΔΔG-Pairs + Helical Band Sweep Validation Report

**Status: FAIL ✗**

Generated: 2025-12-05T06:23:46.031265+00:00

---

## Summary

This report evaluates the DNA breathing dynamics encoding's sensitivity to
thermodynamic perturbations (ΔΔG) across helical frequency bands.

### Configuration

- **Seed**: 42
- **Bootstrap samples**: 1000
- **Permutations**: 500
- **Center sweep**: [10.3, 10.4, 10.5, 10.6, 10.7]
- **Width sweep**: [0.01, 0.02, 0.03, 0.05, 0.06]
- **Temperature**: 310.15 K

### Data

- **Total pairs**: 60000
- **Bin edges (|ΔΔG| kcal/mol)**: [0.0, 0.55, 1.12, 1.83]

---

## Acceptance Criteria Evaluation

| Criterion | Result |
|-----------|--------|
| Primary (|d| ≥ 0.5, q < 0.05, CI excl. 0 in high bin) | FAIL ✗ (0 hits) |
| Robustness (replicates across ≥2 adjacent params) | FAIL ✗ |
| Specificity (off-band/shuffle non-significant) | FAIL ✗ |

### Control Violations

- Off-band violations: 0
- Shuffle violations: 89

---

## Top Results (High Bin, Sorted by |d|)

| Center | Width | Feature | n | d | 95% CI | q |
|--------|-------|---------|---|---|--------|---|
| 10.7 | 0.06 | peak_mag | 20043 | 0.061 | [0.047, 0.076] | 0.0000 |
| 10.7 | 0.05 | peak_mag | 20043 | 0.060 | [0.046, 0.076] | 0.0000 |
| 10.7 | 0.03 | peak_mag | 20043 | 0.057 | [0.045, 0.074] | 0.0000 |
| 10.6 | 0.06 | peak_mag | 20043 | 0.056 | [0.042, 0.071] | 0.0000 |
| 10.7 | 0.02 | peak_mag | 20043 | 0.056 | [0.043, 0.072] | 0.0000 |
| 10.6 | 0.05 | peak_mag | 20043 | 0.055 | [0.042, 0.071] | 0.0000 |
| 10.7 | 0.01 | peak_mag | 20043 | 0.054 | [0.041, 0.070] | 0.0000 |
| 10.6 | 0.03 | peak_mag | 20043 | 0.054 | [0.040, 0.069] | 0.0000 |
| 10.6 | 0.02 | peak_mag | 20043 | 0.053 | [0.039, 0.068] | 0.0000 |
| 10.7 | 0.01 | snr | 20043 | 0.052 | [0.039, 0.068] | 0.0000 |

---

## Trend Test Results (Jonckheere-Terpstra)

| Center | Width | Feature | Low |d̄| | Mid |d̄| | High |d̄| | JT p-value |
|--------|-------|---------|---------|---------|----------|------------|
| 10.3 | 0.01 | peak_mag | 0.018 | 0.004 | 0.040 | 0.3008 |
| 10.3 | 0.01 | band_power | 0.019 | 0.003 | 0.041 | 0.3008 |
| 10.3 | 0.01 | phase_coherence | 0.008 | 0.000 | 0.010 | 0.3008 |
| 10.3 | 0.01 | snr | 0.020 | 0.003 | 0.039 | 0.3008 |
| 10.3 | 0.02 | peak_mag | 0.016 | 0.005 | 0.042 | 0.3008 |
| 10.3 | 0.02 | band_power | 0.019 | 0.003 | 0.041 | 0.3008 |
| 10.3 | 0.02 | phase_coherence | 0.001 | 0.007 | 0.013 | 0.0586 |
| 10.3 | 0.02 | snr | 0.020 | 0.003 | 0.039 | 0.3008 |
| 10.3 | 0.03 | peak_mag | 0.013 | 0.006 | 0.043 | 0.3008 |
| 10.3 | 0.03 | band_power | 0.019 | 0.003 | 0.041 | 0.3008 |
| 10.3 | 0.03 | phase_coherence | 0.002 | 0.008 | 0.012 | 0.0586 |
| 10.3 | 0.03 | snr | 0.020 | 0.003 | 0.039 | 0.3008 |
| 10.3 | 0.05 | peak_mag | 0.006 | 0.006 | 0.045 | 0.0586 |
| 10.3 | 0.05 | band_power | 0.018 | 0.003 | 0.041 | 0.3008 |
| 10.3 | 0.05 | phase_coherence | 0.008 | 0.006 | 0.012 | 0.3008 |

---

## Pre-registration Notes

This validation follows pre-registered acceptance criteria:

1. **Primary**: At least one (center, width, feature) combination in the
   high-ΔΔG bin must achieve |d| ≥ 0.5, FDR-corrected q < 0.05, and
   95% BCa confidence interval excluding zero.

2. **Robustness**: The effect must replicate across at least 2 adjacent
   center frequencies or band widths.

3. **Specificity**: Off-band windows (±20% from f₀) and label-shuffle
   permutation controls must remain non-significant (|d| < 0.2; q ≥ 0.1).

---

*End of report*
