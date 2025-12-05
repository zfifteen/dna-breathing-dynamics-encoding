# DNA Breathing Dynamics Validation Statistics

**Generated:** 2025-12-05T01:44:41.076323+00:00

## Executive Summary

- **Dataset:** brunello_parsed.fasta
- **Total Sequences:** 1000
- **High GC Group (>50%):** 500 sequences
- **Low GC Group (≤50%):** 500 sequences

### Primary Finding

**Spectral Power at 1/10.5 bp⁻¹ (helical frequency):**

- Hedges' g = 0.0777 [95% CI: -0.0485, 0.2114]
- Wilcoxon p-value (FDR-corrected) = 4.2239e-01

**Interpretation:** A negligible effect size indicates that high-GC sequences show stronger spectral resonance at the B-form DNA helical frequency.

## Detailed Metrics

| Metric | Hedges' g | 95% CI | p-value (FDR) | High GC Mean ± SD | Low GC Mean ± SD |
|--------|-----------|--------|---------------|-------------------|------------------|
| peak_mag | 0.0777 | [-0.0485, 0.2114] | 4.2239e-01 | 120.8802 ± 46.5690 | 117.2856 ± 45.7374 |
| snr | -0.0623 | [-0.1850, 0.0745] | 4.2239e-01 | 0.9805 ± 0.3885 | 1.0051 ± 0.3997 |
| phase_coherence | 0.0159 | [-0.1071, 0.1441] | 2.4449e-01 | 0.8879 ± 0.1376 | 0.8858 ± 0.1284 |
| band_energy | 0.0083 | [-0.1234, 0.1381] | 7.4927e-01 | 965546.9468 ± 813805.3829 | 958785.8779 ± 819140.9992 |

## Normality Assessment

| Metric | K-S Statistic | K-S p-value | Normal? |
|--------|---------------|-------------|--------|
| peak_mag | 0.0485 | 1.7575e-02 | No |
| snr | 0.0408 | 6.9554e-02 | Yes |
| phase_coherence | 0.1977 | 1.0032e-34 | No |
| band_energy | 0.1236 | 9.0056e-14 | No |

*Note: Non-parametric Wilcoxon test used regardless of normality.*

## Analysis Parameters

```yaml
dataset: brunello_parsed.fasta
seed: 42
temperature_k: 310.150000
na_concentration_m: 1.000000
at_lifetime: 1.000000
gc_lifetime: 50.000000
helical_period: 10.500000
apply_helical: True
target_frequency: 0.095238
band_width: 0.010000
noise_band_width: 0.030000
czt_points: 256
gc_threshold: 0.500000
max_sequences: 2000
max_per_group: 500
n_bootstrap: 500
confidence_level: 0.950000
analysis_timestamp: 2025-12-05T01:44:41.076323+00:00
```

## Reproducibility

- **Random Seed:** 42
- **Bootstrap Iterations:** 500
- **Confidence Level:** 95%
- **Temperature:** 310.15 K (37°C physiological)

## Methodology

1. **Sequence Encoding:** Complex-valued signal with real part = breathing lifetimes (AT=1.0 ms, GC=50.0 ms) and imaginary part = SantaLucia nearest-neighbor thermodynamic stability.

2. **Spectral Analysis:** Chirp Z-Transform (CZT) focused on 1/10.5 bp⁻¹ helical frequency with helical phase modulation.

3. **Effect Size:** Hedges' g (bias-corrected Cohen's d) with 500-iteration bootstrap for 95% confidence interval.

4. **Significance Testing:** Two-sided Wilcoxon rank-sum test with Benjamini-Hochberg FDR correction for multiple comparisons.
