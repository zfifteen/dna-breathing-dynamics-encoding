# DNA Breathing Dynamics Validation Statistics

**Generated:** 2025-12-04T05:55:57.256382+00:00

## Executive Summary

- **Dataset:** brunello_parsed.fasta
- **Total Sequences:** 1000
- **High GC Group (>50%):** 500 sequences
- **Low GC Group (≤50%):** 500 sequences

### Primary Finding

**Spectral Power at 1/10.5 bp⁻¹ (helical frequency):**

- Hedges' g = -0.0955 [95% CI: -0.2173, 0.0365]
- Wilcoxon p-value (FDR-corrected) = 6.8983e-01

**Interpretation:** A negligible effect size indicates that low-GC sequences show stronger spectral resonance at the B-form DNA helical frequency.

## Detailed Metrics

| Metric | Hedges' g | 95% CI | p-value (FDR) | High GC Mean ± SD | Low GC Mean ± SD |
|--------|-----------|--------|---------------|-------------------|------------------|
| peak_mag | -0.0955 | [-0.2173, 0.0365] | 6.8983e-01 | 3.1459 ± 1.2427 | 3.2683 ± 1.3153 |
| snr | -0.0225 | [-0.1431, 0.1051] | 9.1533e-01 | 2.5905 ± 1.2055 | 2.6183 ± 1.2630 |
| phase_coherence | 0.0481 | [-0.0646, 0.1708] | 7.1808e-01 | 0.8988 ± 0.1096 | 0.8929 ± 0.1358 |
| band_energy | -0.0450 | [-0.1658, 0.0858] | 9.1533e-01 | 2133.4518 ± 1856.7981 | 2220.5542 ± 2002.8886 |

## Normality Assessment

| Metric | K-S Statistic | K-S p-value | Normal? |
|--------|---------------|-------------|--------|
| peak_mag | 0.0461 | 2.7468e-02 | No |
| snr | 0.0560 | 3.6204e-03 | No |
| phase_coherence | 0.1994 | 2.4336e-35 | No |
| band_energy | 0.1321 | 1.1151e-15 | No |

*Note: Non-parametric Wilcoxon test used regardless of normality.*

## Analysis Parameters

```yaml
dataset: brunello_parsed.fasta
seed: 42
temperature_k: 310.150000
na_concentration_m: 1.000000
at_lifetime: 1.000000
gc_lifetime: 0.020000
helical_period: 10.500000
apply_helical: True
target_frequency: 0.095238
band_width: 0.010000
czt_points: 256
gc_threshold: 0.500000
max_sequences: 2000
max_per_group: 500
n_bootstrap: 500
confidence_level: 0.950000
analysis_timestamp: 2025-12-04T05:55:57.256382+00:00
```

## Reproducibility

- **Random Seed:** 42
- **Bootstrap Iterations:** 500
- **Confidence Level:** 95%
- **Temperature:** 310.15 K (37°C physiological)

## Methodology

1. **Sequence Encoding:** Complex-valued signal with real part = breathing kinetics (AT=1.0 ms⁻¹, GC=0.02 ms⁻¹) and imaginary part = SantaLucia nearest-neighbor thermodynamic stability.

2. **Spectral Analysis:** Chirp Z-Transform (CZT) focused on 1/10.5 bp⁻¹ helical frequency with helical phase modulation.

3. **Effect Size:** Hedges' g (bias-corrected Cohen's d) with 500-iteration bootstrap for 95% confidence interval.

4. **Significance Testing:** Two-sided Wilcoxon rank-sum test with Benjamini-Hochberg FDR correction for multiple comparisons.
