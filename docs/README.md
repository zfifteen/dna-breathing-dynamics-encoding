# DNA Breathing Dynamics Encoding - Documentation

Complete documentation for DNA Breathing Dynamics Encoding: using biophysically-grounded base-pair opening rates with fractional-period spectral analysis for CRISPR activity prediction.

## Quick Start

**New to this project?**

1. **[OVERVIEW.md](OVERVIEW.md)** - 5-minute summary of breakthrough findings
2. **[CHARTER.md](../CHARTER.md)** (repo root) - Project mission and rationale
3. **[REPRODUCTION.md](REPRODUCTION.md)** - Reproduce the key results

## Documentation Index

### Core Documentation
- **[OVERVIEW.md](OVERVIEW.md)** - Executive summary and key breakthrough (Cohen's d = +4.130)
- **[THEORY.md](THEORY.md)** - Mathematical foundations and biological basis
- **[IMPLEMENTATION.md](IMPLEMENTATION.md)** - Complete implementation guide (encoders, CZT, ablations)
- **[VALIDATION.md](VALIDATION.md)** - Statistical methodology and scientific rigor
- **[RESULTS.md](RESULTS.md)** - All experimental results and benchmarks
- **[REPRODUCTION.md](REPRODUCTION.md)** - Step-by-step reproduction instructions
- **[DATASETS.md](DATASETS.md)** - Data sources and provenance

## Key Findings

> **Breakthrough (Jan 2025)**: DNA breathing dynamics encoding—based on experimental base pair opening frequencies (AT: 10 MHz, GC: 1 GHz)—achieves Cohen's d = +4.130 over arbitrary encodings for GC-content-affecting mutations. This is the **first biological property** to outperform arbitrary weights in spectral DNA analysis, validating the hypothesis that frequency-native properties are essential for effective spectral encoding. See [RESULTS.md](RESULTS.md) for full results.

## Key Concepts

- **DNA Breathing**: Base-pair opening rates (AT: 1ms, GC: 50ms) from PMC5393899
- **Helical Periodicity**: 10.5 bp/turn requiring CZT for fractional-period analysis
- **Complex Encoding**: Real (kinetics) + Imaginary (thermodynamics) components
- **Statistical Rigor**: Bootstrap CIs, permutation tests, FDR correction

## Reading Paths

### Quick (15 minutes)
```
OVERVIEW.md → RESULTS.md → REPRODUCTION.md
```

### Complete (60 minutes)
```
OVERVIEW.md → THEORY.md → IMPLEMENTATION.md → VALIDATION.md → RESULTS.md
```

### Research Focus
```
THEORY.md → VALIDATION.md → RESULTS.md → DATASETS.md
```

## Scientific Gates

All work satisfies 8 mandatory gates (from CHARTER.md):
1. Human DNA only (A/C/G/T/N)
2. No fabrication (real sequences)
3. Fail-fast validation
4. Biophysical anchoring
5. Dimensionless encoding
6. Geometric resolution (θ′(n,k) with k≈0.3)
7. Discrete domain
8. Full reproducibility

---

**Project Status**: Complete validation package
**Confidence Level**: 10/10 (per CHARTER.md)
**Last Updated**: 2025-11-05
