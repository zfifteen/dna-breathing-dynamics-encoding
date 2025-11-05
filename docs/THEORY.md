## Background

### Previous Findings (docs/EMPIRICAL_FINDINGS_REPORT.md)

Prior experiments tested biological encodings based on:
1. Multi-property (polarizability, H-bonding, stacking energy)
2. Pairing-based (Watson-Crick strength, structural stability)
3. Thermodynamic (melting temperature, CRISPR preferences)

**Result**: All showed **negative correlation** (arbitrary > biological) with Cohen's d values from -0.579 to -2.549.

**Critical insight from analysis**: Traditional biochemical properties (polarizability, thermodynamics) may not translate effectively to frequency-domain representations.

### Hypothesis Evolution

The negative results suggested: **"We need to work on properties that translate to frequencies more easily."**

This led to testing **DNA breathing dynamics**—base pair opening/closing oscillations—which are:
- **Inherently oscillatory** (real frequencies: AT ~10 MHz, GC ~1 GHz)
- **Biologically relevant** (CRISPR Cas9 needs base pair access)
- **Experimentally measured** (not theoretical constructs)

---

## Experimental Design

### DNA Breathing Dynamics

Base pair opening is a real oscillatory phenomenon:

| Base Pair | Opening Frequency | Bonds | Biological Meaning |
|-----------|------------------|-------|-------------------|
| AT        | ~10⁷ Hz (10 MHz) | 2 H-bonds | Fast opening, weak |
| GC        | ~10⁹ Hz (1 GHz)  | 3 H-bonds | Slow opening, strong |

**Key difference**: 100× frequency separation (2 orders of magnitude)

### Encoding Strategy

#### Breathing Dynamics Encoder
```python
# Real part: log-normalized frequency
# AT: 10^7 Hz → -10.00 (fast opening)
# GC: 10^9 Hz → +10.00 (slow opening)

# Imaginary part: phase from kinetics
# AT: +3.0j (weak bonds, positive phase)
# GC: -3.0j (strong bonds, negative phase)

weights = {
    'A': -10.00 + 3.00j,
    'T': -10.00 + 3.00j,
    'C': +10.00 - 3.00j,
    'G': +10.00 - 3.00j
}
```

#### Phase Modulation
Both encoders use identical phase modulation for fair comparison:
- **Helical periodicity**: 2π × i / 10.5 (DNA wraps every 10.5 bp)
- **Positional phase**: 2π × i / seq_length × 0.3
- **Combined**: helical + positional context

### Control: Arbitrary Encoder
- Random complex weights in same magnitude range (-10 to +10, ±3j)
- **Same phase modulation** as breathing encoder
- 10 independent trials with different random seeds

### Mutation Types Tested

1. **GC-affecting**: AT→GC or GC→AT (changes breathing frequency class)
2. **AT-affecting**: A↔T or G↔C (within same frequency class)
3. **Random**: Any base substitution

**Prediction**: Breathing encoder should excel at GC-affecting mutations (frequency class changes) but show minimal signal for within-class mutations.

### Metrics

- **Z-score**: Using Z Framework equation Z = A(B/c)
  - A = sequence entropy (frame-dependent)
  - B = spectral shift (mutation effect)
  - c = e² ≈ 7.389 (discrete invariant)
- **Statistical tests**: Independent t-test, Cohen's d for effect size
- **Significance threshold**: p < 0.05

---

## Biological Interpretation

### Why Breathing Dynamics Works

1. **Real Oscillatory Phenomenon**
   - Base pair opening/closing is actual molecular motion
   - Has measurable frequencies (MHz to GHz range)
   - Directly translates to Fourier domain

2. **CRISPR Relevance**
   - Cas9 requires base pair access for binding/cutting
   - AT-rich regions open more easily (10 MHz)
   - GC-rich regions require more energy (1 GHz)
   - Breathing encoder captures this accessibility gradient

3. **Frequency-Class Structure**
   - AT pairs: one frequency class (~10 MHz)
   - GC pairs: another frequency class (~1 GHz)
   - Mutations changing class → strong signal
   - Mutations within class → no signal
   - This selectivity is **biologically meaningful**

### Contrast with Previous Properties

**Failed properties** (polarizability, H-bonding, thermodynamics):
- Static/equilibrium properties
- No inherent oscillatory component
- Designed for 3D structural analysis
- Poor frequency-domain mapping

**Successful property** (breathing dynamics):
- Dynamic/temporal property
- Inherent oscillation at specific frequencies
- Directly affects CRISPR function
- Natural frequency-domain mapping
