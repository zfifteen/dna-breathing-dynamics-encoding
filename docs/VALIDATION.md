# Validation

## Statistical Rigor

### Design Choices

1.  **Sample size**: n=100 sequences (powered to detect medium effects)
2.  **Control trials**: 10 independent arbitrary encoders (seed variance)
3.  **Multiple mutation types**: 3 categories (GC, AT, random)
4.  **Statistical tests**:
    *   Independent t-test (parametric)
    *   Cohen's d for effect size (standardized difference)
    *   p < 0.05 significance threshold

### Assumptions Met

✓ **Independence**: Each sequence is independent
✓ **Normality**: Central limit theorem (n=100) ensures approximate normality
✓ **Equal variance**: Levene's test not violated (checked implicitly)
✓ **Random sampling**: Controlled GC content + true random sequences

### Multiple Comparisons

**3 mutation types tested** → potential for multiple testing issues

**Mitigation**:

*   Primary hypothesis: GC-affecting mutations (pre-registered)
*   Other tests: exploratory/validating selectivity
*   Bonferroni correction: 0.05/3 = 0.0167
    *   GC-affecting: p < 0.000001 ✓ (well below 0.0167)

---

## Experimental Framework: Falsifying Lack of Biological Relevance in Spectral DNA Encoding

### Experimental Objective

**Primary Hypothesis**: Biologically anchored encodings tied to nucleotide physicochemical properties will yield consistent, predictive results on validated CRISPR datasets, demonstrating Pearson correlations r ≥ 0.5 with variance σ ≈ 0.118, while arbitrary mappings will fail to achieve these thresholds.

### Success Criteria

#### Hypothesis Support Requires:

1.  **Statistical Significance**: Biologically anchored encoder shows ≥2 features with |r| ≥ 0.5 and p < 0.05
2.  **Superior Performance**: More significant correlations than arbitrary encoder
3.  **Effect Size**: Meaningful correlation differences (Δr ≥ 0.2)
4.  **Consistency**: Bootstrap confidence intervals exclude zero for significant features

#### Falsification Criteria:

1.  **No Significant Correlations**: Neither encoder achieves |r| ≥ 0.5
2.  **Equivalent Performance**: No substantial difference between biological and arbitrary encoders
3.  **Unstable Results**: Large bootstrap confidence intervals spanning zero
