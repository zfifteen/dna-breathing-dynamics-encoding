# Integration Plan for DNA Breathing Dynamics CRISPR Tool

## Overview

I've forked both FlashFry and CRISPOR repositories to integrate my DNA breathing dynamics scoring method. This document outlines the roadmap for implementing and validating the breathing dynamics features in these established CRISPR design tools.

## Phase 1: Quick Win with FlashFry (This Week)

### Create External Scorer Wrapper

- **Script**: `integrations/flashfry/add_breathing_score.py`
- **Functionality**:
  - Reads FlashFry tab-delimited output
  - Computes breathing dynamics score for each guide
  - Outputs new column with scores

### Test on Brunello Dataset

- Index Brunello sequences in FlashFry
- Run integrated pipeline
- Compare breathing scores vs existing FlashFry metrics

## Phase 2: CRISPOR Integration (Next Month)

### Modify crisporWebsite Fork

- **crisporEffScores.py**: Add `scoreBreathingDynamics()` function
- **crispor.py**: Call function in main pipeline
- **crispor.html**: Add column for breathing score in output

## Phase 3: Validation Study (3-6 Months)

### With Both Integrations

1. Generate predictions for 50 guides:
   - 20 high-risk seed GC changes
   - 20 low-risk PAM changes
   - 10 controls
2. Partner with CRISPR lab for wet-lab testing
3. Measure correlation between predictions and actual off-targets
4. Publish results or bioRxiv preprint

## Immediate To-Do

1. **Create add_breathing_score.py** (2-4 hours)
2. **Test on 10 sample sequences**
3. **Document usage**
4. **Run full pipeline on Brunello subset**
5. **Reach out to CRISPR labs for validation**

## Key Findings to Validate

- **Seed region sensitivity**: Cohen's d > 2.0
- **PAM tolerance**: d ~0.7-1.1
- **Off-target reduction**: 20-30% reduction in off-target predictions

## Repositories Ready

- `zfifteen/FlashFry` (forked)
- `zfifteen/crisporWebsite` (forked)
- `zfifteen/dna-breathing-dynamics-encoding` (main)

## Technical Approach

### Integration with wave-crispr-signal Framework

The implementation creates a unified integration that bridges the DNA breathing dynamics encoding with the wave-crispr-signal spectral framework:

1. **Read FlashFry TSV output**
2. **Compute Z-score** using `PhaseWeightedScorecard` from wave-crispr-signal
3. **Optionally compute CZT-based breathing features** from dna-breathing-dynamics
4. **Output augmented TSV** with new columns

### Scoring Modes

- **z_score**: Phase-weighted Z-invariant only (fast, ~0.5ms/guide)
- **breathing**: CZT breathing features only (~2ms/guide)  
- **both**: All features (~2.5ms/guide)

### Validation Metrics (from Grok collaboration)

**Round 1 Technical Validation:**
- Start with phase-weighted Z-score first (fast, ~0.5ms/guide)
- Track: ΔAUROC/AUPRC vs FlashFry baseline
- Cohen's d on GC-mutated subsamples
- Spearman r between breathing vs efficiency
- Permutation p-values
- False pos/neg rates in low-GC tails

**Test Sequence Stratification:**
- 2 low-GC (<40%)
- 3 mid (40-60%)
- 2 high (>60%)
- Vary PAM contexts (NGG, NGA)
- Include 3 with known off-targets

**Round 2 Ablation Study Design:**

For Brunello subsample (n=500):

1. **Metric Computation Pipeline:**
   - Load Brunello subset with GC stratification
   - Compute both phase-weighted Z-scores and full CZT features
   - Calculate ΔAUROC/AUPRC vs FlashFry baseline
   - Compute Cohen's d on GC-mutated subsamples
   - Spearman r correlations
   - Permutation testing framework

2. **GC Stratification Logic:**
   - Bin n=500 guides into low/mid/high GC bins
   - Ensure balanced representation
   - Edge case handling for boundary sequences

3. **Comparison Framework:**
   - Fair evaluation of phase-only vs full-CZT performance
   - Thresholds: full CZT features worth 4x compute cost if ΔAUROC>0.05
   - Statistical power: n=500 sufficient for d>0.3 at 80% power

## Next Step

Draft `add_breathing_score.py` implementation

---

*Generated: January 7, 2026*  
*Status: Planning Phase*  
*Collaboration: Perplexity Space + Grok technical validation*
