# Scenario 1: GC-Bias in DNA Breathing Features (Scaled Validation)

## Methods
- Data: 1000 randomly sampled sequences (500 high GC >50%, 500 low GC <40%) from the full Brunello library, len=20 ATGC.
- Groups: high_gc vs low_gc (labeled in FASTA).
- Run: `python src/dna_breathing_gist.py --input data/processed/brunello_gc1000.fasta --output results/scenarios/scenario1_gc_bias/brunello_gc1000_features.csv --num-shuffles 50 --bootstrap 1000 --num-perm 500 --seed 42`
- Output: `results/scenarios/scenario1_gc_bias/brunello_gc1000_features.csv`, `results/scenarios/scenario1_gc_bias/run_log_gc1000.txt`
- Plots: (To be generated: hists for peak_mag/snr/coherence)

## Results
| Metric | Cohen's d | CI Low | CI High | p_perm |
|--------|-----------|--------|---------|--------|
| peak_mag | -0.1833 | -0.2960 | -0.0640 | 0.0020 |
| snr | -0.0261 | -0.1449 | 0.1100 | 0.6920 |
| phase_coherence | 0.1461 | 0.0274 | 0.2714 | 0.0120 |

- Interpretation: GC bias persists for Peak Mag and Coherence but with much smaller effect size than the initial small-sample study. The large effect in n=20 was likely overestimated or due to sampling bias. The negative d for Peak Mag confirms high GC suppresses breathing amplitude.

## Plots
- (Placeholder for plots to be generated based on n=1000 data)

## Conclusion
Significant biophysical GC-bias in CZT features (Peak Mag and Coherence) is validated with a larger dataset. The effect sizes are smaller (d=-0.18 for Peak Mag) compared to initial small-scale observations, highlighting the importance of larger sample sizes for robust effect estimation. Further investigation into the difference in effect sizes between small and large samples would be beneficial.