# Testing Guidelines for DNA Breathing Dynamics CZT Feature Extractor

## Core Requirement: Real Data Only
- **No Synthetic Data**: All tests, scenarios, and validations must use real, publicly available datasets (e.g., Brunello or GeCKO CRISPR sgRNA libraries from Addgene/DepMap/NCBI). Synthetic, simulated, or generated sequences (e.g., random AT/GC patterns, mutated placeholders) are strictly invalid and must be discarded immediately.
- **Enforcement**: Always refuse any plan, execution, or tool call involving synthetic data. If a fetch fails or data is unavailable, halt and request manual real data (e.g., "Download from Addgene and provide TSV"). Prioritize public sources; validate seqs (len=20bp, ATGC only, biological context like hg38 CRISPR spacers).
- **Rationale**: Synthetic data undermines biophysical validity (e.g., artificial GC bias misses real guide dynamics). Real data ensures relevance for CRISPR applications (GC affects stability/breathing).

## Validation Framework (AC1–AC6 Integration)
- **AC1 (Encoding)**: Assert real_part means differ by group (e.g., low-GC < high-GC due to lifetimes 5ms vs 25ms); no IUPAC unless in source.
- **AC2 (CZT)**: Peak_freq ≈0.0952 Hz (1/10.5); band within ±0.01.
- **AC3 (Features)**: SNR/peak_mag >20 in real; QC_flag=1 for >95% seqs.
- **AC4 (Controls)**: Shuffles/scrambles ↓ peak_mag >10%; p<0.01 vs orig (KS test).
- **AC5 (Stats)**: Bootstrap CI excludes 0; power ≥0.8 (tt_ind_solve_power).
- **AC6 (Z-Hypothesis)**: Corr(coherence, GC%) |r|>0.4 (negative expected).

## Scaling & Significance
- **n Minimum**: ≥300/group (power~0.8 for d=0.5, α=0.05); balance groups (e.g., high/low GC).
- **Stats**: Cohen's d≥0.5, p<0.05 (perm test); FDR (Benjamini-Hochberg) for multiples (peak_mag, SNR, coherence).
- **Params**: --num-shuffles=50+, --bootstrap=1000, --num-perm=500; --seed 42 for repro.

## Workflow
1. **Source Real Data**: Tools (webfetch/task) for TSV/FASTA; manual if fails (e.g., Addgene download).
2. **Prep**: Subsample balanced (awk/pandas); save data/processed/[name].fasta.
3. **Power Check**: src/power_analysis.py (if exists) or sim; confirm n sufficient.
4. **Run**: python src/dna_breathing_gist.py ... | tee results/scenarios/[name]/log.txt
5. **Validate**: src/plot_results.py for PNGs/FDR; assert ACs in REPL (e.g., pd.read_csv assert p<0.05).
6. **Report**: Update reports/scenarios/[name].md with tables/plots; cross in analysis/large_overview.md.
7. **Tests**: pytest src/tests/ (add unit for ACs); lint black/ruff/mypy src/.

## Non-Compliance
- Discard invalid outputs (e.g., rm synthetic CSVs).
- Log refusal: "Refused: Synthetic data invalid per guidelines."

Last Updated: 2025-12-02