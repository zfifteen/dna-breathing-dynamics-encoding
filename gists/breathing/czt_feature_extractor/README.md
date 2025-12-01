# DNA Breathing Dynamics CZT Feature Extractor Gist

This gist implements a biophysical encoding of DNA sequences into complex signals using base-pair breathing lifetimes (real part) and thermodynamic stability (imag part from SantaLucia-inspired ΔG), with optional helical phase modulation (2π/10.5 bp). It performs Chirp Z-Transform (CZT) zoomed to the helical harmonic (~1/10.5 cycles/bp) using scipy.signal.czt, extracts features (peak magnitude/phase, peak-to-skirt ratio, phase coherence, band energy, SNR, QC flag), generates controls (dinuc-approximate shuffles, phase scramble, label permutations), and computes statistics (Cohen's d with bootstrap CI, permutation p-values).

## Requirements
- Python 3.12+
- NumPy
- SciPy (>=1.12 for signal.czt; pip install scipy)

No external data; uses user-provided real sequences (FASTA/TSV/CSV).

## Usage
```
python dna_breathing_gist.py --input sequences.fasta --output features.csv --seed 42 [options]
```

### Key Options
- `--input`: Input file (FASTA: >header\nseq; TSV/CSV: seq,label columns)
- `--output`: Output CSV with features per seq + stats row
- `--groups groupA,groupB`: Groups for Cohen's d (if labels present)
- `--seed 42`: RNG seed for reproducibility
- `--at-lifetime 5 --gc-lifetime 25`: Lifetimes (ms)
- `--helical-period 10.5 --apply-helical`: Helical modulation (default on)
- `--iupac-policy average`: 'average' or 'mask' for ambiguity
- `--band-width 0.01`: CZT zoom band (± around 1/10.5)
- `--apply-taper --taper-k 0.3`: Golden-ratio taper
- `--num-shuffles 100`: Number of shuffles for null
- `--bootstrap-samples 1000 --num-perm 100`: For CI and p-values

### Example
Assume `sample.fasta`:
```
>seq1
ATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTA
```

Run:
```
python dna_breathing_gist.py --input sample.fasta --output test.csv --seed 42 --num-shuffles 5 --bootstrap-samples 10 --num-perm 10
```

Outputs `test.csv` with features (peak_mag, etc.) and stats row (cohens_d, p_perm). Prints summary.

## Features
- Peak magnitude/index/freq/phase in helical band
- Peak-to-skirt ratio (off-band mean)
- Phase coherence (|mean(exp(iφ))| in band)
- Band energy (sum |s|^2)
- SNR (band/off-band mean)
- QC flag (peak near edge)

## Controls & Stats
- Shuffles: Base-count preserving (seeded)
- Phase scramble: Random phases on spectrum
- Permutation test on primary feature (peak_mag)
- Cohen's d with 95% bootstrap CI

## Reproducibility
Deterministic with seed; <5s on small sets (e.g., 10 seqs ~20bp). Modular pipeline.

For full validation, use real CRISPR guides with GC-affecting mutations to test selectivity.

See code docstrings for details.