# DNA Breathing Dynamics CZT Feature Extractor Gist

See original README content here, but with updated folder structure.

## Folder Structure
```
czt_feature_extractor/
├── README.md
├── TESTING_GUIDELINES.md  # Strict real-data policy
├── data/
│   ├── raw/
│   │   └── Brunello_library_v1.1_sequences.txt.zip
│   ├── processed/
│   │   ├── sample.fasta
│   │   ├── real_brunello.fasta
│   │   ├── real_crispr.fasta
│   │   └── mutated_brunello.fasta
│   └── archive/
│       └── brunello_sample.fasta
├── results/
│   ├── scenarios/
│   │   ├── scenario1_gc_bias/
│   │   │   └── brunello_features.csv
│   │   ├── scenario2_mutational/
│   │   │   └── mutational_features.csv
│   │   └── scenario3_controls/
│   │       └── control_features.csv
│   ├── archive/
│   │   ├── features.csv
│   │   ├── features_test.csv
│   │   └── real_features.csv
│   └── summary/
├── reports/
│   ├── scenarios/
│   │   ├── scenario1_gc_bias.md
│   │   ├── scenario2_mutational.md
│   │   └── scenario3_controls.md
│   ├── analysis/
│   └── archive/
└── src/
    └── dna_breathing_gist.py
```

## Usage
Run from root: `python src/dna_breathing_gist.py --input data/processed/real_brunello.fasta --output results/scenarios/scenario1_gc_bias/output.csv [options]`

For large tests, place inputs in data/processed/, outputs in results/scenarios/[name]/.

**Testing Policy**: See TESTING_GUIDELINES.md – real data only, no synthetic. All tests must validate AC1–AC6.

Original docstrings/details unchanged.