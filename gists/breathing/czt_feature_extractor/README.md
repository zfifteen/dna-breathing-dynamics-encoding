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

## Methodology

### Dinucleotide Shuffling Strategy
To preserve the thermodynamic stability profile of DNA sequences (driven largely by nearest-neighbor interactions), this tool uses a **graph-theoretic approach** for shuffling:

1.  **Eulerian Path Analysis**: The sequence is modeled as a directed multigraph where nodes are nucleotides and edges are transitions.
2.  **Strict Validation**: Before shuffling, the graph is checked for **Eulerian conditions** (balanced in/out degrees).
    *   **Behavior Change**: Non-Eulerian inputs (often caused by gaps or invalid characters) now **raise a `ValueError`** with detailed imbalance reports. They do *not* fallback to mononucleotide shuffles or original sequences, preventing statistical artifacts in spectral baselines.
3.  **Reproducibility**: Valid shuffles are generated deterministically using the provided `--seed`.

Original docstrings/details unchanged.