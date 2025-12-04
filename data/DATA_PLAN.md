# Data Plan

This document captures the operational plan for the `data/` directory in this repository.

## Quick Start

Run the data download script to obtain, parse, curate, and validate all datasets:

```bash
./data/download_data.sh
```

This will:
1. Create `data/raw/` directory (gitignored)
2. Download raw Brunello library data (~9MB TSV)
3. Parse to FASTA format (~77,000 sequences)
4. Curate and subsample to 1000 sequences (seed=42)
5. Create a minimal seed/sample dataset for testing
6. Validate all datasets

## Goals
- Seed `data/` with curated, provenance-tracked subsets of real biological sequence datasets.
- Commit only small curated artifacts (<= 1000 sequences, <= 5 MB per file).
- Keep raw archives outside git under `data/raw/` (gitignored).
- Provide reproducible scripts to download, curate, and validate datasets.

## Constraints
- Max committed file size: 5 MB.
- Fixed sample size per dataset: 1000 sequences (or fewer if the dataset has fewer than 1000 available sequences after filtering).
- Sequence requirements: 20nt length, ACGT bases only.
- Single repository license: MIT (applies to code, scripts, and curated artifacts).

## Directory Structure

```
data/
├── .gitignore              # Ignores raw/ directory
├── DATA_PLAN.md            # This file
├── README.md               # Data usage documentation
├── datasets.yml            # Dataset catalog with metadata
├── download_data.sh        # Main script to obtain all data
├── raw/                    # Raw downloads (gitignored)
│   ├── broadgpp-brunello-library-contents.txt
│   └── brunello_parsed.fasta
├── human/
│   └── brunello/
│       ├── METADATA.md
│       ├── DOWNLOAD_INSTRUCTIONS.md
│       └── sequences.fasta # 1000 curated sgRNAs
└── seed/
    └── sample/
        ├── METADATA.md
        └── sequences.fasta # 3 sgRNAs for CI testing
```

## Execution Steps

### Automated (Recommended)

```bash
# From repository root
./data/download_data.sh
```

### Manual Steps

If you need to run steps manually:

#### Step 1: Download raw Brunello data

```bash
mkdir -p data/raw
curl -fSL -A "dna-breathing-dynamics-encoding/data-downloader" \
    -o data/raw/broadgpp-brunello-library-contents.txt \
    "https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt"
```

**Fallback:** If download fails, copy from repository:
```bash
cp gists/breathing/czt_feature_extractor/data/processed/brunello.fasta data/raw/brunello_parsed.fasta
```

#### Step 2: Parse raw TSV to FASTA

The raw file is a TSV (tab-separated values) with columns including:
- Target Gene ID
- Target Gene Symbol  
- sgRNA Target Sequence (20nt)
- Target Context Sequence (30nt)

The file uses carriage returns (`\r`) as line separators. Parse it using:

```bash
python3 -c "
import csv
from pathlib import Path

raw = Path('data/raw/broadgpp-brunello-library-contents.txt')
fasta = Path('data/raw/brunello_parsed.fasta')

content = raw.read_text().replace('\r\n', '\n').replace('\r', '\n')
lines = content.strip().split('\n')
reader = csv.reader(lines, delimiter='\t')
header = next(reader)

# Find column indices
seq_idx = header.index('sgRNA Target Sequence')  # Column 6
gene_idx = header.index('Target Gene Symbol')     # Column 1
id_idx = header.index('Target Gene ID')           # Column 0

seqs = []
for row in reader:
    if len(row) > seq_idx:
        seq = row[seq_idx].strip().upper()
        if len(seq) == 20 and all(c in 'ATGC' for c in seq):
            gene = row[gene_idx] if len(row) > gene_idx else 'Unknown'
            gene_id = row[id_idx] if len(row) > id_idx else 'Unknown'
            seqs.append((gene.replace(' ', '_'), gene_id, seq))

with open(fasta, 'w') as f:
    for gene, gene_id, seq in seqs:
        f.write(f'>{gene}|{gene_id}\n{seq}\n')

print(f'Parsed {len(seqs)} sequences')
"
```

#### Step 3: Curate and subsample

Use the canonical curation script:

```bash
python scripts/curate_and_subsample.py --dataset human/brunello --max-seqs 1000 --seed 42
```

#### Step 4: Create seed/sample dataset

```bash
head -6 data/human/brunello/sequences.fasta > data/seed/sample/sequences.fasta
```

#### Step 5: Validate

```bash
python scripts/validate_dataset.py --dataset human/brunello
python scripts/validate_dataset.py --dataset seed/sample
```

## Dataset Details

### Dataset 1: human/brunello

- **Name:** Brunello sgRNA library (curated subset)
- **Source:** [Addgene Brunello Library](https://www.addgene.org/pooled-library/broadgpp-brunello/)
- **Paper:** Doench et al. (2016) Nature Biotechnology
- **DOI:** 10.1038/nbt.3437
- **Raw URL:** `https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt`
- **Raw size:** ~9 MB (77,441 sgRNAs after parsing)
- **Curated:** 1000 sequences (seed=42, 20nt, ACGT only)
- **Output:** `data/human/brunello/sequences.fasta` (~35KB)
- **License:** See source (Addgene / paper supplement)

### Dataset 2: seed/sample

- **Name:** Minimal test dataset
- **Source:** Derived from human/brunello
- **Size:** 3 sequences
- **Purpose:** CI testing and script validation
- **Output:** `data/seed/sample/sequences.fasta`
- **License:** MIT

## Planned Datasets (Future)

The following datasets are defined in `datasets.yml` as planned:

1. **human/doench2016** — Guide-efficiency dataset from same paper
2. **human/depmap_subsample** — DepMap-derived sample (CC-BY 4.0)
3. **mouse/gecko_v2** — GeCKO v2 mouse library

## Reproducibility

- **Random seed:** All subsampling uses seed=42 for reproducibility
- **Commands logged:** Exact curation commands are recorded in METADATA.md
- **SHA256 checksums:** Stored in METADATA.md for verification
- **Validation:** Run `scripts/validate_dataset.py` to verify constraints

## Validation Criteria

All curated datasets must pass:

1. ✅ Maximum 1000 sequences
2. ✅ Maximum 5 MB file size
3. ✅ All sequences exactly 20 nucleotides
4. ✅ Only A/C/G/T bases (no ambiguity codes)
5. ✅ No duplicate sequences
6. ✅ Valid FASTA format

## Troubleshooting

### Download fails with 403

Addgene may block automated downloads. Solutions:
1. Use the repository fallback: `cp gists/.../brunello.fasta data/raw/`
2. Download manually in a browser and place in `data/raw/`

### Validation fails

Check the error message from `validate_dataset.py`:
- "Too many sequences": Increase `--max-seqs` or check subsampling
- "Invalid bases": Source data may have ambiguity codes; use `--relaxed`
- "Wrong length": Check `--seq-length` parameter

### Script permissions

```bash
chmod +x data/download_data.sh
```