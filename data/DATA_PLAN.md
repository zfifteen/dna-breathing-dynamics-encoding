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

The raw file is a tab-separated format with (sequence, gene_name) pairs. Parse it:

```bash
python3 -c "
import sys
from pathlib import Path

raw = Path('data/raw/broadgpp-brunello-library-contents.txt')
fasta = Path('data/raw/brunello_parsed.fasta')

content = raw.read_text().strip().split('\t')
seqs = []
for i in range(0, len(content), 2):
    seq = content[i].upper()
    gene = content[i+1] if i+1 < len(content) else 'Unknown'
    if len(seq) == 20 and all(c in 'ACGT' for c in seq):
        seqs.append((gene.replace(' ', '_'), seq))

with open(fasta, 'w') as f:
    for i, (gene, seq) in enumerate(seqs, 1):
        f.write(f'>{gene}|{i}\n{seq}\n')

print(f'Parsed {len(seqs)} sequences')
"
```

#### Step 3: Curate and subsample

```bash
python scripts/curate_and_subsample.py --dataset human/brunello --max-seqs 1000 --seed 42
```

Or manually:

```bash
python3 -c "
import random
from pathlib import Path

random.seed(42)
fasta_in = Path('data/raw/brunello_parsed.fasta')
fasta_out = Path('data/human/brunello/sequences.fasta')

# Read sequences
records = []
header, seq = None, []
for line in fasta_in.read_text().splitlines():
    if line.startswith('>'):
        if header: records.append((header[1:], ''.join(seq)))
        header, seq = line, []
    else: seq.append(line)
if header: records.append((header[1:], ''.join(seq)))

# Filter and subsample
valid = [(h, s.upper()) for h, s in records if len(s) == 20 and all(c in 'ACGT' for c in s)]
if len(valid) > 1000: valid = random.sample(valid, 1000)

fasta_out.parent.mkdir(parents=True, exist_ok=True)
fasta_out.write_text('\n'.join(f'>{h}\n{s}' for h, s in valid) + '\n')
print(f'Wrote {len(valid)} sequences')
"
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