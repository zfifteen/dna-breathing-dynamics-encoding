#!/usr/bin/env bash
#
# download_data.sh
#
# This script executes the data plan to download, extract, parse, curate, 
# and validate all datasets for the dna-breathing-dynamics-encoding project.
#
# Usage:
#   ./data/download_data.sh
#
# Prerequisites:
#   - Python 3.8+ with standard library
#   - curl or wget for downloads
#   - Internet connection for downloading raw files
#
# The script will:
#   1. Create necessary directories
#   2. Download raw data files to data/raw/
#   3. Parse raw files to FASTA format
#   4. Curate and subsample to <=1000 sequences per dataset
#   5. Validate the curated datasets
#
# Note: data/raw/ is gitignored. Only curated files in data/<dataset>/ are committed.

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root is the parent of the data directory
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=============================================="
echo "DNA Breathing Dynamics Data Download Script"
echo "=============================================="
echo "Project root: $PROJECT_ROOT"
echo ""

# Change to project root
cd "$PROJECT_ROOT"

# Create required directories
echo "[Step 1] Creating directories..."
mkdir -p data/raw
mkdir -p data/human/brunello
mkdir -p data/seed/sample
echo "  ✓ Directories created"
echo ""

# ==============================================================================
# DATASET 1: human/brunello - Brunello sgRNA Library
# ==============================================================================
echo "[Step 2] Processing human/brunello dataset..."
echo "  Source: Addgene Brunello Library (Doench et al. 2016)"
echo "  DOI: 10.1038/nbt.3437"
echo ""

BRUNELLO_RAW="data/raw/broadgpp-brunello-library-contents.txt"
BRUNELLO_URL="https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt"
BRUNELLO_FASTA="data/raw/brunello_parsed.fasta"

# Check if we already have the raw file in the gists directory (from repo)
REPO_BRUNELLO="gists/breathing/czt_feature_extractor/data/raw/broadgpp-brunello-library-contents.txt"

if [ -f "$REPO_BRUNELLO" ] && [ -s "$REPO_BRUNELLO" ]; then
    echo "  → Using existing Brunello raw file from repository..."
    cp "$REPO_BRUNELLO" "$BRUNELLO_RAW"
elif [ -f "$BRUNELLO_RAW" ] && [ -s "$BRUNELLO_RAW" ]; then
    echo "  → Raw file already exists, skipping download"
else
    echo "  → Downloading Brunello library contents..."
    curl -fSL -A "dna-breathing-dynamics-encoding/data-downloader" \
         -o "$BRUNELLO_RAW" "$BRUNELLO_URL" || {
        echo "  ⚠ Download failed. Attempting from processed FASTA in repo..."
        # Fallback: use pre-processed FASTA from repo if available
        REPO_PROCESSED="gists/breathing/czt_feature_extractor/data/processed/brunello.fasta"
        if [ -f "$REPO_PROCESSED" ]; then
            cp "$REPO_PROCESSED" "$BRUNELLO_FASTA"
            echo "  → Using pre-processed FASTA from $REPO_PROCESSED"
        else
            echo "  ✗ No fallback available. Please download manually."
            exit 1
        fi
    }
fi

# Parse the raw TSV to FASTA if we have the raw file and no parsed FASTA yet
if [ -f "$BRUNELLO_RAW" ] && [ -s "$BRUNELLO_RAW" ] && [ ! -f "$BRUNELLO_FASTA" ]; then
    echo "  → Parsing raw TSV to FASTA format..."
    python3 << 'PYTHON_PARSER'
import sys
from pathlib import Path

raw_path = Path("data/raw/broadgpp-brunello-library-contents.txt")
fasta_path = Path("data/raw/brunello_parsed.fasta")

if not raw_path.exists():
    print("  ✗ Raw file not found")
    sys.exit(1)

# Read and parse the TSV-like file (may be single line with tab-separated values)
with open(raw_path, 'r') as f:
    content = f.read()

# Split by tabs and process pairs of (sequence, gene_name)
parts = content.strip().split('\t')
valid_seqs = []
i = 0
seq_num = 0

while i < len(parts):
    # Each record should have: sequence, then gene info
    seq = parts[i].strip().upper() if i < len(parts) else ""
    gene = parts[i+1].strip() if i+1 < len(parts) else ""
    
    # Handle empty or missing gene names
    if not gene:
        gene = "Unknown"
    
    # Validate: 20bp, ACGT only
    if len(seq) == 20 and all(c in "ACGT" for c in seq):
        seq_num += 1
        header = f">{gene.replace(' ', '_')}|{seq_num}"
        valid_seqs.append((header, seq))
    
    i += 2  # Move to next pair

with open(fasta_path, 'w') as f:
    for header, seq in valid_seqs:
        f.write(f"{header}\n{seq}\n")

print(f"  ✓ Parsed {len(valid_seqs)} valid 20nt sequences to {fasta_path}")
PYTHON_PARSER
fi

# Check if we have the parsed FASTA from repo fallback
if [ ! -f "$BRUNELLO_FASTA" ]; then
    REPO_PROCESSED="gists/breathing/czt_feature_extractor/data/processed/brunello.fasta"
    if [ -f "$REPO_PROCESSED" ]; then
        cp "$REPO_PROCESSED" "$BRUNELLO_FASTA"
        echo "  → Using pre-processed FASTA from repository"
    fi
fi

# Curate and subsample to 1000 sequences
if [ -f "$BRUNELLO_FASTA" ]; then
    echo "  → Curating and subsampling to 1000 sequences..."
    python3 << 'PYTHON_CURATE'
import random
import hashlib
from pathlib import Path

random.seed(42)  # Reproducibility

fasta_in = Path("data/raw/brunello_parsed.fasta")
fasta_out = Path("data/human/brunello/sequences.fasta")

# Read all sequences
records = []
header = None
seq_lines = []
with open(fasta_in, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines)))
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        records.append((header, "".join(seq_lines)))

print(f"  → Read {len(records)} sequences from {fasta_in}")

# Filter: 20nt, ACGT only
filtered = [(h, s.upper()) for h, s in records 
            if len(s) == 20 and all(c in "ACGT" for c in s.upper())]
print(f"  → {len(filtered)} sequences pass validation")

# Subsample to 1000
if len(filtered) > 1000:
    filtered = random.sample(filtered, 1000)
    print(f"  → Subsampled to 1000 sequences (seed=42)")

# Write output
fasta_out.parent.mkdir(parents=True, exist_ok=True)
with open(fasta_out, 'w') as f:
    for header, seq in filtered:
        f.write(f">{header}\n{seq}\n")

# Compute SHA256
h = hashlib.sha256()
with open(fasta_out, 'rb') as f:
    for chunk in iter(lambda: f.read(8192), b''):
        h.update(chunk)
sha256 = h.hexdigest()

print(f"  ✓ Wrote {len(filtered)} sequences to {fasta_out}")
print(f"  → SHA256: {sha256}")
PYTHON_CURATE
else
    echo "  ✗ No parsed FASTA found. Check download."
    exit 1
fi

echo ""

# ==============================================================================
# DATASET 2: seed/sample - Minimal test dataset
# ==============================================================================
echo "[Step 3] Processing seed/sample dataset..."
echo "  Source: Derived from human/brunello (first 3 sequences)"
echo ""

# Create minimal sample for CI/tests
python3 << 'PYTHON_SAMPLE'
from pathlib import Path

source = Path("data/human/brunello/sequences.fasta")
dest = Path("data/seed/sample/sequences.fasta")

if not source.exists():
    print("  ✗ Source file not found. Run Brunello curation first.")
    exit(1)

# Read first 3 sequences
records = []
header = None
seq_lines = []
with open(source, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines)))
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        records.append((header, "".join(seq_lines)))

# Take first 3
sample = records[:3]

dest.parent.mkdir(parents=True, exist_ok=True)
with open(dest, 'w') as f:
    for header, seq in sample:
        f.write(f">{header}\n{seq}\n")

print(f"  ✓ Created sample with {len(sample)} sequences at {dest}")
PYTHON_SAMPLE

echo ""

# ==============================================================================
# VALIDATION
# ==============================================================================
echo "[Step 4] Validating all datasets..."
echo ""

# Validate using the validate_dataset.py script if available
if [ -f "scripts/validate_dataset.py" ]; then
    echo "  → Validating human/brunello..."
    python3 scripts/validate_dataset.py --dataset human/brunello || {
        echo "  ⚠ Validation failed for human/brunello"
    }
    
    echo ""
    echo "  → Validating seed/sample..."
    python3 scripts/validate_dataset.py --dataset seed/sample || {
        echo "  ⚠ Validation failed for seed/sample"  
    }
else
    # Basic validation if script not available
    echo "  → Basic validation (validate_dataset.py not found)..."
    for dataset in "data/human/brunello/sequences.fasta" "data/seed/sample/sequences.fasta"; do
        if [ -f "$dataset" ]; then
            count=$(grep -c "^>" "$dataset" || echo 0)
            size=$(wc -c < "$dataset" | tr -d ' ')
            echo "  ✓ $dataset: $count sequences, $size bytes"
        else
            echo "  ✗ $dataset: NOT FOUND"
        fi
    done
fi

echo ""
echo "=============================================="
echo "Data download and curation complete!"
echo "=============================================="
echo ""
echo "Curated datasets ready in:"
echo "  - data/human/brunello/sequences.fasta (1000 sgRNAs)"
echo "  - data/seed/sample/sequences.fasta (3 sgRNAs for testing)"
echo ""
echo "Raw files stored in data/raw/ (gitignored)"
echo ""
echo "To validate datasets manually:"
echo "  python scripts/validate_dataset.py --dataset human/brunello"
echo "  python scripts/validate_dataset.py --dataset seed/sample"
echo ""
