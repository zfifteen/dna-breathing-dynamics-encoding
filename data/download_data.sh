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
#   4. Curate and subsample to <=1000 sequences per dataset (using scripts/curate_and_subsample.py)
#   5. Validate the curated datasets
#
# Note: data/raw/ is gitignored. Only curated files in data/<dataset>/ are committed.

set -euo pipefail

# Ensure BASH_SOURCE[0] is set
if [[ -z "${BASH_SOURCE[0]:-}" ]]; then
  echo "Error: BASH_SOURCE[0] is not set. Please run this script with Bash, not sh or another shell." >&2
  exit 1
fi

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
import csv
import sys
from pathlib import Path

raw_path = Path("data/raw/broadgpp-brunello-library-contents.txt")
fasta_path = Path("data/raw/brunello_parsed.fasta")

if not raw_path.exists():
    print("  ✗ Raw file not found")
    sys.exit(1)

try:
    valid_seqs = []
    with open(raw_path, "r") as f:
        # The file uses carriage returns (\r) as line separators
        content = f.read().replace('\r\n', '\n').replace('\r', '\n')
    
    # Parse as CSV with tab delimiter
    lines = content.strip().split('\n')
    reader = csv.reader(lines, delimiter='\t')
    header = next(reader)
    
    # Find column indices
    # Expected columns: Target Gene ID, Target Gene Symbol, Target Transcript, 
    # Genomic Sequence, Position, Strand, sgRNA Target Sequence, Target Context Sequence, ...
    try:
        seq_idx = header.index("sgRNA Target Sequence")
        gene_idx = header.index("Target Gene Symbol")
        id_idx = header.index("Target Gene ID")
    except ValueError as e:
        # Fallback indices if header doesn't match exactly
        print(f"  ⚠ Warning: Column header not found ({e}), using fallback indices")
        seq_idx = 6   # sgRNA Target Sequence
        gene_idx = 1  # Target Gene Symbol
        id_idx = 0    # Target Gene ID
        # Validate fallback indices are reasonable
        if len(header) <= seq_idx:
            print(f"  ✗ Error: TSV has only {len(header)} columns, expected at least {seq_idx+1}")
            sys.exit(1)
    
    for row_num, row in enumerate(reader, 1):
        if len(row) > seq_idx:
            seq = row[seq_idx].strip().upper()
            # Validate sequence: 20bp and only ACGT (matching curate_and_subsample.py)
            if len(seq) == 20 and all(c in "ACGT" for c in seq):
                gene = row[gene_idx].strip() if len(row) > gene_idx else ""
                sg_id = row[id_idx].strip() if len(row) > id_idx else str(row_num)
                # Handle empty gene names
                if not gene:
                    gene = "Unknown"
                # Clean gene name for FASTA header
                gene = gene.replace(" ", "_")
                fasta_header = f">{gene}|{sg_id}"
                valid_seqs.append((fasta_header, seq))
    
    with open(fasta_path, "w") as f:
        for header, seq in valid_seqs:
            f.write(header + "\n" + seq + "\n")
    
    print(f"  ✓ Parsed {len(valid_seqs)} valid 20nt sequences to {fasta_path}")

except Exception as e:
    print(f"  ✗ Error during parsing: {e}")
    sys.exit(1)
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

# Curate and subsample to 1000 sequences using the canonical script
if [ -f "$BRUNELLO_FASTA" ]; then
    echo "  → Curating and subsampling to 1000 sequences..."
    echo "    (using scripts/curate_and_subsample.py)"
    python3 scripts/curate_and_subsample.py --dataset human/brunello --max-seqs 1000 --seed 42
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

# Create minimal sample for CI/tests by taking first 3 sequences from brunello
BRUNELLO_OUT="data/human/brunello/sequences.fasta"
SAMPLE_OUT="data/seed/sample/sequences.fasta"

if [ -f "$BRUNELLO_OUT" ]; then
    # Extract first 3 sequences (6 lines: 3 headers + 3 sequences)
    head -6 "$BRUNELLO_OUT" > "$SAMPLE_OUT"
    count=$(grep -c "^>" "$SAMPLE_OUT" || echo 0)
    echo "  ✓ Created sample with $count sequences at $SAMPLE_OUT"
else
    echo "  ✗ Source file not found. Run Brunello curation first."
    exit 1
fi

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
