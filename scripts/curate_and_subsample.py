#!/usr/bin/env python3
"""
Curate and subsample raw datasets into small curated artifacts.

Usage: python scripts/curate_and_subsample.py --dataset <dataset_id>

This script expects a `data/<dataset>/METADATA.md` file describing the raw filename
and `data/raw/<raw_filename>` to exist. It will write curated outputs into
`data/<dataset>/` as `sequences.fasta` and `labels.csv` (if labels are present).

The script enforces these constraints:
- At most 1000 sequences (configurable via --max-seqs)
- Only ATGC and common IUPAC codes allowed (basic filter)
- Output files must be <= 5 MB (checked after writing)

The script is intentionally small and avoids external heavy deps; it uses
Biopython if available for robust FASTA parsing, otherwise falls back to a
simple parser.

"""

import argparse
import os
import random
import textwrap
import hashlib
from pathlib import Path

# Allowed bases (IUPAC basic set)
ALLOWED = set(list("ATGCatgcNRYKMSWBDHVrykmswbdhv"))

DEFAULT_MAX_SEQS = 1000
DEFAULT_MAX_BYTES = 5 * 1024 * 1024  # 5 MB


def sha256_of_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fasta_simple(path):
    """Simple FASTA reader returning list of (header, seq)"""
    records = []
    header = None
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records


def write_fasta(records, path):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            # wrap at 80
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def filter_records(records):
    good = []
    for header, seq in records:
        if not seq:
            continue
        if any(ch not in ALLOWED for ch in seq):
            continue
        good.append((header, seq.upper()))
    return good


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--max-seqs", type=int, default=DEFAULT_MAX_SEQS)
    parser.add_argument("--max-bytes", type=int, default=DEFAULT_MAX_BYTES)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)

    base = Path("data") / args.dataset
    meta = base / "METADATA.md"
    if not meta.exists():
        print(f"METADATA.md not found for dataset {args.dataset} at {meta}")
        raise SystemExit(2)

    # parse raw filename from METADATA.md (very small parser)
    raw_filename = None
    with open(meta, "r") as f:
        for line in f:
            if line.lower().startswith("raw_filename:"):
                raw_filename = line.split(":", 1)[1].strip()
                break
    if not raw_filename:
        print("raw_filename not specified in METADATA.md")
        raise SystemExit(2)

    raw_path = Path("data") / "raw" / raw_filename
    if not raw_path.exists():
        print(f"Raw file not found at {raw_path}. Place the raw archive in data/raw/")
        raise SystemExit(2)

    # For now assume raw file is a plain FASTA; future: support zip/tar extracts
    records = read_fasta_simple(raw_path)
    print(f"Read {len(records)} records from {raw_path}")

    records = filter_records(records)
    print(f"{len(records)} records remain after filtering by allowed bases")

    if len(records) == 0:
        print("No valid records after filtering; aborting")
        raise SystemExit(2)

    # subsample
    if len(records) > args.max_seqs:
        records = random.sample(records, args.max_seqs)
        print(f"Subsampled to {len(records)} sequences")

    out_dir = base
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fasta = out_dir / "sequences.fasta"
    write_fasta(records, out_fasta)

    # optional labels.csv generation: try to parse header if it contains a label
    # This is simplified; real scripts should consult dataset-specific parsers.

    size = out_fasta.stat().st_size
    if size > args.max_bytes:
        print(f"Output file size {size} exceeds max bytes {args.max_bytes}; aborting")
        raise SystemExit(2)

    print(f"Wrote curated FASTA to {out_fasta} ({size} bytes)")
    print("Done")


if __name__ == "__main__":
    main()

