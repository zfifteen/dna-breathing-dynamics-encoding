#!/usr/bin/env python3
"""
Validate curated datasets for size, sequence count, and provenance.

Usage: python scripts/validate_dataset.py --dataset <dataset_id>

Checks performed:
- data/<dataset>/METADATA.md exists and contains required fields
- Committed FASTA exists and has <= 1000 sequences
- Committed FASTA size <= 5 MB
- Committed files listed in METADATA.md exist and their sha256 match (if provided)

Returns non-zero exit on validation failures.
"""

import argparse
from pathlib import Path
import sys
import re
import hashlib


REQUIRED_META_FIELDS = ["raw_filename", "source_url", "license"]
MAX_SEQS = 1000
MAX_BYTES = 5 * 1024 * 1024


def sha256_of_file(path: Path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_metadata(meta_path: Path):
    meta = {}
    with open(meta_path, "r") as f:
        for line in f:
            if ":" in line:
                k, v = line.split(":", 1)
                meta[k.strip()] = v.strip()
    return meta


def count_fasta_records(path: Path):
    c = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                c += 1
    return c


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    args = parser.parse_args()

    base = Path("data") / args.dataset
    meta = base / "METADATA.md"
    if not meta.exists():
        print(f"METADATA.md missing at {meta}")
        sys.exit(2)

    meta_kv = parse_metadata(meta)
    for field in REQUIRED_META_FIELDS:
        if field not in meta_kv:
            print(f"Required metadata field '{field}' missing in {meta}")
            sys.exit(2)

    fasta = base / "sequences.fasta"
    if not fasta.exists():
        print(f"Curated fasta missing: {fasta}")
        sys.exit(2)

    recs = count_fasta_records(fasta)
    if recs > MAX_SEQS:
        print(f"Too many sequences ({recs}) in {fasta}; max allowed is {MAX_SEQS}")
        sys.exit(2)

    size = fasta.stat().st_size
    if size > MAX_BYTES:
        print(f"File too large ({size} bytes) for {fasta}; max is {MAX_BYTES} bytes")
        sys.exit(2)

    # optional: check committed_sha256 in metdata
    if "committed_sha256" in meta_kv:
        expected = meta_kv["committed_sha256"]
        actual = sha256_of_file(fasta)
        if actual != expected:
            print(f"SHA256 mismatch for {fasta}: expected {expected} got {actual}")
            sys.exit(2)

    print("Validation passed")
    sys.exit(0)


if __name__ == "__main__":
    main()

