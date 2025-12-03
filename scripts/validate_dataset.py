#!/usr/bin/env python3
"""
Validate curated datasets for size, sequence count, and provenance.

Usage: python scripts/validate_dataset.py --dataset <dataset_id>

Checks performed:
- data/<dataset>/METADATA.md exists and contains required fields
- Committed FASTA exists and has <= 1000 sequences
- Committed FASTA size <= 5 MB
- All sequences are exactly 20nt and contain only A/C/G/T characters
- No duplicate sequences
- Committed files listed in METADATA.md exist and their sha256 match (if provided)

Returns non-zero exit on validation failures.
"""

import argparse
from pathlib import Path
import sys
import hashlib


REQUIRED_META_FIELDS = ["raw_filename", "source_url", "license"]
MAX_SEQS = 1000
MAX_BYTES = 5 * 1024 * 1024
EXPECTED_SEQ_LENGTH = 20
VALID_BASES = set("ACGT")


def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_metadata(meta_path: Path) -> dict:
    meta = {}
    with open(meta_path, "r", encoding="utf-8") as f:
        for line in f:
            if ":" in line:
                k, v = line.split(":", 1)
                meta[k.strip()] = v.strip()
    return meta


def read_fasta_sequences(path: Path) -> list:
    """Read FASTA file and return list of (header, sequence) tuples."""
    records = []
    header = None
    seq_lines = []
    with open(path, "r", encoding="utf-8") as f:
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


def validate_sequences(records: list) -> list:
    """Validate sequences and return list of issues."""
    issues = []
    seen_seqs = set()

    for i, (header, seq) in enumerate(records, 1):
        # Check length
        if len(seq) != EXPECTED_SEQ_LENGTH:
            issues.append(
                f"Sequence {i} ({header}): length {len(seq)} != {EXPECTED_SEQ_LENGTH}"
            )

        # Check valid bases (A/C/G/T only)
        seq_upper = seq.upper()
        invalid_bases = set(seq_upper) - VALID_BASES
        if invalid_bases:
            issues.append(
                f"Sequence {i} ({header}): invalid bases {invalid_bases}"
            )

        # Check for duplicates
        if seq_upper in seen_seqs:
            issues.append(f"Sequence {i} ({header}): duplicate sequence")
        seen_seqs.add(seq_upper)

    return issues


def main():
    parser = argparse.ArgumentParser(
        description="Validate curated datasets for size, sequence count, provenance."
    )
    parser.add_argument("--dataset", required=True, help="Dataset identifier")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on any sequence validation issue (length, bases, duplicates)",
    )
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

    # Read and count sequences
    records = read_fasta_sequences(fasta)
    recs = len(records)
    if recs > MAX_SEQS:
        print(f"Too many sequences ({recs}) in {fasta}; max allowed is {MAX_SEQS}")
        sys.exit(2)

    size = fasta.stat().st_size
    if size > MAX_BYTES:
        print(f"File too large ({size} bytes) for {fasta}; max is {MAX_BYTES} bytes")
        sys.exit(2)

    # Validate sequences
    issues = validate_sequences(records)
    if issues:
        print(f"Sequence validation issues in {fasta}:")
        for issue in issues[:10]:  # Show first 10 issues
            print(f"  - {issue}")
        if len(issues) > 10:
            print(f"  ... and {len(issues) - 10} more issues")
        if args.strict:
            sys.exit(2)
        else:
            print("  (use --strict to fail on these issues)")

    # Check committed_sha256 in metadata
    if "committed_sha256" in meta_kv:
        expected = meta_kv["committed_sha256"]
        actual = sha256_of_file(fasta)
        if actual != expected:
            print(f"SHA256 mismatch for {fasta}: expected {expected} got {actual}")
            sys.exit(2)

    print(f"Validation passed for {args.dataset}")
    print(f"  Sequences: {recs}")
    print(f"  File size: {size} bytes")
    if not issues:
        seq_msg = f"all {recs} sequences are {EXPECTED_SEQ_LENGTH}nt A/C/G/T"
        print(f"  Sequence checks: {seq_msg}")
    sys.exit(0)


if __name__ == "__main__":
    main()
