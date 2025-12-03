#!/usr/bin/env python3
"""
Curate and subsample raw datasets into small curated artifacts.

Usage: python scripts/curate_and_subsample.py --dataset <dataset_id>

This script expects a `data/<dataset>/METADATA.md` file describing the raw filename
and `data/raw/<raw_filename>` to exist. It will write curated outputs into
`data/<dataset>/` as `sequences.fasta` and `labels.csv` (if labels are present).

The script enforces these constraints:
- At most 1000 sequences (configurable via --max-seqs)
- Sequences must be exactly 20nt (configurable via --seq-length)
- Only A/C/G/T bases allowed (strict mode)
- Output files must be <= 5 MB (checked after writing)

The script is intentionally small and avoids external heavy deps; it uses
Biopython if available for robust FASTA parsing, otherwise falls back to a
simple parser.

"""

import argparse
import hashlib
import random
from pathlib import Path

# Allowed bases for strict mode (ACGT only as per data plan)
ALLOWED_STRICT = set("ACGTacgt")
# Allowed bases for relaxed mode (includes IUPAC ambiguity codes)
ALLOWED_RELAXED = set(list("ATGCatgcNRYKMSWBDHVrykmswbdhv"))

DEFAULT_MAX_SEQS = 1000
DEFAULT_MAX_BYTES = 5 * 1024 * 1024  # 5 MB
DEFAULT_SEQ_LENGTH = 20


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
                f.write(seq[i:i + 80] + "\n")


def filter_records(records, allowed_bases, seq_length=None):
    """Filter records by allowed bases and optionally exact sequence length."""
    good = []
    for header, seq in records:
        if not seq:
            continue
        if any(ch not in allowed_bases for ch in seq):
            continue
        if seq_length is not None and len(seq) != seq_length:
            continue
        good.append((header, seq.upper()))
    return good


def main():
    parser = argparse.ArgumentParser(
        description="Curate and subsample raw datasets into small curated artifacts."
    )
    parser.add_argument("--dataset", required=True, help="Dataset identifier")
    parser.add_argument(
        "--max-seqs",
        type=int,
        default=DEFAULT_MAX_SEQS,
        help=f"Maximum number of sequences (default: {DEFAULT_MAX_SEQS})",
    )
    parser.add_argument(
        "--max-bytes",
        type=int,
        default=DEFAULT_MAX_BYTES,
        help=f"Maximum output file size in bytes (default: {DEFAULT_MAX_BYTES})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible subsampling (default: 42)",
    )
    parser.add_argument(
        "--seq-length",
        type=int,
        default=DEFAULT_SEQ_LENGTH,
        help=(
            f"Required sequence length (default: {DEFAULT_SEQ_LENGTH}). "
            "Set to 0 to disable length filtering."
        ),
    )
    parser.add_argument(
        "--relaxed",
        action="store_true",
        help="Allow IUPAC ambiguity codes (N, R, Y, etc.) in addition to ACGT",
    )
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
        print(f"Raw file not found at {raw_path}. Place raw file in data/raw/")
        raise SystemExit(2)

    # For now assume raw file is a plain FASTA; future: support zip/tar extracts
    records = read_fasta_simple(raw_path)
    print(f"Read {len(records)} records from {raw_path}")

    # Filter by bases and length
    allowed_bases = ALLOWED_RELAXED if args.relaxed else ALLOWED_STRICT
    seq_length = args.seq_length if args.seq_length > 0 else None
    records = filter_records(records, allowed_bases, seq_length)
    print(f"{len(records)} records remain after filtering")
    if seq_length:
        bases_msg = "ACGT+IUPAC" if args.relaxed else "ACGT only"
        print(f"  (required length: {seq_length}nt, bases: {bases_msg})")

    if len(records) == 0:
        print("No valid records after filtering; aborting")
        raise SystemExit(2)

    # subsample
    if len(records) > args.max_seqs:
        records = random.sample(records, args.max_seqs)
        print(f"Subsampled to {len(records)} sequences (seed={args.seed})")

    out_dir = base
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fasta = out_dir / "sequences.fasta"
    write_fasta(records, out_fasta)

    size = out_fasta.stat().st_size
    if size > args.max_bytes:
        print(f"Output file size {size} exceeds max bytes {args.max_bytes}; aborting")
        raise SystemExit(2)

    # Compute SHA256 for reproducibility
    file_hash = sha256_of_file(out_fasta)

    print(f"Wrote curated FASTA to {out_fasta}")
    print(f"  Sequences: {len(records)}")
    print(f"  Size: {size} bytes")
    print(f"  SHA256: {file_hash}")
    print()
    print("Update METADATA.md with:")
    print(f"  committed_sha256: {file_hash}")
    print()
    print("Reproducibility command:")
    cmd_parts = [
        f"python scripts/curate_and_subsample.py --dataset {args.dataset}",
        f"--max-seqs {args.max_seqs}",
        f"--seed {args.seed}",
    ]
    if args.seq_length != DEFAULT_SEQ_LENGTH:
        cmd_parts.append(f"--seq-length {args.seq_length}")
    if args.relaxed:
        cmd_parts.append("--relaxed")
    print("  " + " \\\n    ".join(cmd_parts))
    print("Done")


if __name__ == "__main__":
    main()
