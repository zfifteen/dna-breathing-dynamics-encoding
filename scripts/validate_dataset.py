#!/usr/bin/env python3
from Bio import SeqIO
from pathlib import Path
import os


def validate_dataset(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    lengths = [len(r.seq) for r in records]
    bases = all(all(b in "ACGT" for b in r.seq) for r in records)
    duplicates = len(records) != len(set(r.id for r in records))
    print(
        f"Validation for {fasta_path}: {len(records)} sequences, all 20nt: {all(l == 20 for l in lengths)}, valid bases: {bases}, no duplicates: {not duplicates}"
    )


if __name__ == "__main__":
    for fasta in Path("data").rglob("*.fasta"):
        validate_dataset(fasta)
