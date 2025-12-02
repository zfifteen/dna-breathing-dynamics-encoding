#!/usr/bin/env python3
import csv
import re
from collections import Counter

tsv_path = "data/raw/Brunello_library_v1.1_sequences.txt"
fasta_path = "data/processed/depmap_subsample.fasta"

with open(tsv_path, "r") as f:
    reader = csv.reader(f, delimiter="\t")
    sequences = []
    for row in reader:
        if len(row) < 2:
            continue
        seq = row[0].strip().upper()
        gene = row[1].strip()
        if len(seq) != 20:
            continue
        if not re.match(r"^[ATGC]+$", seq):
            continue
        gc_count = seq.count("G") + seq.count("C")
        gc_percent = gc_count / len(seq)
        if gc_percent > 0.5 or gc_percent < 0.4:
            group = "high_gc" if gc_percent > 0.5 else "low_gc"
            header = f">{group}_{gene}"
            sequences.append((header, seq))

# Write FASTA
with open(fasta_path, "w") as f:
    for header, seq in sequences:
        f.write(header + "\n")
        f.write(seq + "\n")

print(f"Subsampled {len(sequences)} sequences to {fasta_path}")
