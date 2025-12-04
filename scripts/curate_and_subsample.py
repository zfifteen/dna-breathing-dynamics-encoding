#!/usr/bin/env python3
import yaml
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import numpy as np
from pathlib import Path
import os

random.seed(42)
np.random.seed(42)


def load_datasets():
    with open("data/datasets.yml", "r") as f:
        return yaml.safe_load(f)["datasets"]


def curate_and_subsample(dataset):
    name = dataset["name"]
    sample_size = dataset["sample_size"]
    raw_dir = Path("data/raw")
    processed_dir = Path(f"data/{name}")
    processed_dir.mkdir(parents=True, exist_ok=True)

    # Assume TSV with columns like Gene, Sequence, Score
    filename = [
        f for f in os.listdir(raw_dir) if name.replace("/", "_") in f or "txt" in f
    ][0]  # Adjust based on filename
    filepath = raw_dir / filename
    df = pd.read_csv(filepath, sep="\t")

    # Filter 20nt ACGT
    df = df[df["Sequence"].str.len() == 20]
    df = df[df["Sequence"].str.match("^[ACGT]+$")]

    # Subsample
    if len(df) > sample_size:
        df = df.sample(n=sample_size)

    # Write FASTA
    fasta_path = processed_dir / "sequences.fasta"
    with open(fasta_path, "w") as f:
        for i, row in df.iterrows():
            seq = Seq(row["Sequence"])
            record = SeqRecord(seq, id=f"{row.get('Gene', i)}", description="")
            SeqIO.write(record, f, "fasta")

    print(f"Curated and subsampled {name} to {len(df)} sequences")


if __name__ == "__main__":
    datasets = load_datasets()
    for ds in datasets:
        curate_and_subsample(ds)
