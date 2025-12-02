import requests
import zipfile
import io
import csv
from typing import List


def download_and_extract(
    url: str = "https://media.addgene.org/files/Brunello_library_v1.1_sequences.txt.zip",
    extract_to: str = "temp_brunello.txt",
) -> str:
    """Download ZIP and extract the TXT file."""
    response = requests.get(url)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        z.extractall()
    # Assume the extracted file is Brunello_library_v1.1_sequences.txt
    return "Brunello_library_v1.1_sequences.txt"


def parse_tsv_to_fasta(
    tsv_file: str, fasta_file: str = "../data/processed/brunello.fasta"
) -> None:
    """Parse TSV, filter valid 20bp DNA seqs, output FASTA."""
    valid_seqs = []
    with open(tsv_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) >= 3:  # At least Gene, ID, Sequence
                seq = (
                    row[2].strip().upper()
                )  # Assuming sequence is in column 2 (0-based)
                if len(seq) == 20 and all(c in "ACGT" for c in seq):
                    header = (
                        f">{row[0]}|{row[1]}|{row[2][:10]}..."  # Gene|ID|partial seq
                    )
                    valid_seqs.append((header, seq))

    with open(fasta_file, "w") as f:
        for header, seq in valid_seqs:
            f.write(header + "\n" + seq + "\n")
    print(f"Output FASTA with {len(valid_seqs)} valid sequences to {fasta_file}")


if __name__ == "__main__":
    tsv = download_and_extract()
    parse_tsv_to_fasta(tsv)
