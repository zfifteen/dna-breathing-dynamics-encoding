import csv
import os


def parse_local_tsv(fasta_file="data/processed/brunello.fasta"):
    # Parse local TSV, filter valid 20bp DNA seqs, output FASTA.
    tsv_file = "data/raw/Brunello_library_v1.1_sequences.txt"
    if not os.path.exists(tsv_file):
        print(f"TSV file not found: {tsv_file}")
        return

    valid_seqs = []
    with open(tsv_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)  # Skip header
        for row_num, row in enumerate(reader, 1):
            if len(row) >= 3:
                gene = row[0].strip()
                sg_id = row[1].strip()
                seq = row[2].strip().upper()
                if len(seq) == 20 and all(c in "ACGT" for c in seq):
                    header_str = f">{sg_id} | {gene}"
                    valid_seqs.append((header_str, seq))

    with open(fasta_file, "w") as f:
        for header_str, seq in valid_seqs:
            f.write(header_str + "\n" + seq + "\n")
    print(f"Generated {fasta_file} with {len(valid_seqs)} valid sgRNA sequences.")


if __name__ == "__main__":
    parse_local_tsv()
