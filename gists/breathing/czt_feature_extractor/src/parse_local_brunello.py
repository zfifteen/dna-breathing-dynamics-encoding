import csv
import sys

def parse_brunello_txt(
    txt_path: str = "data/raw/broadgpp-brunello-library-contents.txt",
    fasta_file: str = "data/processed/brunello.fasta",
):
    valid_seqs = []
    try:
        with open(txt_path, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            
            # Find the column index for "sgRNA Target Sequence"
            try:
                seq_idx = header.index("sgRNA Target Sequence")
                gene_idx = header.index("Target Gene Symbol")
                id_idx = header.index("Target Gene ID")
            except ValueError:
                # Fallback indices if header doesn't match exactly
                seq_idx = 6
                gene_idx = 1
                id_idx = 0
                
            for row_num, row in enumerate(reader, 1):
                if len(row) > seq_idx:
                    seq = row[seq_idx].strip().upper()
                    # Validate sequence: 20bp and only ATGC
                    if len(seq) == 20 and all(c in "ATGC" for c in seq):
                        gene = row[gene_idx] if len(row) > gene_idx else "Unknown"
                        sg_id = row[id_idx] if len(row) > id_idx else str(row_num)
                        # Clean gene name for FASTA header
                        gene = gene.replace(" ", "_")
                        fasta_header = f">{gene}|{sg_id}"
                        valid_seqs.append((fasta_header, seq))
                        
        with open(fasta_file, "w") as f:
            for header, seq in valid_seqs:
                f.write(header + "\n" + seq + "\n")
                
        print(f"Output FASTA with {len(valid_seqs)} valid sequences to {fasta_file}")
        
    except FileNotFoundError:
        print(f"Error: File not found at {txt_path}")
        sys.exit(1)

if __name__ == "__main__":
    parse_brunello_txt()