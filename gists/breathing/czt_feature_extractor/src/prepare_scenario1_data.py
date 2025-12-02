import random
from typing import List, Tuple

def read_fasta(fasta_path: str) -> List[Tuple[str, str]]:
    seqs = []
    header = None
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
            else:
                if header:
                    seqs.append((header, line))
                    header = None
    return seqs

def prepare_gc_dataset(
    input_fasta: str = "data/processed/brunello.fasta",
    output_fasta: str = "data/processed/brunello_gc1000.fasta",
    n_per_group: int = 500,
    seed: int = 42,
):
    random.seed(seed)
    all_seqs = read_fasta(input_fasta)
    
    high_gc = []
    low_gc = []
    
    for header, seq in all_seqs:
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        # Use original cutoffs: High > 50%, Low < 40%
        if gc > 0.50:
            high_gc.append((header, seq))
        elif gc < 0.40:
            low_gc.append((header, seq))
            
    print(f"Found {len(high_gc)} high GC (>50%) and {len(low_gc)} low GC (<40%) sequences.")
    
    if len(high_gc) < n_per_group or len(low_gc) < n_per_group:
        print(f"Warning: Not enough sequences to satisfy n={n_per_group}. Using max available.")
        n_high = min(len(high_gc), n_per_group)
        n_low = min(len(low_gc), n_per_group)
    else:
        n_high = n_per_group
        n_low = n_per_group
        
    sampled_high = random.sample(high_gc, n_high)
    sampled_low = random.sample(low_gc, n_low)
    
    with open(output_fasta, "w") as f:
        for h, s in sampled_high:
            # Tag header with group for analysis script: >ID|Label
            # Replace existing pipes in ID to avoid confusion
            clean_h = h.replace(">", "").replace("|", "_")
            f.write(f">{clean_h}|high_gc\n{s}\n")
        for h, s in sampled_low:
            clean_h = h.replace(">", "").replace("|", "_")
            f.write(f">{clean_h}|low_gc\n{s}\n")
            
    print(f"Created {output_fasta} with {n_high} high GC and {n_low} low GC sequences.")

if __name__ == "__main__":
    prepare_gc_dataset()
