id: human_brunello
name: Brunello sgRNA library (subset)
species: Homo sapiens
taxon_id: 9606
source_url: https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt
source_doi: 10.1038/nbt.3437
source_date: "2025-12-03"
raw_filename: brunello_parsed.fasta
raw_source_file: broadgpp-brunello-library-contents.txt
committed_files:
  - sequences.fasta
committed_sha256: 4874d56903c09bbd1290d643ebbfd8a0c70d228bf64c9e0ff5c457aa743ede95
citation: "Doench et al., 2016, Nature Biotechnology. DOI: 10.1038/nbt.3437"
license: "see source (Addgene / paper supplement)"
curation_command: |
  ./data/download_data.sh
  # Or manually:
  # python scripts/curate_and_subsample.py --dataset human/brunello --max-seqs 1000 --seed 42
notes: |
  Raw data: broadgpp-brunello-library-contents.txt (~9MB, 77,441 sgRNAs)
  Downloaded from Addgene Brunello library page or copied from repo fallback.
  Parsed TSV to FASTA, filtered to 20nt ACGT sequences, subsampled 1000 (seed=42).
  Curated subset contains 1000 representative sgRNAs for testing and validation.
