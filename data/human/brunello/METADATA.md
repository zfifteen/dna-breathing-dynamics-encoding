id: human_brunello
name: Brunello sgRNA library (subset)
species: Homo sapiens
taxon_id: 9606
source_url: https://media.addgene.org/files/Brunello_library_v1.1_sequences.txt.zip
source_doi: 10.1038/nbt.3437
source_date: "2025-12-03"
raw_filename: brunello_parsed.fasta
raw_zip_filename: Brunello_library_v1.1_sequences.txt.zip
raw_sha256: "pending (official download blocked; used repo-processed FASTA)"
committed_files:
  - sequences.fasta
committed_sha256: 4874d56903c09bbd1290d643ebbfd8a0c70d228bf64c9e0ff5c457aa743ede95
citation: "Doench et al., 2016, Nature Biotechnology. DOI: 10.1038/nbt.3437"
license: "see source (Addgene / paper supplement)"
notes: |
  Official download from Addgene URL returned HTTP 403 (attempted with user-agent header on 2025-12-03).
  Used existing processed FASTA from repo (gists/.../data/processed/brunello.fasta) as raw input.
  Downloaded from Addgene supplemental file; parsed to include only 20bp ACGT sgRNA sequences.
  Curated subsample of up to 1000 sequences created with scripts/curate_and_subsample.py (seed=42).
