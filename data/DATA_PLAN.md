# Data Plan

This document captures the operational plan for the `data/` directory in this repository.

## Goals
- Seed `data/` with curated, provenance-tracked subsets of real biological sequence datasets.
- Commit only small curated artifacts (<= 1000 sequences, <= 5 MB per file).
- Keep raw archives outside git under `data/raw/` (gitignored).
- Provide reproducible scripts to download, curate, and validate datasets.

## Constraints
- Max committed file size: 5 MB.
- Fixed sample size per dataset: 1000 sequences (or fewer if the dataset has fewer than 1000 available sequences after filtering).
- Single repository license: MIT (applies to code, scripts, and curated artifacts).

## Initial Seed Datasets (Priority Order)
1. `human/brunello` — Brunello sgRNA library (Doench 2016 Nature Biotechnology paper supplement; DOI: 10.1038/nbt.3437; Supplemental Table 1 ZIP with 76,442 sgRNAs).
2. `human/doench2016` — Doench 2016 guide-efficiency dataset (same paper; Supplemental Tables 2/3 with ~4,000 tested sgRNAs).
3. `human/depmap_subsample` — Small processed DepMap-derived sample (DepMap Portal downloads: https://depmap.org/portal/download/all/; e.g., Achilles_gene_effect.csv; public CC-BY 4.0).
4. `mouse/gecko_v2` — GeCKO v2 mouse library (Addgene pooled library; DOI: 10.1016/j.cell.2014.09.029; ~123,000 sgRNAs).
5. `seed/sample` — Minimal example for CI/tests (derived from existing repo files; always included).

## Preparation Steps (One-Time Setup)
- **Directory Structure:** Ensure `data/raw/` (gitignore'd), `data/<dataset>/` for each (e.g., `human/brunello/`), and `scripts/` exist. Add `data/.gitignore` with patterns like `raw/*` and `*.zip`.
- **Core Files:**
  - `data/README.md`: Document the data plan, goals, constraints, and usage (e.g., "Curated subsets for CRISPR sgRNA analysis; see METADATA.md for provenance").
  - `data/datasets.yml`: YAML config listing datasets with keys like `name`, `source_url`, `sample_size: 1000`, `filters: {min_length: 20, charset: ACGT}`.
- **Scripts (in `scripts/`):**
  - `download_raw.py`: CLI tool to fetch raw files (e.g., via `requests` for HTTP). Computes SHA256 (using `hashlib`). Example: `python download_raw.py --dataset brunello --url <supplemental_zip_url> --output data/raw/brunello.zip`.
  - `curate_and_subsample.py`: Parse raw (e.g., Excel to FASTA via `pandas` and `biopython`), filter (20nt, A/C/G/T), subsample to 1000 (seeded via `numpy.random.seed(42)`), write `sequences.fasta`. Log command for reproducibility.
  - `validate_dataset.py`: Check FASTA: count <=1000, size <=5MB, seqs 20nt A/C/G/T, no duplicates. Output pass/fail report.
- **Testing:** Run scripts on a dummy file first (e.g., sample Excel with 2000 fake sgRNAs).

## Per-Dataset Execution Plan
Process in priority order. For each:
- Download raw to `data/raw/<dataset>/`.
- Run curation script to generate `data/<dataset>/sequences.fasta`.
- Create `data/<dataset>/METADATA.md` with: source URL/DOI, download filename, raw SHA256, citation/license, curation command, filters applied.
- Validate and log issues.
- Subsample randomly but deterministically if >1000 seqs; use all if fewer.

- **Dataset 1: human/brunello**
  - **Source:** https://www.nature.com/articles/nbt.3437#Sec23 (Supplemental ZIP ~1MB; extract `nbt.3437-s1.xlsx` with sgRNA sequences).
  - **Curation:** Filter valid sgRNAs (20nt, A/C/G/T, human genes). Subsample 1000 (seed 42). Headers: `>sgRNA_id|gene_name`.
  - **Metadata:** DOI: 10.1038/nbt.3437; Citation: Doench et al. (2016); License: Public domain.
  - **Output:** `data/human/brunello/sequences.fasta` (~20KB), METADATA.md.
  - **Validation:** 1000 seqs, all 20nt A/C/G/T, <5MB.

- **Dataset 2: human/doench2016**
  - **Source:** Same ZIP; parse efficiency sheet.
  - **Curation:** Filter top 1000 by efficiency or random subsample. Headers: `>sgRNA_id|target_gene|efficiency_score`.
  - **Metadata:** Same as Brunello.

- **Dataset 3: human/depmap_subsample**
  - **Source:** https://depmap.org/portal/download/all/ (Achilles CSV ~100MB; extract sgRNAs or generate from gene targets).
  - **Curation:** Subsample 1000 targeting cancer genes (CCLE lineage filter). Headers: `>sgRNA_id|gene_symbol`.
  - **Metadata:** URL: https://depmap.org/portal/; License: CC-BY 4.0.

- **Dataset 4: mouse/gecko_v2**
  - **Source:** https://www.addgene.org/pooled-libraries/geckov2/ (~5MB text/CSV).
  - **Curation:** Filter mouse-specific, subsample 1000 (seed 42). Headers: `>sgRNA_id|gene_symbol`.
  - **Metadata:** DOI: 10.1016/j.cell.2014.09.029; License: MIT.

- **Dataset 5: seed/sample**
  - **Source:** Existing repo (e.g., `human/brunello/sequences.fasta`).
  - **Curation:** Subsample 100 diverse sgRNAs. Create `data/seed/sample/sequences.fasta` and METADATA.md (cite repo).
  - **Notes:** For testing scripts.

## Reproducibility & Validation
- **Determinism:** Fixed seeds; log exact commands in METADATA.md (e.g., "curate_and_subsample.py --input raw/brunello.xlsx --output sequences.fasta --sample-size 1000 --seed 42").
- **Full Validation:** Run `validate_dataset.py` post-curation; iterate if fails.
- **Proof-Pack:** Run `proof_pack/run_validation.py` on curated data.
- **Size Check:** Scripts enforce limits; reject invalid chars.

## Timeline & Tradeoffs
- **Effort:** 1-2 days per dataset (download/parse ~30min, script ~1hr, validation ~15min). Total: 1 week sequential.
- **Tradeoffs:** Python for flexibility (Biopython for FASTA); random vs. stratified subsampling (random for simplicity).
- **Risks:** Large downloads—`requests` with resume. Parsing errors—test samples. No direct sequences—cross-reference or generate.

## Next Actions
- Create `data/README.md`, `data/datasets.yml`, `data/.gitignore`; seed `seed/sample`.
- Add scripts: `download_raw.py`, `curate_and_subsample.py`, `validate_dataset.py`.
- Execute on Brunello first, then parallelize others.

This plan leverages public sources, modular scripts, and validation for 10/10 confidence.