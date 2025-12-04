# Data Plan

- [x] Create data directory structure with README.md, datasets.yml, .gitignore
- [x] Add scripts: download_raw.py, curate_and_subsample.py, validate_dataset.py
- [x] Implement Brunello dataset (human/brunello) with 1000 sequences
- [x] Add seed/sample dataset for CI/tests
- [x] Add DATA_PLAN.md documenting operational plan
- [x] Enhance download_raw.py with SHA256 hash computation
- [x] Enhance curate_and_subsample.py with strict 20nt ACGT filtering
- [x] Enhance validate_dataset.py with duplicate and sequence validation
- [x] Update datasets.yml with sample_size and filters specification
- [x] Add UTF-8 encoding to file operations
- [x] Add User-Agent header and explicit SSL verification to downloads
- [x] Extract FASTA_LINE_WIDTH constant for maintainability
- [x] Update DATA_PLAN.md with accurate execution steps and Quick Start guide
- [x] Create comprehensive download_data.sh that executes the full data plan
- [x] Fix parsing logic: Correctly extract Gene Symbol (column 1) instead of Context Sequence (column 7)
- [x] Unify execution protocol: Call scripts/curate_and_subsample.py instead of embedded HEREDOC
- [x] Add BASH_SOURCE validation for shell compatibility
- [x] Add proper error handling for TSV parsing
- [x] human/doench2016 dataset
- [x] human/depmap_subsample dataset
- [x] mouse/gecko_v2 dataset

## Quick Start
./data/download_data.sh
