# data/

This directory stores curated, small subsets of real biological sequence datasets used by this project.

Policy
- No Git LFS. Raw archives go into `data/raw/` and are gitignored.
- Max committed file size: 5 MB per file.
- Curated artifacts will contain at most 1000 sequences per dataset (unless the dataset legitimately contains fewer after filtering).
- Single repository license: MIT.

Layout
- data/
  - datasets.yml           # machine-readable catalog
  - INDEX.md               # human-readable index (created by you)
  - DATA_PLAN.md           # operational plan for the data directory
  - SOURCES.md             # canonical sources (created previously)
  - raw/                   # raw downloaded archives (gitignored)
  - human/                 # human datasets (brunello, etc.)
  - mouse/                 # mouse datasets
  - seed/sample/           # tiny sample dataset for CI and examples

How to add a new dataset
1. Add an entry to `data/datasets.yml` with metadata fields.
2. Add a directory `data/<dataset_id>/` and create a `METADATA.md` file with source URL, citation, and license.
3. Place raw artifact in `data/raw/` (do not commit). Compute sha256 and write to `METADATA.md`.
4. Run `scripts/curate_and_subsample.py --dataset <dataset_id>` to produce curated artifacts.
5. Run `scripts/validate_dataset.py --dataset <dataset_id>` to validate.
6. Commit curated artifacts and `METADATA.md` (do not commit `data/raw/`).

Reproducibility
- All curator operations should be recorded in the dataset `METADATA.md` with exact command lines, seed used for sampling, and raw file checksum.

Contact
- If you have questions about dataset licensing or provenance, open an issue or consult `SOURCES.md`.

