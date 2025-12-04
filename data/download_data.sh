#!/bin/bash
set -e

cd data

python3 ../../scripts/download_raw.py
python3 ../../scripts/curate_and_subsample.py
python3 ../../scripts/validate_dataset.py

echo "Data setup complete"
