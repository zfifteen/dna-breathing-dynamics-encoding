#!/usr/bin/env python3
"""
Small helper to download raw dataset archives into data/raw/.

Usage: python scripts/download_raw.py --url <url> --out <filename>

Note: This script does not bypass site terms or authentication. If a resource
requires manual acceptance (e.g., DepMap), perform that in a browser and place
the file into data/raw/ manually.
"""

import argparse
import hashlib
import sys
from pathlib import Path

import requests


def sha256_of_file(path: Path) -> str:
    """Compute SHA256 hash of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="Download raw dataset archives into data/raw/."
    )
    parser.add_argument("--url", required=True, help="URL to download from")
    parser.add_argument("--out", required=True, help="filename under data/raw/")
    parser.add_argument(
        "--dataset",
        help=(
            "Dataset identifier (e.g., brunello). "
            "If provided, creates data/raw/<dataset>/ subdirectory."
        ),
    )
    args = parser.parse_args()

    if args.dataset:
        out_path = Path("data") / "raw" / args.dataset / args.out
    else:
        out_path = Path("data") / "raw" / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {args.url} -> {out_path}")
    headers = {"User-Agent": "dna-breathing-dynamics-encoding/data-downloader"}
    r = requests.get(
        args.url,
        stream=True,
        timeout=60,
        verify=True,
        headers=headers,
    )
    if r.status_code != 200:
        print(f"Failed to download: HTTP {r.status_code}")
        sys.exit(2)

    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)

    # Compute and display SHA256
    file_hash = sha256_of_file(out_path)
    file_size = out_path.stat().st_size

    print(f"Downloaded to {out_path}")
    print(f"  Size: {file_size} bytes")
    print(f"  SHA256: {file_hash}")
    print()
    print("Add this to your METADATA.md:")
    print(f"  raw_sha256: {file_hash}")


if __name__ == "__main__":
    main()
