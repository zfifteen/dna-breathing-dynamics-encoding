#!/usr/bin/env python3
"""
Small helper to download raw dataset archives into data/raw/.

Usage: python scripts/download_raw.py --url <url> --out <filename>

Note: This script does not bypass site terms or authentication. If a resource
requires manual acceptance (e.g., DepMap), perform that in a browser and place
the file into data/raw/ manually.
"""

import argparse
import requests
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--out", required=True, help="filename under data/raw/")
    args = parser.parse_args()

    out_path = Path("data") / "raw" / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {args.url} -> {out_path}")
    r = requests.get(args.url, stream=True)
    if r.status_code != 200:
        print(f"Failed to download: {r.status_code}")
        sys.exit(2)

    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)

    print(f"Downloaded to {out_path}")


if __name__ == "__main__":
    main()

