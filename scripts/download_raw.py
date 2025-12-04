#!/usr/bin/env python3
import yaml
import requests
import os
from pathlib import Path


def load_datasets():
    with open("data/datasets.yml", "r") as f:
        return yaml.safe_load(f)["datasets"]


def download_raw(dataset):
    name = dataset["name"]
    source = dataset["source"]
    raw_dir = Path("data/raw")
    raw_dir.mkdir(exist_ok=True)
    filename = source.split("/")[-1]
    filepath = raw_dir / filename
    if not filepath.exists():
        response = requests.get(source, headers={"User-Agent": "Mozilla/5.0"})
        response.raise_for_status()
        with open(filepath, "wb") as f:
            f.write(response.content)
        print(f"Downloaded {filename}")
    else:
        print(f"{filename} already exists")


if __name__ == "__main__":
    datasets = load_datasets()
    for ds in datasets:
        download_raw(ds)
