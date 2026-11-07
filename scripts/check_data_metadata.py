#!/usr/bin/env python3
"""Verify that every simulated dataset still matches its metadata.

Checks, one per line of output, for every file in ``datasets/sample``:

1. a metadata YAML exists for it and names the same file;
2. the SHA-256 in the metadata matches the bytes on disk;
3. the row count in the metadata matches the file;
4. the metadata declares ``synthetic: true`` — this repository has no measured data,
   and a dataset that does not say so would be indistinguishable from a copied one;
5. the generator recorded in the metadata exists and is the script that is supposed
   to produce it;
6. a seed is recorded whenever the error model adds noise (a simulated dataset with
   no seed cannot be regenerated, only trusted);
7. every column has a declared unit (or is a label such as ``replicate``);
8. no CSV is orphaned and no metadata points at a missing file.

Exit code 0 if everything matches, 1 otherwise. ``--quiet`` prints only failures.

Usage::

    python scripts/check_data_metadata.py
    python scripts/check_data_metadata.py --quiet
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "datasets" / "sample"
METADATA = ROOT / "datasets" / "metadata"

LABEL_COLUMNS = {"replicate", "n"}
GENERATOR = Path("scripts/make_datasets.py")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_one(csv_path: Path) -> list[str]:
    problems: list[str] = []
    meta_path = METADATA / f"{csv_path.stem}.yml"
    if not meta_path.exists():
        return [f"{csv_path.name}: no metadata file at datasets/metadata/{csv_path.stem}.yml"]

    try:
        meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        return [f"{meta_path.name}: not valid YAML ({exc})"]
    if not isinstance(meta, dict):
        return [f"{meta_path.name}: metadata is not a mapping"]

    if meta.get("file") != f"datasets/sample/{csv_path.name}":
        problems.append(
            f"{meta_path.name}: 'file' is {meta.get('file')!r}, expected "
            f"'datasets/sample/{csv_path.name}'"
        )

    actual_sha = sha256(csv_path)
    if meta.get("sha256") != actual_sha:
        problems.append(
            f"{csv_path.name}: sha256 {actual_sha[:16]}… does not match metadata "
            f"{str(meta.get('sha256'))[:16]}… — regenerate or update the metadata"
        )

    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if meta.get("rows") != len(rows):
        problems.append(f"{csv_path.name}: {len(rows)} data rows, metadata says {meta.get('rows')}")

    if meta.get("synthetic") is not True:
        problems.append(f"{csv_path.name}: metadata must declare 'synthetic: true'")

    generator = meta.get("generator")
    if generator != str(GENERATOR):
        problems.append(f"{csv_path.name}: generator is {generator!r}, expected {str(GENERATOR)!r}")
    elif not (ROOT / GENERATOR).exists():
        problems.append(f"{csv_path.name}: generator {GENERATOR} does not exist")

    error_model = meta.get("error_model") or {}
    adds_noise = error_model.get("name") not in (None, "none")
    if adds_noise and meta.get("seed") is None:
        problems.append(f"{csv_path.name}: a noisy dataset with no seed cannot be regenerated")

    units = meta.get("units") or {}
    header = list(rows[0].keys()) if rows else list(meta.get("columns") or [])
    for column in header:
        if column in LABEL_COLUMNS:
            continue
        if column not in units:
            problems.append(f"{csv_path.name}: column {column!r} has no declared unit")

    if meta.get("teaches") in (None, ""):
        problems.append(f"{csv_path.name}: metadata does not say what the dataset teaches")

    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quiet", action="store_true", help="print only failures")
    args = parser.parse_args()

    csv_files = sorted(SAMPLE.glob("*.csv"))
    yml_files = sorted(METADATA.glob("*.yml"))
    failures: list[str] = []

    if not csv_files:
        print(f"FAIL datasets/sample is empty (looked in {SAMPLE})")
        return 1

    for csv_path in csv_files:
        problems = check_one(csv_path)
        failures.extend(problems)
        if not args.quiet:
            status = "OK  " if not problems else "FAIL"
            print(f"{status} {csv_path.name:<32} {len(problems)} problem(s)")

    csv_stems = {p.stem for p in csv_files}
    for meta_path in yml_files:
        if meta_path.stem not in csv_stems:
            failures.append(f"{meta_path.name}: metadata without a dataset file")

    extra = [p.name for p in SAMPLE.iterdir() if p.suffix not in {".csv", ".md"}]
    if extra:
        failures.append(f"unexpected files in datasets/sample: {extra}")

    if failures:
        print()
        for problem in failures:
            print(f"FAIL {problem}")
        print(f"\n{len(failures)} problem(s) across {len(csv_files)} dataset(s)")
        return 1

    total_rows = sum(len(list(csv.DictReader(p.open(encoding="utf-8")))) for p in csv_files)
    print(
        f"\n{len(csv_files)}/{len(csv_files)} datasets match their metadata "
        f"({total_rows} rows, all synthetic, all sha256 verified)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
