#!/usr/bin/env python3
"""Generate every simulated dataset in ``datasets/sample`` and its metadata.

All data in this repository are synthetic. The generator is the record of what the
numbers mean: concentrations, the true parameters, the error model, and the seed.
Nothing here is copied from a textbook figure, because a figure cannot be re-derived
and a seed can.

Each generator returns ``(rows, columns, notes)`` and declares its own error model:

``proportional_plus_absolute``
    ``v_obs = v_true·(1 + η) + ε`` with ``η ~ N(0, (cv)²)`` and ``ε ~ N(0, sd_abs²)``.
    This is the usual shape of an enzyme assay: a constant relative error from
    pipetting and a floor from the detector.
``absolute``
    ``y_obs = y_true + ε`` with ``ε ~ N(0, sd²)``. Constant variance, which is what
    ordinary least squares assumes — included so the weighting module can show what
    happens when that assumption is false.

Run with ``python scripts/make_datasets.py``. The script is deterministic: running it
twice produces byte-identical CSV files, which ``scripts/check_data_metadata.py``
verifies against the committed SHA-256 sums.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "datasets" / "sample"
METADATA = ROOT / "datasets" / "metadata"

GENERATOR = "scripts/make_datasets.py"

# Units, stated once. Every file uses them; no file converts anything.
UNITS = {
    "S": "µM (substrate)",
    "I": "µM (inhibitor)",
    "v": "µmol min^-1 (mg protein)^-1 (initial rate)",
    "sd": "same units as the value it qualifies",
    "concentration": "µM (standard concentration)",
    "signal": "absorbance units at 340 nm",
    "replicate": "integer label, 1-based",
}


def _mm(s: float, vmax: float, km: float) -> float:
    return vmax * s / (km + s)


def _competitive(s: float, vmax: float, km: float, i: float, ki: float) -> float:
    return vmax * s / (km * (1.0 + i / ki) + s)


def _hill(s: float, vmax: float, k_half: float, n: float) -> float:
    return vmax * s**n / (k_half**n + s**n)


def proportional_plus_absolute(rng, truth, *, cv: float, sd_abs: float):
    """The documented error model: multiplicative noise plus a detector floor."""
    return truth * (1.0 + cv * rng.standard_normal(np.shape(truth))) + sd_abs * rng.standard_normal(
        np.shape(truth)
    )


# ---------------------------------------------------------------------------
# generators
# ---------------------------------------------------------------------------


def michaelis_menten() -> tuple[list[dict], list[str], dict]:
    """A single [S] series with 3 replicates at each substrate concentration."""
    truth = {"vmax": 10.0, "km": 25.0}
    seed = 20261107
    rng = np.random.default_rng(seed)
    S_values = [1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, 200.0]
    cv, sd_abs = 0.03, 0.05
    rows: list[dict] = []
    for s in S_values:
        for rep in (1, 2, 3):
            v_true = _mm(s, truth["vmax"], truth["km"])
            v_obs = float(proportional_plus_absolute(rng, v_true, cv=cv, sd_abs=sd_abs))
            rows.append({"S": s, "replicate": rep, "v": round(v_obs, 5)})
    meta = {
        "title": "Michaelis-Menten substrate series, three replicates",
        "teaches": "the primary fit; replicate spread; why the highest [S] carries the "
        "most information about Vmax and the lowest about Km",
        "true_parameters": truth,
        "error_model": {
            "name": "proportional_plus_absolute",
            "cv": cv,
            "sd_abs": sd_abs,
            "formula": "v_obs = v_true*(1 + N(0, cv^2)) + N(0, sd_abs^2)",
        },
        "seed": seed,
        "substrate_concentrations": S_values,
        "replicates": 3,
    }
    return rows, ["S", "replicate", "v"], meta


def michaelis_menten_mean() -> tuple[list[dict], list[str], dict]:
    """The same data collapsed to a mean and a sample SD per [S] — CLI-ready."""
    rows, _, meta = michaelis_menten()
    by_s: dict[float, list[float]] = {}
    for row in rows:
        by_s.setdefault(row["S"], []).append(row["v"])
    out = []
    for s in sorted(by_s):
        values = np.array(by_s[s])
        sd = float(values.std(ddof=1))
        out.append(
            {"S": s, "v": round(float(values.mean()), 5), "sd": round(sd, 5), "n": len(values)}
        )
    meta = dict(meta)
    meta["title"] = "Michaelis-Menten series, mean and sample SD per substrate concentration"
    meta["teaches"] = (
        "the same data in the shape a fitting tool wants, with the replicate "
        "spread carried as a weight instead of being thrown away"
    )
    meta["derived_from"] = "michaelis_menten.csv (mean of three replicates; SD with n-1)"
    return out, ["S", "v", "sd", "n"], meta


def cooperative() -> tuple[list[dict], list[str], dict]:
    """A sigmoid series where the Hill coefficient is genuinely 2."""
    truth = {"vmax": 12.0, "k_half": 30.0, "n": 2.0}
    seed = 424242
    rng = np.random.default_rng(seed)
    S_values = [2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 100.0, 200.0]
    cv, sd_abs = 0.03, 0.08
    rows = []
    for s in S_values:
        v_true = _hill(s, truth["vmax"], truth["k_half"], truth["n"])
        v_obs = float(proportional_plus_absolute(rng, v_true, cv=cv, sd_abs=sd_abs))
        rows.append({"S": s, "v": round(v_obs, 5)})
    meta = {
        "title": "A sigmoid substrate series with n = 2",
        "teaches": "when a Hill fit earns its third parameter, and that sigmoidicity alone "
        "does not demonstrate subunit cooperativity",
        "true_parameters": truth,
        "error_model": {
            "name": "proportional_plus_absolute",
            "cv": cv,
            "sd_abs": sd_abs,
            "formula": "v_obs = v_true*(1 + N(0, cv^2)) + N(0, sd_abs^2)",
        },
        "seed": seed,
        "substrate_concentrations": S_values,
    }
    return rows, ["S", "v"], meta


def inhibition_series() -> tuple[list[dict], list[str], dict]:
    """A competitive inhibitor at Ki = 8 µM, [I] up to 16 µM: a decisive experiment."""
    truth = {"vmax": 10.0, "km": 25.0, "ki": 8.0, "mechanism": "competitive"}
    seed = 112233
    rng = np.random.default_rng(seed)
    S_values = [2.0, 5.0, 10.0, 25.0, 100.0]
    I_values = [0.0, 2.0, 8.0, 16.0]
    cv, sd_abs = 0.03, 0.05
    rows = []
    for s in S_values:
        for i in I_values:
            v_true = _competitive(s, truth["vmax"], truth["km"], i, truth["ki"])
            v_obs = float(proportional_plus_absolute(rng, v_true, cv=cv, sd_abs=sd_abs))
            rows.append({"S": s, "I": i, "v": round(v_obs, 5)})
    meta = {
        "title": "Competitive inhibition, [I] spanning 0 to 2·Ki",
        "teaches": "typing a mechanism by global fitting; the ΔAICc that justifies calling "
        "competitive the best of the three candidates",
        "true_parameters": truth,
        "error_model": {
            "name": "proportional_plus_absolute",
            "cv": cv,
            "sd_abs": sd_abs,
            "formula": "v_obs = v_true*(1 + N(0, cv^2)) + N(0, sd_abs^2)",
        },
        "seed": seed,
        "substrate_concentrations": S_values,
        "inhibitor_concentrations": I_values,
    }
    return rows, ["S", "I", "v"], meta


def inhibition_weak() -> tuple[list[dict], list[str], dict]:
    """The same mechanism, but a badly designed experiment: [I]/Ki never exceeds 0.004."""
    truth = {"vmax": 10.0, "km": 25.0, "ki": 500.0, "mechanism": "competitive"}
    seed = 998877
    rng = np.random.default_rng(seed)
    S_values = [2.0, 5.0, 10.0, 25.0, 100.0]
    I_values = [0.0, 1.0, 2.0]
    cv, sd_abs = 0.03, 0.05
    rows = []
    for s in S_values:
        for i in I_values:
            v_true = _competitive(s, truth["vmax"], truth["km"], i, truth["ki"])
            v_obs = float(proportional_plus_absolute(rng, v_true, cv=cv, sd_abs=sd_abs))
            rows.append({"S": s, "I": i, "v": round(v_obs, 5)})
    meta = {
        "title": "Weak competitive inhibition that the experiment cannot type",
        "teaches": "the most important negative result in this repository: all three "
        "mechanisms fit equally well, so the honest report is 'unidentified', "
        "not the nominal winner",
        "true_parameters": truth,
        "error_model": {
            "name": "proportional_plus_absolute",
            "cv": cv,
            "sd_abs": sd_abs,
            "formula": "v_obs = v_true*(1 + N(0, cv^2)) + N(0, sd_abs^2)",
        },
        "seed": seed,
        "substrate_concentrations": S_values,
        "inhibitor_concentrations": I_values,
        "known_problem": "[I] never exceeds 0.004·Ki, so Ki is not identified; the file is "
        "here so a learner can watch a model comparison say 'indistinguishable'",
    }
    return rows, ["S", "I", "v"], meta


def calibration_standards() -> tuple[list[dict], list[str], dict]:
    """Standards with proportional error: the case where weighting matters."""
    truth = {"slope": 0.0485, "intercept": 0.011}
    seed = 5150
    rng = np.random.default_rng(seed)
    concentrations = [0.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0]
    cv, sd_blank = 0.04, 0.0025
    rows = []
    for c in concentrations:
        signal_true = truth["slope"] * c + truth["intercept"]
        reps = [
            float(
                signal_true * (1.0 + cv * rng.standard_normal()) + sd_blank * rng.standard_normal()
            )
            for _ in range(3)
        ]
        for rep, value in enumerate(reps, start=1):
            rows.append({"concentration": c, "replicate": rep, "signal": round(value, 6)})
    meta = {
        "title": "Calibration standards, three replicates, strength-proportional error",
        "teaches": "calibration, LOD/LOQ, and why an unweighted fit is the wrong default when "
        "the error grows with the signal",
        "true_parameters": truth,
        "error_model": {
            "name": "proportional_plus_absolute",
            "cv": cv,
            "sd_abs": sd_blank,
            "formula": "signal_obs = signal_true*(1 + N(0, cv^2)) + N(0, sd_blank^2)",
        },
        "seed": seed,
        "standard_concentrations": concentrations,
        "replicates": 3,
    }
    return rows, ["concentration", "replicate", "signal"], meta


def calibration_mean() -> tuple[list[dict], list[str], dict]:
    """The calibration collapsed to mean and SD per standard — CLI-ready."""
    rows, _, meta = calibration_standards()
    by_c: dict[float, list[float]] = {}
    for row in rows:
        by_c.setdefault(row["concentration"], []).append(row["signal"])
    out = []
    for c in sorted(by_c):
        values = np.array(by_c[c])
        out.append(
            {
                "concentration": c,
                "signal": round(float(values.mean()), 6),
                "sd": round(float(values.std(ddof=1)), 6),
            }
        )
    meta = dict(meta)
    meta["title"] = "Calibration standards, mean and sample SD per concentration"
    meta["teaches"] = "LOD and LOQ computed from measured replicate SD rather than a rule of thumb"
    meta["derived_from"] = "calibration_standards.csv (mean and n-1 SD of three replicates)"
    return out, ["concentration", "signal", "sd"], meta


def malformed_example() -> tuple[list[dict], list[str], dict]:
    """A file that is wrong on purpose, so the CLI's refusals can be taught."""
    rows = [
        {"S": 1.0, "v": 1.05},
        {"S": 2.5, "v": 1.92},
        {"S": 5.0, "v": "3.4x"},  # not a number
        {"S": 10.0, "v": -0.4},  # physically impossible
        {"S": 20.0, "v": 6.61},
        {"S": 40.0, "v": 7.98},
    ]
    meta = {
        "title": "A deliberately malformed file",
        "teaches": "the CLI's exit codes: 2 for a file that cannot be read (the '3.4x' cell), "
        "and a fit that refuses negative rates rather than fitting them",
        "true_parameters": {},
        "error_model": {
            "name": "none",
            "formula": "noise added by hand to make the defects obvious to a reader",
        },
        "seed": None,
        "known_problem": "one non-numeric cell and one negative rate; both are intentional",
    }
    return rows, ["S", "v"], meta


GENERATORS = {
    "michaelis_menten": michaelis_menten,
    "michaelis_menten_mean": michaelis_menten_mean,
    "cooperative_hill": cooperative,
    "inhibition_series": inhibition_series,
    "inhibition_weak": inhibition_weak,
    "calibration_standards": calibration_standards,
    "calibration_mean": calibration_mean,
    "malformed_example": malformed_example,
}


def write_dataset(name: str, rows, columns, meta) -> dict:
    csv_path = SAMPLE / f"{name}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    payload = csv_path.read_bytes()
    sha = hashlib.sha256(payload).hexdigest()
    metadata = {
        "name": name,
        "file": f"datasets/sample/{name}.csv",
        "columns": columns,
        "units": {c: UNITS[c] for c in columns if c in UNITS},
        "rows": len(rows),
        "sha256": sha,
        "generator": GENERATOR,
        "synthetic": True,
        "licence": "CC0-1.0 (data generated by this script; no third-party content)",
        **meta,
    }
    yml = METADATA / f"{name}.yml"
    yml.write_text(_to_yaml(metadata), encoding="utf-8")
    return {"name": name, "rows": len(rows), "sha256": sha}


def _to_yaml(obj, indent: int = 0) -> str:
    """Small YAML writer: the metadata is simple, and PyYAML is a dev dependency only."""
    pad = "  " * indent
    if isinstance(obj, dict):
        if not obj:
            return f"{pad}{{}}\n"
        out = []
        for key, value in obj.items():
            if isinstance(value, (dict, list)) and value:
                out.append(f"{pad}{key}:\n{_to_yaml(value, indent + 1)}")
            elif isinstance(value, list):
                out.append(f"{pad}{key}: []\n")
            elif value is None:
                out.append(f"{pad}{key}: null\n")
            elif isinstance(value, bool):
                out.append(f"{pad}{key}: {'true' if value else 'false'}\n")
            elif isinstance(value, (int, float)):
                out.append(f"{pad}{key}: {value}\n")
            else:
                out.append(f"{pad}{key}: {_quote(str(value))}\n")
        return "".join(out)
    if isinstance(obj, list):
        out = []
        for item in obj:
            if isinstance(item, (dict, list)):
                out.append(f"{pad}-\n{_to_yaml(item, indent + 1)}")
            elif isinstance(item, (int, float)) and not isinstance(item, bool):
                out.append(f"{pad}- {item}\n")
            else:
                out.append(f"{pad}- {_quote(str(item))}\n")
        return "".join(out)
    return f"{pad}{_quote(str(obj))}\n"


def _quote(text: str) -> str:
    if text == "" or any(ch in text for ch in ":#{}[]&*!|>'\"%@`\n") or text.strip() != text:
        return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return text


def main() -> int:
    SAMPLE.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)
    written = []
    for name, generator in GENERATORS.items():
        rows, columns, meta = generator()
        written.append(write_dataset(name, rows, columns, meta))
    index = {"datasets": written, "count": len(written)}
    (METADATA / "index.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
    width = max(len(w["name"]) for w in written)
    for item in written:
        print(f"{item['name']:<{width}}  {item['rows']:>3} rows  {item['sha256'][:16]}…")
    print(
        f"\n{len(written)} datasets written to datasets/sample/ with metadata in datasets/metadata/"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
