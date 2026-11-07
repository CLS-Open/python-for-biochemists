"""The committed datasets: metadata, provenance, and what each one is for.

These tests run against the files in ``datasets/sample`` — the same files a learner
opens — so that a regeneration which changes the numbers cannot quietly invalidate a
claim made in the modules.
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from kinetics.calibration import fit_calibration
from kinetics.cli import EXIT_UNREADABLE, main
from kinetics.fitting import compare_linearisations, fit_michaelis_menten, fit_nonlinear
from kinetics.inhibition import compare_inhibition_models
from kinetics.models import hill

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "datasets" / "sample"
METADATA = ROOT / "datasets" / "metadata"

CSV_FILES = sorted(SAMPLE.glob("*.csv"))


def read_rows(name: str) -> list[dict[str, str]]:
    with (SAMPLE / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def metadata(name: str) -> dict:
    return yaml.safe_load((METADATA / f"{Path(name).stem}.yml").read_text(encoding="utf-8"))


def test_datasets_exist():
    assert len(CSV_FILES) >= 6


@pytest.mark.parametrize("path", CSV_FILES, ids=lambda p: p.stem)
def test_every_dataset_declares_itself_synthetic(path):
    meta = metadata(path.name)
    assert meta["synthetic"] is True
    assert meta["generator"] == "scripts/make_datasets.py"
    assert meta["licence"].startswith("CC0")
    assert meta["teaches"]


@pytest.mark.parametrize("path", CSV_FILES, ids=lambda p: p.stem)
def test_row_count_matches_metadata(path):
    assert len(read_rows(path.name)) == metadata(path.name)["rows"]


@pytest.mark.parametrize("path", CSV_FILES, ids=lambda p: p.stem)
def test_every_numeric_column_has_a_unit(path):
    meta = metadata(path.name)
    for column in read_rows(path.name)[0]:
        if column in {"replicate", "n"}:
            continue
        assert column in meta["units"], f"{path.name}: {column} has no declared unit"


def test_the_metadata_checker_passes_on_the_committed_files():
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "check_data_metadata.py")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_index_lists_every_dataset():
    index = json.loads((METADATA / "index.json").read_text(encoding="utf-8"))
    assert index["count"] == len(CSV_FILES)
    assert {item["name"] for item in index["datasets"]} == {p.stem for p in CSV_FILES}


# ---------------------------------------------------------------------------
# what each dataset is supposed to demonstrate
# ---------------------------------------------------------------------------


def test_michaelis_menten_dataset_recovers_its_recorded_truth():
    rows = read_rows("michaelis_menten.csv")
    S = [float(r["S"]) for r in rows]
    v = [float(r["v"]) for r in rows]
    truth = metadata("michaelis_menten.csv")["true_parameters"]
    fit = fit_michaelis_menten(S, v)
    assert fit.params["vmax"] == pytest.approx(truth["vmax"], rel=0.05)
    assert fit.params["km"] == pytest.approx(truth["km"], rel=0.10)
    assert fit.r_squared > 0.99


def test_replicate_dataset_has_three_measurements_per_concentration():
    rows = read_rows("michaelis_menten.csv")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["S"]] = counts.get(row["S"], 0) + 1
    assert set(counts.values()) == {3}


def test_mean_dataset_is_the_mean_of_the_replicates():
    raw = read_rows("michaelis_menten.csv")
    means = read_rows("michaelis_menten_mean.csv")
    by_s: dict[str, list[float]] = {}
    for row in raw:
        by_s.setdefault(row["S"], []).append(float(row["v"]))
    for row in means:
        values = by_s[row["S"]]
        assert float(row["v"]) == pytest.approx(sum(values) / len(values), abs=1e-5)
        assert float(row["n"]) == 3


def test_cooperative_dataset_is_actually_sigmoid():
    rows = read_rows("cooperative_hill.csv")
    S = [float(r["S"]) for r in rows]
    v = [float(r["v"]) for r in rows]
    truth = metadata("cooperative_hill.csv")["true_parameters"]
    hill_fit = fit_nonlinear(S, v, hill, (12.0, 30.0, 1.0), ("vmax", "k_half", "n"), name="hill")
    assert hill_fit.params["n"] == pytest.approx(truth["n"], rel=0.25)
    mm = fit_michaelis_menten(S, v)
    assert hill_fit.aicc < mm.aicc  # the sigmoid data prefer the Hill model


def test_inhibition_dataset_types_as_competitive():
    rows = [
        {"S": float(r["S"]), "I": float(r["I"]), "v": float(r["v"])}
        for r in read_rows("inhibition_series.csv")
    ]
    comparison = compare_inhibition_models(rows)
    assert comparison.best == "competitive"
    assert comparison.delta_aicc["noncompetitive"] > 10.0


def test_weak_inhibition_dataset_is_the_one_that_cannot_be_typed():
    rows = [
        {"S": float(r["S"]), "I": float(r["I"]), "v": float(r["v"])}
        for r in read_rows("inhibition_weak.csv")
    ]
    comparison = compare_inhibition_models(rows)
    assert comparison.delta_aicc[comparison.best] == 0.0
    assert min(d for m, d in comparison.delta_aicc.items() if m != comparison.best) < 2.0
    assert any("indistinguishable" in note for note in comparison.notes)


def test_calibration_dataset_supports_both_weightings_and_they_disagree():
    rows = read_rows("calibration_mean.csv")
    x = [float(r["concentration"]) for r in rows]
    y = [float(r["signal"]) for r in rows]
    sd = [float(r["sd"]) for r in rows]
    weighted = fit_calibration(x, y, weighted=True, replicate_sd=sd)
    unweighted = fit_calibration(x, y)
    truth = metadata("calibration_standards.csv")["true_parameters"]
    assert abs(weighted.slope - truth["slope"]) < abs(unweighted.slope - truth["slope"])
    assert weighted.lod < unweighted.lod


def test_malformed_dataset_is_rejected_by_the_cli():
    assert main(["mm", str(SAMPLE / "malformed_example.csv")]) == EXIT_UNREADABLE


def test_the_linearisation_trap_is_visible_in_the_committed_data():
    """On this dataset the nonlinear fit is closer to the truth than Lineweaver-Burk."""
    rows = read_rows("michaelis_menten_mean.csv")
    S = [float(r["S"]) for r in rows]
    v = [float(r["v"]) for r in rows]
    results = {r.model: r for r in compare_linearisations(S, v)}
    truth = 10.0
    assert abs(results["michaelis_menten"].params["vmax"] - truth) < abs(
        results["lineweaver_burk"].params["vmax"] - truth
    )
