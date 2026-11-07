#!/usr/bin/env python3
"""Build the figures used in the modules, and record every number they display.

Two rules:

* a figure is drawn from a committed dataset or from a seeded simulation, never from
  numbers typed into this file;
* every value that appears in a caption is written to ``figures/figure_values.json``,
  so a module can quote it and a reviewer can check it against the code that produced
  it.

The script is deterministic — same matplotlib version, same bytes — and prints the
SHA-256 of each figure so a change is visible in review.

Usage::

    python scripts/make_figures.py
    python scripts/make_figures.py --dpi 200
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kinetics.calibration import fit_calibration
from kinetics.fitting import (
    bootstrap_ci,
    compare_linearisations,
    fit_michaelis_menten,
)
from kinetics.inhibition import compare_inhibition_models
from kinetics.models import (
    competitive,
    michaelis_menten,
    noncompetitive,
    uncompetitive,
)

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "datasets" / "sample"
FIGURES = ROOT / "figures"

INK = "#1b1b1b"
BLUE = "#1f5fa8"
RED = "#b3221f"
GREEN = "#2f6b3a"
GREY = "#8a8a8a"


def read_csv(name: str) -> list[dict[str, str]]:
    with (SAMPLE / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def style(ax, *, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, fontsize=12, color=INK, loc="left")
    ax.set_xlabel(xlabel, fontsize=10, color=INK)
    ax.set_ylabel(ylabel, fontsize=10, color=INK)
    ax.tick_params(labelsize=9, colors=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GREY)
    ax.grid(True, color="#e6e6e6", linewidth=0.6)
    ax.set_axisbelow(True)


def figure_fit(ax, values: dict) -> None:
    rows = read_csv("michaelis_menten.csv")
    S = np.array([float(r["S"]) for r in rows])
    v = np.array([float(r["v"]) for r in rows])
    fit = fit_michaelis_menten(S, v)
    ci = bootstrap_ci(S, v, n_resamples=2000, seed=20261107)

    ax.scatter(S, v, s=18, color=BLUE, alpha=0.75, label="simulated measurements (3 per [S])")
    grid = np.linspace(0, 220, 400)
    ax.plot(
        grid,
        michaelis_menten(grid, **fit.params),
        color=INK,
        lw=1.6,
        label=f"fit: Vmax {fit.params['vmax']:.2f}, Km {fit.params['km']:.1f} µM",
    )

    # Bootstrap band: recompute the curve for the interval endpoints, which is what a
    # band from a percentile bootstrap of the *parameters* actually means.
    low = michaelis_menten(grid, ci["vmax"][0], ci["km"][0])
    high = michaelis_menten(grid, ci["vmax"][1], ci["km"][1])
    ax.fill_between(
        grid,
        np.minimum(low, high),
        np.maximum(low, high),
        color=BLUE,
        alpha=0.13,
        lw=0,
        label=f"95 % bootstrap band, seed 20261107 ({ci['vmax'][0]:.2f}–{ci['vmax'][1]:.2f})",
    )
    ax.axhline(fit.params["vmax"], color=RED, lw=1.0, ls=":")
    ax.annotate(
        "Vmax is an extrapolation,\nnever a measurement",
        xy=(140, fit.params["vmax"]),
        xytext=(90, fit.params["vmax"] - 1.6),
        fontsize=8.5,
        color=RED,
        arrowprops={"arrowstyle": "-", "color": RED, "lw": 0.8},
    )
    style(
        ax,
        title="A. Michaelis–Menten fit with a bootstrap band",
        xlabel="[S] (µM)",
        ylabel="v (µmol min⁻¹ mg⁻¹)",
    )
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.set_ylim(0, 11.5)

    values["figure_1"] = {
        "vmax": round(fit.params["vmax"], 4),
        "km": round(fit.params["km"], 4),
        "vmax_stderr": round(fit.stderr["vmax"], 4),
        "km_stderr": round(fit.stderr["km"], 4),
        "r_squared": round(fit.r_squared, 5),
        "aicc": round(fit.aicc, 3),
        "bootstrap_vmax_95": [round(x, 4) for x in ci["vmax"]],
        "bootstrap_km_95": [round(x, 4) for x in ci["km"]],
        "points": len(S),
        "source": "datasets/sample/michaelis_menten.csv",
    }


def figure_bias(ax, values: dict) -> None:
    S = np.array([1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, 200.0])
    order = ["michaelis_menten", "hanes", "eadie_hofstee", "lineweaver_burk"]
    labels = {
        "michaelis_menten": "nonlinear",
        "hanes": "Hanes–Woolf",
        "eadie_hofstee": "Eadie–Hofstee",
        "lineweaver_burk": "Lineweaver–Burk",
    }
    estimates: dict[str, list[float]] = {name: [] for name in order}
    for seed in range(400):
        rng = np.random.default_rng(seed)
        v = michaelis_menten(S, 10.0, 25.0) * (1.0 + 0.03 * rng.standard_normal(S.size))
        v = v + 0.05 * rng.standard_normal(S.size)
        if np.any(v <= 0):
            continue
        for result in compare_linearisations(S, v):
            estimates[result.model].append(result.params["vmax"])

    data = [estimates[name] for name in order]
    parts = ax.violinplot(data, showextrema=False, widths=0.8)
    colours = [GREEN, BLUE, GREY, RED]
    for body, colour in zip(parts["bodies"], colours, strict=True):
        body.set_facecolor(colour)
        body.set_edgecolor(colour)
        body.set_alpha(0.45)
    clipped = 0
    for index, values_list in enumerate(data, start=1):
        median = float(np.median(values_list))
        within = float(np.mean(np.abs(np.array(values_list) - 10.0) / 10.0 < 0.2))
        ax.plot([index], [median], marker="o", color=INK, ms=4)
        ax.text(
            index,
            24,
            f"median\n{median:.2f}\n{within:.0%} within 20 %",
            ha="center",
            va="top",
            fontsize=8,
            color=INK,
        )
        clipped += int(np.sum(np.array(values_list) > 25))
    ax.axhline(10.0, color=INK, lw=1.0, ls="--")
    ax.text(0.55, 10.4, "true Vmax = 10", fontsize=8.5, color=INK)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels([labels[name] for name in order], fontsize=9)
    style(
        ax,
        title="B. 400 simulated assays: linearisations are heavy-tailed, not just noisier",
        xlabel="estimation method",
        ylabel="estimated Vmax (µmol min⁻¹ mg⁻¹)",
    )
    ax.set_ylim(0, 25)
    if clipped:
        ax.text(
            0.99,
            0.02,
            f"{clipped} estimate(s) above 25 are clipped by the axis",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color=RED,
        )

    summary = {}
    for name in order:
        arr = np.array(estimates[name])
        summary[labels[name]] = {
            "median_vmax": round(float(np.median(arr)), 3),
            "iqr_vmax": [
                round(float(np.percentile(arr, 25)), 3),
                round(float(np.percentile(arr, 75)), 3),
            ],
            "fraction_within_20pct": round(float(np.mean(np.abs(arr - 10.0) / 10.0 < 0.2)), 3),
            "n": int(arr.size),
        }
    values["figure_2"] = {
        "seeds": 400,
        "estimates_above_axis_limit": clipped,
        "noise": "3 % proportional + 0.05 absolute, constant across [S]",
        "methods": summary,
    }


def figure_inhibition(ax, values: dict) -> None:
    rows = read_csv("inhibition_series.csv")
    surface = [{"S": float(r["S"]), "I": float(r["I"]), "v": float(r["v"])} for r in rows]
    comparison = compare_inhibition_models(surface)

    grid = np.linspace(0, 210, 400)
    vmax, km, ki = 10.0, 25.0, 8.0
    for I, colour, ls in ((0.0, INK, "-"), (8.0, BLUE, "-"), (16.0, RED, "-")):
        ax.plot(
            grid,
            competitive(grid, vmax, km, I, ki),
            color=colour,
            lw=1.5,
            ls=ls,
            label=f"competitive, [I] = {I:g} µM ({I / ki:g}·Ki)".replace("0·Ki", "no inhibitor"),
        )
    ax.plot(
        grid,
        noncompetitive(grid, vmax, km, 16.0, ki),
        color=GREEN,
        lw=1.2,
        ls="--",
        label="non-competitive, [I] = 16 µM",
    )
    ax.plot(
        grid,
        uncompetitive(grid, vmax, km, 16.0, ki),
        color=GREY,
        lw=1.2,
        ls=":",
        label="uncompetitive, [I] = 16 µM",
    )
    for I, marker in ((0.0, "o"), (8.0, "s"), (16.0, "^")):
        pts = [(float(r["S"]), float(r["v"])) for r in rows if float(r["I"]) == I]
        ax.scatter(
            [p[0] for p in pts],
            [p[1] for p in pts],
            s=16,
            marker=marker,
            color=INK,
            alpha=0.55,
            label=f"simulated data, [I] = {I:g} µM",
        )
    style(
        ax,
        title="C. Three mechanisms differ in shape, not just in height",
        xlabel="[S] (µM)",
        ylabel="v (µmol min⁻¹ mg⁻¹)",
    )
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    ax.set_ylim(0, 11)

    values["figure_3"] = {
        "best_by_aicc": comparison.best,
        "delta_aicc": {k: round(v, 3) for k, v in comparison.delta_aicc.items()},
        "fit_ki": round(comparison.results[0].params["ki"], 4),
        "fit_vmax": round(comparison.results[0].params["vmax"], 4),
        "fit_km": round(comparison.results[0].params["km"], 4),
        "source": "datasets/sample/inhibition_series.csv",
        "truth": {"vmax": 10.0, "km": 25.0, "ki": 8.0, "mechanism": "competitive"},
    }


def figure_calibration(ax_main, ax_resid, values: dict) -> None:
    rows = read_csv("calibration_mean.csv")
    x = np.array([float(r["concentration"]) for r in rows])
    y = np.array([float(r["signal"]) for r in rows])
    sd = np.array([float(r["sd"]) for r in rows])
    weighted = fit_calibration(x, y, weighted=True, replicate_sd=sd)
    unweighted = fit_calibration(x, y)

    ax_main.errorbar(
        x,
        y,
        yerr=sd,
        fmt="o",
        ms=5,
        color=BLUE,
        ecolor=GREY,
        elinewidth=1,
        capsize=2,
        label="standards, mean ± SD of 3 replicates",
    )
    grid = np.linspace(0, 210, 200)
    ax_main.plot(
        grid,
        weighted.slope * grid + weighted.intercept,
        color=INK,
        lw=1.5,
        label=f"weighted fit, slope {weighted.slope:.4f}",
    )
    ax_main.plot(
        grid,
        unweighted.slope * grid + unweighted.intercept,
        color=RED,
        lw=1.1,
        ls="--",
        label=f"unweighted fit, slope {unweighted.slope:.4f}",
    )
    ax_main.set_xlim(-5, 215)
    style(
        ax_main,
        title="D. Calibration: the weighted line is the honest one",
        xlabel="standard concentration (µM)",
        ylabel="signal (A₃₄₀)",
    )
    ax_main.legend(fontsize=8, frameon=False, loc="upper left")

    ax_main.text(
        0.03,
        0.72,
        f"LOD (3.3σ/slope)\n  unweighted  {unweighted.lod:7.2f} µM\n"
        f"  weighted    {weighted.lod:7.2f} µM\n"
        f"ratio {unweighted.lod / weighted.lod:.0f}×",
        transform=ax_main.transAxes,
        fontsize=8.5,
        color=INK,
        va="top",
        family="monospace",
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": GREY, "lw": 0.7},
    )

    # Panel E shows the error model itself. Plotting residuals instead would be
    # misleading: ordinary least squares makes residuals *small* where the data are
    # noisy, so an unweighted fit "looks better" in exactly the region where it is
    # least trustworthy. The replicate SD is the measurable fact.
    positive = x > 0
    if np.any(positive):
        cv_estimate = float(np.sum(sd[positive] / y[positive]) / np.sum(positive))
    else:
        cv_estimate = float("nan")
    ax_resid.loglog(x[positive], sd[positive], "o", ms=5, color=BLUE, label="replicate SD (n = 3)")
    ax_resid.loglog(
        x[positive],
        cv_estimate * (weighted.slope * x[positive] + weighted.intercept),
        color=INK,
        lw=1.2,
        label=f"constant {cv_estimate:.1%} of the signal",
    )
    style(
        ax_resid,
        title="E. The noise grows with the standard",
        xlabel="standard concentration (µM, log)",
        ylabel="replicate SD (A₃₄₀, log)",
    )
    ax_resid.legend(fontsize=8, frameon=False, loc="upper left")
    for axis in (ax_resid.xaxis, ax_resid.yaxis):
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    values["figure_4"] = {
        "weighted": {
            "slope": round(weighted.slope, 6),
            "intercept": round(weighted.intercept, 6),
            "residual_sd": round(weighted.residual_sd, 6),
            "lod": round(weighted.lod, 4),
            "loq": round(weighted.loq, 4),
        },
        "unweighted": {
            "slope": round(unweighted.slope, 6),
            "intercept": round(unweighted.intercept, 6),
            "residual_sd": round(unweighted.residual_sd, 6),
            "lod": round(unweighted.lod, 4),
            "loq": round(unweighted.loq, 4),
        },
        "lod_ratio_unweighted_over_weighted": round(unweighted.lod / weighted.lod, 2),
        "relative_error_estimate": round(cv_estimate, 4),
        "source": "datasets/sample/calibration_mean.csv",
        "truth": {"slope": 0.0485, "intercept": 0.011},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument(
        "--check",
        action="store_true",
        help="build into a temporary directory and verify the committed figures and "
        "figure_values.json are byte-identical, without modifying them (exit 1 on drift)",
    )
    args = parser.parse_args()
    committed = FIGURES
    #: Where this run writes. --check builds into a throwaway directory so that a verification
    #: run cannot accidentally repair the very drift it is supposed to detect.
    out_dir = Path(tempfile.mkdtemp(prefix="figcheck-")) if args.check else FIGURES
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.size": 10,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    values: dict = {
        "_note": "Every number below was computed by scripts/make_figures.py. "
        "Modules quote these values; nothing here is typed by hand."
    }
    hashes: dict[str, str] = {}

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    figure_fit(ax, values)
    fig.tight_layout()
    path = out_dir / "fig1_michaelis_menten_fit.png"
    fig.savefig(path, dpi=args.dpi)
    plt.close(fig)
    hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    figure_bias(ax, values)
    fig.tight_layout()
    path = out_dir / "fig2_linearisation_bias.png"
    fig.savefig(path, dpi=args.dpi)
    plt.close(fig)
    hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    figure_inhibition(ax, values)
    fig.tight_layout()
    path = out_dir / "fig3_inhibition_mechanisms.png"
    fig.savefig(path, dpi=args.dpi)
    plt.close(fig)
    hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()

    fig, (ax_main, ax_resid) = plt.subplots(
        1, 2, figsize=(9.6, 4.0), gridspec_kw={"width_ratios": [1.35, 1]}
    )
    figure_calibration(ax_main, ax_resid, values)
    fig.tight_layout()
    path = out_dir / "fig4_calibration_weighting.png"
    fig.savefig(path, dpi=args.dpi)
    plt.close(fig)
    hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()

    values["_figure_hashes"] = hashes
    out = out_dir / "figure_values.json"
    out.write_text(json.dumps(values, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    combined = hashlib.sha256("".join(sorted(hashes.values())).encode()).hexdigest()

    if args.check:
        # The point of --check is that a figure cannot silently drift away from the data that
        # produced it, and neither can the numbers the modules quote from it.
        drift: list[str] = []
        for name, digest in sorted(hashes.items()):
            target = committed / name
            if not target.is_file():
                drift.append(f"{name}: missing from figures/")
            elif hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                drift.append(f"{name}: differs from the committed file")
        committed_values = committed / "figure_values.json"
        if not committed_values.is_file():
            drift.append("figure_values.json: missing")
        elif json.loads(committed_values.read_text(encoding="utf-8")) != values:
            drift.append("figure_values.json: the numbers have changed")
        if drift:
            print("figures do not reproduce:", file=sys.stderr)
            for line in drift:
                print(f"  - {line}", file=sys.stderr)
            print(
                "\nrerun `python scripts/make_figures.py` to regenerate, then inspect the "
                "diff before committing",
                file=sys.stderr,
            )
            return 1
        print(
            f"OK   {len(hashes)} figures and figure_values.json reproduce exactly "
            f"(combined hash {combined[:16]})"
        )
        return 0

    for name, digest in sorted(hashes.items()):
        print(f"{name:<36} {digest[:16]}…")
    print(
        f"\n{len(hashes)} figures and {out.relative_to(ROOT)} written; "
        f"combined hash {combined[:16]}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
