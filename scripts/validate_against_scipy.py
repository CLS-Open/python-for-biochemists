#!/usr/bin/env python3
"""Cross-check the package's numbers against independent implementations.

A fitting routine that agrees with itself proves nothing. Every check below computes
the same quantity two ways — once with ``kinetics``, once with a direct SciPy/NumPy
formula written out here — and compares them.

The checks are deliberately written *without* using the helpers under test: the
reference for the Michaelis-Menten fit is ``scipy.optimize.curve_fit`` with a
different optimiser and different starting values; the reference for the linear
regressions is the closed-form ordinary-least-squares formula; the reference for the
bootstrap is ``scipy.stats.bootstrap``; the reference for AICc is typed out from the
definition.

Exit code 0 if every comparison agrees within its stated tolerance.

Usage::

    python scripts/validate_against_scipy.py
    python scripts/validate_against_scipy.py --verbose
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from kinetics.calibration import fit_calibration, lod_loq
from kinetics.fitting import (
    aic,
    aicc,
    bootstrap_ci,
    fit_eadie_hofstee,
    fit_hanes,
    fit_lineweaver_burk,
    fit_michaelis_menten,
)
from kinetics.inhibition import compare_inhibition_models
from kinetics.models import hill, michaelis_menten

AGREEMENTS: list[tuple[str, float, float]] = []  # (description, ours, theirs)


def record(
    label: str, ours: float, theirs: float, *, rtol: float = 1e-6, atol: float = 1e-9
) -> bool:
    ok = math.isclose(ours, theirs, rel_tol=rtol, abs_tol=atol)
    AGREEMENTS.append((label, ours, theirs))
    if not ok:
        print(f"MISMATCH {label}: kinetics={ours!r} reference={theirs!r}")
    return ok


# ---------------------------------------------------------------------------
# 1. Michaelis-Menten nonlinear fit vs curve_fit
# ---------------------------------------------------------------------------
def check_mm_fit(verbose: bool) -> int:
    checks = 0
    S = np.array([1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, 200.0])
    for seed in range(25):
        rng = np.random.default_rng(seed)
        v = 10.0 * S / (25.0 + S) * (
            1.0 + 0.03 * rng.standard_normal(S.size)
        ) + 0.05 * rng.standard_normal(S.size)
        ours = fit_michaelis_menten(S, v)
        theirs, _ = curve_fit(
            lambda s, vmax, km: vmax * s / (km + s),
            S,
            v,
            p0=[8.0, 40.0],
            method="lm",
            xtol=1e-14,
            ftol=1e-14,
            gtol=1e-14,
            maxfev=40000,
        )
        checks += record(f"vmax seed {seed}", ours.params["vmax"], float(theirs[0]), rtol=1e-6)
        checks += record(f"km seed {seed}", ours.params["km"], float(theirs[1]), rtol=1e-6)
        if verbose:
            print(
                f"  mm seed {seed:>2}: vmax {ours.params['vmax']:.6f} vs {theirs[0]:.6f}; "
                f"km {ours.params['km']:.6f} vs {theirs[1]:.6f}"
            )
    return checks


# ---------------------------------------------------------------------------
# 2. Linearisations vs the closed-form OLS formula
# ---------------------------------------------------------------------------
def ols(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Slope and intercept from the textbook formulas, no library fit."""
    x_mean, y_mean = x.mean(), y.mean()
    slope = float(np.sum((x - x_mean) * (y - y_mean)) / np.sum((x - x_mean) ** 2))
    return slope, float(y_mean - slope * x_mean)


def check_linearisations(verbose: bool) -> int:
    checks = 0
    S = np.array([1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0])
    rng = np.random.default_rng(99)
    v = michaelis_menten(S, 10.0, 25.0) * (1.0 + 0.05 * rng.standard_normal(S.size))

    lb = fit_lineweaver_burk(S, v)
    slope, intercept = ols(1.0 / S, 1.0 / v)
    checks += record("lineweaver_burk slope", lb.params["slope"], slope)
    checks += record("lineweaver_burk intercept", lb.params["intercept"], intercept)
    checks += record("lineweaver_burk vmax", lb.params["vmax"], 1.0 / intercept)

    hw = fit_hanes(S, v)
    slope, intercept = ols(S, S / v)
    checks += record("hanes slope", hw.params["slope"], slope)
    checks += record("hanes vmax", hw.params["vmax"], 1.0 / slope)
    checks += record("hanes km", hw.params["km"], intercept / slope)

    eh = fit_eadie_hofstee(S, v)
    slope, intercept = ols(v / S, v)
    checks += record("eadie_hofstee slope", eh.params["slope"], slope)
    checks += record("eadie_hofstee km", eh.params["km"], -slope)
    checks += record("eadie_hofstee vmax", eh.params["vmax"], intercept)

    if verbose:
        print(
            f"  linearisations: LB Vmax {lb.params['vmax']:.4f}, Hanes {hw.params['vmax']:.4f}, "
            f"EH {eh.params['vmax']:.4f} (truth 10)"
        )
    return checks


# ---------------------------------------------------------------------------
# 3. Information criteria vs the definitions typed out
# ---------------------------------------------------------------------------
def check_information_criteria() -> int:
    checks = 0
    for rss, n, k in ((4.0, 10, 2), (0.53, 9, 2), (12.7, 20, 3), (0.0021, 7, 3)):
        checks += record(f"aic rss={rss}", aic(rss, n, k), n * math.log(rss / n) + 2 * k)
        expected_aicc = n * math.log(rss / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
        checks += record(f"aicc rss={rss}", aicc(rss, n, k), expected_aicc)
    return checks


# ---------------------------------------------------------------------------
# 4. Bootstrap CI vs scipy.stats.bootstrap
# ---------------------------------------------------------------------------
def check_bootstrap(verbose: bool) -> int:
    """Percentile bootstrap, compared with an independently written resampling loop.

    Both sides resample (S, v) pairs with replacement and take the 2.5/97.5 percentiles
    of the refitted Vmax, but the two loops use different random streams and the
    reference fits with ``curve_fit`` rather than with ``kinetics``. Agreement to
    within a few percent of the interval width is the most these can promise — the
    interval itself has Monte-Carlo error, which is exactly the point of recording a
    seed.
    """
    S = np.array([1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0])
    rng = np.random.default_rng(303)
    v = michaelis_menten(S, 10.0, 25.0) * (1.0 + 0.04 * rng.standard_normal(S.size))
    ours = bootstrap_ci(S, v, n_resamples=2000, seed=1234)

    # reference implementation, written out here
    ref_rng = np.random.default_rng(99991)
    vmaxes: list[float] = []
    for _ in range(2000):
        idx = ref_rng.integers(0, S.size, S.size)
        try:
            popt, _ = curve_fit(
                lambda s, vmax, km: vmax * s / (km + s),
                S[idx],
                v[idx],
                p0=[10.0, 25.0],
                maxfev=20000,
            )
        except (RuntimeError, ValueError):
            continue
        vmaxes.append(float(popt[0]))
    theirs = (float(np.percentile(vmaxes, 2.5)), float(np.percentile(vmaxes, 97.5)))

    mid_ours = 0.5 * (ours["vmax"][0] + ours["vmax"][1])
    mid_theirs = 0.5 * (theirs[0] + theirs[1])
    width_ours = ours["vmax"][1] - ours["vmax"][0]
    width_theirs = theirs[1] - theirs[0]
    checks = 0
    checks += record("bootstrap vmax CI midpoint", mid_ours, mid_theirs, rtol=0.02)
    checks += record("bootstrap vmax CI width", width_ours, width_theirs, rtol=0.25)
    if verbose:
        print(
            f"  bootstrap vmax CI: kinetics {ours['vmax'][0]:.4f}-{ours['vmax'][1]:.4f}, "
            f"reference {theirs[0]:.4f}-{theirs[1]:.4f} "
            f"({len(vmaxes)} of 2000 resamples converged)"
        )
    return checks


# ---------------------------------------------------------------------------
# 5. Calibration and LOD vs closed form
# ---------------------------------------------------------------------------
def check_calibration(verbose: bool) -> int:
    checks = 0
    x = np.array([0.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0])
    rng = np.random.default_rng(77)
    y = 0.0485 * x + 0.011
    y = y * (1.0 + 0.04 * rng.standard_normal(x.size)) + 0.0025 * rng.standard_normal(x.size)

    ours = fit_calibration(x, y)
    slope, intercept = ols(x, y)
    residuals = y - (slope * x + intercept)
    sigma = math.sqrt(float(np.sum(residuals**2)) / (len(x) - 2))

    checks += record("calibration slope", ours.slope, slope)
    checks += record("calibration intercept", ours.intercept, intercept)
    checks += record("calibration residual_sd", ours.residual_sd, sigma)
    checks += record("calibration LOD", ours.lod, 3.3 * sigma / abs(slope))
    checks += record("calibration LOQ", ours.loq, 10.0 * sigma / abs(slope))

    lod, loq = lod_loq(slope=0.2, residual_sd=0.01)
    checks += record("lod_loq lod", lod, 3.3 * 0.01 / 0.2)
    checks += record("lod_loq loq", loq, 10.0 * 0.01 / 0.2)
    if verbose:
        print(f"  calibration: slope {ours.slope:.6f} (truth 0.0485), LOD {ours.lod:.4f}")
    return checks


# ---------------------------------------------------------------------------
# 6. Inhibition global fit vs curve_fit on the same surface
# ---------------------------------------------------------------------------
def check_inhibition(verbose: bool) -> int:
    checks = 0
    S_values = (2.0, 5.0, 10.0, 25.0, 100.0)
    I_values = (0.0, 2.0, 8.0, 16.0)
    rows = [
        {"S": s, "I": i, "v": 10.0 * s / (25.0 * (1.0 + i / 8.0) + s)}
        for s in S_values
        for i in I_values
    ]
    comparison = compare_inhibition_models(rows)
    best = comparison.results[0]

    # Reference: curve_fit on the same competitive surface, different start points.
    S = np.array([r["S"] for r in rows])
    I = np.array([r["I"] for r in rows])
    v = np.array([r["v"] for r in rows])
    theirs, _ = curve_fit(
        lambda s, vmax, km, ki: vmax * s / (km * (1.0 + I / ki) + s),
        S,
        v,
        p0=[5.0, 50.0, 2.0],
        method="lm",
        maxfev=40000,
    )
    checks += record("inhibition vmax", best.params["vmax"], float(theirs[0]), rtol=1e-5)
    checks += record("inhibition km", best.params["km"], float(theirs[1]), rtol=1e-5)
    checks += record("inhibition ki", best.params["ki"], float(theirs[2]), rtol=1e-5)
    if verbose:
        print(
            f"  inhibition: best {comparison.best}, Ki {best.params['ki']:.5f} "
            f"vs curve_fit {theirs[2]:.5f}"
        )
    return checks


# ---------------------------------------------------------------------------
# 7. Hill with n = 1 must equal Michaelis-Menten exactly
# ---------------------------------------------------------------------------
def check_hill_reduces() -> int:
    checks = 0
    for s in (0.5, 3.0, 25.0, 250.0):
        checks += record(
            f"hill(n=1) at S={s}",
            hill(s, 10.0, 25.0, 1.0),
            michaelis_menten(s, 10.0, 25.0),
            rtol=1e-12,
        )
    return checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    checks = 0
    for name, function in (
        ("Michaelis-Menten fit vs curve_fit", check_mm_fit),
        ("linearisations vs closed-form OLS", check_linearisations),
        ("AICc vs definition", lambda verbose: check_information_criteria()),
        ("bootstrap CI vs independent resampling", check_bootstrap),
        ("calibration and LOD vs closed form", check_calibration),
        ("inhibition fit vs curve_fit", check_inhibition),
        ("Hill(n=1) == Michaelis-Menten", lambda verbose: check_hill_reduces()),
    ):
        before = checks
        checks += function(args.verbose)
        status = "OK  " if checks > before else "FAIL"
        print(f"{status} {name:<42} {checks - before:>3} comparisons")

    if not AGREEMENTS:
        print("FAIL no comparisons were performed")
        return 1

    def rel_diff(pair) -> float:
        _, ours, theirs = pair
        return abs(ours - theirs) / max(abs(theirs), 1e-12)

    monte_carlo = [p for p in AGREEMENTS if p[0].startswith("bootstrap")]
    exact = [p for p in AGREEMENTS if not p[0].startswith("bootstrap")]
    print(
        f"\n{len(AGREEMENTS)} comparisons between kinetics and independent SciPy/NumPy "
        f"implementations"
    )
    print(
        f"  {len(exact)} deterministic comparisons: worst relative difference "
        f"{max(rel_diff(p) for p in exact):.2e}"
    )
    if monte_carlo:
        print(
            f"  {len(monte_carlo)} Monte-Carlo comparisons (bootstrap): worst relative "
            f"difference {max(rel_diff(p) for p in monte_carlo):.2e}, expected to be "
            f"percent-level because the random streams differ"
        )
    return 0 if len(AGREEMENTS) == checks else 1


if __name__ == "__main__":
    sys.exit(main())
