"""Calibration lines and detection limits.

A calibration is the one place in this package where the answer is usually a
regression slope, so the dangerous failure is not a crash — it is a tidy R² that
hides curvature, or a LOD computed with σ from a different dataset.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from kinetics.calibration import CalibrationResult, fit_calibration, lod_loq

SLOPE, INTERCEPT = 0.5, 0.01


def test_exact_standards_recover_the_line():
    x = np.array([0.0, 1.0, 2.0, 5.0, 10.0, 20.0])
    y = SLOPE * x + INTERCEPT
    res = fit_calibration(x, y)
    assert res.slope == pytest.approx(SLOPE, rel=1e-9)
    assert res.intercept == pytest.approx(INTERCEPT, rel=1e-9)
    assert res.r_squared == pytest.approx(1.0)
    assert res.residual_sd == pytest.approx(0.0, abs=1e-12)


def test_noiseless_calibration_has_a_numerically_zero_detection_limit():
    """σ = 0 gives LOD ≈ 0 — arithmetic, not a measurement claim.

    The residual σ on perfectly linear points is not exactly zero in floating point
    (measured here: about 1.7e-16), so LOD comes out around 6e-15 rather than 0. A
    detection limit of 1e-14 µM is a statement about double precision, which is why
    the docstring of ``lod_loq`` says the σ must come from real replicates.
    """
    x = np.array([0.0, 1.0, 2.0, 5.0, 10.0])
    res = fit_calibration(x, SLOPE * x + INTERCEPT)
    assert 0.0 <= res.lod < 1e-9
    assert 0.0 <= res.loq < 1e-8
    assert res.loq > res.lod


def test_slope_and_intercept_match_polyfit(calibration_data):
    x, y = calibration_data()
    res = fit_calibration(x, y)
    slope, intercept = np.polyfit(x, y, 1)
    assert res.slope == pytest.approx(slope)
    assert res.intercept == pytest.approx(intercept)


def test_standard_errors_match_the_textbook_formulas(calibration_data):
    x, y = calibration_data(noise=0.05, seed=3)
    res = fit_calibration(x, y)
    residuals = y - (res.slope * x + res.intercept)
    dof = len(x) - 2
    sigma = math.sqrt(float(np.sum(residuals**2)) / dof)
    sxx = float(np.sum((x - x.mean()) ** 2))
    assert res.residual_sd == pytest.approx(sigma)
    assert res.stderr_slope == pytest.approx(sigma / math.sqrt(sxx))
    assert res.stderr_intercept == pytest.approx(
        sigma * math.sqrt(1.0 / len(x) + x.mean() ** 2 / sxx)
    )


def test_lod_and_loq_are_the_iupac_multiples_of_sigma_over_slope(calibration_data):
    x, y = calibration_data()
    res = fit_calibration(x, y)
    assert res.lod == pytest.approx(3.3 * res.residual_sd / res.slope)
    assert res.loq == pytest.approx(10.0 * res.residual_sd / res.slope)
    assert res.loq > res.lod


def test_loq_is_three_times_lod_with_the_default_multipliers():
    lod, loq = lod_loq(slope=0.25, residual_sd=0.02)
    assert lod == pytest.approx(3.3 * 0.02 / 0.25)
    assert loq == pytest.approx(10.0 * 0.02 / 0.25)


def test_alternative_conventions_are_reachable_and_explicit():
    """Some fields use k = 3 and 10; changing it is a stated decision, so it is a knob."""
    lod, loq = lod_loq(slope=0.25, residual_sd=0.02, k_lod=3.0, k_loq=6.0)
    assert lod == pytest.approx(3.0 * 0.02 / 0.25)
    assert loq == pytest.approx(6.0 * 0.02 / 0.25)


def test_zero_slope_has_no_detection_limit():
    with pytest.raises(ValueError, match="zero-slope"):
        lod_loq(slope=0.0, residual_sd=0.01)


def test_negative_sigma_is_refused():
    with pytest.raises(ValueError, match="cannot be negative"):
        lod_loq(slope=1.0, residual_sd=-0.01)


def test_two_standards_are_refused_because_sigma_is_undefined():
    with pytest.raises(ValueError, match="at least three standards"):
        fit_calibration([1.0, 2.0], [0.5, 1.0])


def test_mismatched_lengths_are_refused():
    with pytest.raises(ValueError, match="same length"):
        fit_calibration([1.0, 2.0, 3.0], [1.0, 2.0])


def test_a_single_concentration_cannot_have_a_slope():
    with pytest.raises(ValueError, match="same concentration"):
        fit_calibration([5.0, 5.0, 5.0], [1.0, 1.1, 0.9])


def test_nan_in_a_standard_is_refused_not_dropped(calibration_data):
    x, y = calibration_data()
    y = np.array(y, dtype=float)
    y[2] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        fit_calibration(x, y)


def test_high_r_squared_comes_with_a_warning_because_it_is_not_linearity(calibration_data):
    x, y = calibration_data(noise=0.01)
    res = fit_calibration(x, y)
    assert res.r_squared > 0.995
    assert any("not evidence of linearity" in note for note in res.notes)


def test_a_non_zero_blank_is_flagged_when_it_is_large_enough(calibration_data):
    x, y = calibration_data(intercept=2.0, noise=0.02)
    res = fit_calibration(x, y)
    assert any(
        "intercept is more than three standard errors from zero" in note for note in res.notes
    )


def test_proportional_weighting_changes_the_fit_when_error_grows_with_signal():
    """Heteroscedastic standards: the weighted fit must follow the low standards."""
    x = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    y = np.array([0.26, 0.49, 1.02, 2.44, 5.20, 9.60])  # low standards read high
    sd = np.array([0.01, 0.02, 0.04, 0.10, 0.20, 0.40])
    unweighted = fit_calibration(x, y)
    weighted = fit_calibration(x, y, weighted=True, replicate_sd=sd)
    assert weighted.weighted is True
    assert abs(weighted.slope - 0.5) < abs(unweighted.slope - 0.5)
    assert any("replicate standard deviations" in note for note in weighted.notes)


def test_proportional_rule_requires_positive_signals():
    with pytest.raises(ValueError, match="positive signals"):
        fit_calibration([0.0, 1.0, 2.0], [0.0, 1.0, 2.0], weighted=True)


def test_weighted_r_squared_carries_its_caveat(calibration_data):
    x, y = calibration_data()
    res = fit_calibration(x, y, weighted=True)
    assert any("unweighted residuals" in note for note in res.notes)


def test_inverting_the_calibration_round_trips(calibration_data):
    x, y = calibration_data()
    res = fit_calibration(x, y)
    signal = res.slope * 7.5 + res.intercept
    assert res.concentration(signal) == pytest.approx(7.5)


def test_result_serialises_with_none_for_non_finite_limits():
    res = CalibrationResult(
        slope=1.0,
        intercept=0.0,
        stderr_slope=0.1,
        stderr_intercept=0.1,
        r_squared=0.99,
        residual_sd=0.0,
        n=5,
    )
    payload = res.as_dict()
    assert payload["lod"] == 0.0 and payload["loq"] == 0.0

    res_bad = CalibrationResult(
        slope=float("nan"),
        intercept=0.0,
        stderr_slope=0.1,
        stderr_intercept=0.1,
        r_squared=0.99,
        residual_sd=float("nan"),
        n=5,
    )
    assert res_bad.as_dict()["lod"] is None


def test_weighted_residual_sd_is_in_signal_units_not_a_scale_factor():
    """Regression test for a units bug that produced an LOD of 66 µM.

    ``sqrt(weighted_rss/dof)`` with weights ``1/sd²`` is dimensionless — it came out
    near 1.0 and was then divided by the slope, inflating the LOD by two orders of
    magnitude. The reported σ must be on the scale of the replicate SDs.
    """
    x = np.array([0.0, 5.0, 10.0, 25.0, 50.0, 100.0])
    sd = np.array([0.0025, 0.003, 0.006, 0.015, 0.03, 0.06])
    y = 0.05 * x + 0.01
    rng = np.random.default_rng(4)
    y = y + sd * rng.standard_normal(x.size)
    res = fit_calibration(x, y, replicate_sd=sd)
    assert 0.001 < res.residual_sd < 0.01
    assert res.lod < 1.0  # concentration units, not tens of µM
    assert any("scale factor" in note for note in res.notes)


def test_weighting_matters_by_an_order_of_magnitude_on_heteroscedastic_standards():
    """The repository's calibration data: 4 % proportional error over 0–200 µM.

    Measured: unweighted slope 0.047635 with LOD 6.36 µM; weighted slope 0.049013
    with LOD 0.109 µM. An unweighted fit is not "the safe default" here — it reports a
    detection limit 58× too high, because the residuals are dominated by the strongest
    standards, which are the ones with the largest noise.
    """
    import csv
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "datasets/sample/calibration_mean.csv"
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    x = [float(r["concentration"]) for r in rows]
    y = [float(r["signal"]) for r in rows]
    sd = [float(r["sd"]) for r in rows]

    unweighted = fit_calibration(x, y)
    weighted = fit_calibration(x, y, weighted=True, replicate_sd=sd)
    assert weighted.lod < unweighted.lod / 10
    assert abs(weighted.slope - 0.0485) < abs(unweighted.slope - 0.0485)
