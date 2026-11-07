"""Error paths and edge cases that would otherwise be untested code.

Each test here corresponds to a guard in the package. They are collected in one file
so that the guard list is readable in one place: every message a user can see is
asserted somewhere, which is also how the docs stay honest about what the tool
refuses to do.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from kinetics.calibration import CalibrationResult, fit_calibration
from kinetics.fitting import (
    aic,
    bootstrap_ci,
    fit_eadie_hofstee,
    fit_hanes,
    fit_michaelis_menten,
    fit_nonlinear,
)
from kinetics.models import hill, michaelis_menten, mixed_inhibition
from kinetics.report import results_to_markdown


def test_aic_needs_more_points_than_parameters():
    with pytest.raises(ValueError, match="more data points than parameters"):
        aic(1.0, 2, 2)


def test_fit_nonlinear_checks_lengths_and_finiteness():
    with pytest.raises(ValueError, match="same length"):
        fit_nonlinear([1.0, 2.0, 3.0], [1.0, 2.0], michaelis_menten, (10.0, 1.0), ("vmax", "km"))
    with pytest.raises(ValueError, match="must be finite"):
        fit_nonlinear(
            [1.0, 2.0, float("nan")], [1.0, 2.0, 3.0], michaelis_menten, (10.0, 1.0), ("vmax", "km")
        )
    with pytest.raises(ValueError, match="need more points"):
        fit_nonlinear([1.0, 2.0], [1.0, 2.0], michaelis_menten, (10.0, 1.0), ("vmax", "km"))


def test_fit_nonlinear_checks_sigma_length():
    with pytest.raises(ValueError, match="one value per data point"):
        fit_nonlinear(
            [1.0, 2.0, 5.0, 10.0],
            [1.0, 1.8, 3.4, 5.0],
            michaelis_menten,
            (10.0, 25.0),
            ("vmax", "km"),
            sigma=[0.1, 0.1],
        )


def test_linear_fit_needs_three_points():
    with pytest.raises(ValueError, match="fewer than three points"):
        fit_hanes([1.0, 2.0], [1.0, 2.0])


def test_a_line_through_every_point_has_no_residual_variance():
    """Exactly linear data: the slope's standard error is undefined, and says so.

    This calls the private helper because no public transformation produces exactly
    collinear points from real kinetics — Eadie-Hofstee on exact Michaelis-Menten data
    leaves a residual variance of order 1e-27, which yields a standard error of 5e-14
    rather than a NaN. That is also the honest answer: a standard error of 5e-14 is
    "this fit has no noise to estimate from", not a precise measurement.
    """
    from kinetics.fitting import _linear_fit

    res = _linear_fit([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0], "exact", [])
    assert math.isnan(res.stderr["slope"])
    assert any("floating-point noise" in note for note in res.notes)

    S = np.array([1.0, 2.0, 5.0, 10.0])
    transformed = fit_eadie_hofstee(S, michaelis_menten(S, 10.0, 25.0))
    assert math.isnan(transformed.stderr["slope"])

    noisy = fit_eadie_hofstee(S, michaelis_menten(S, 10.0, 25.0) + np.array([0.1, -0.1, 0.1, -0.1]))
    assert math.isfinite(noisy.stderr["slope"])


def test_eadie_hofstee_refuses_zero_substrate():
    with pytest.raises(ValueError, match="divides by .S."):
        fit_eadie_hofstee([0.0, 1.0, 2.0], [0.0, 1.0, 2.0])


def test_bootstrap_reports_when_every_resample_fails():
    """Two-parameter fit on resamples of a perfectly flat response cannot converge."""
    with pytest.raises(ValueError, match="every bootstrap resample failed"):
        bootstrap_ci([1.0, 2.0, 3.0, 4.0], [5.0, 5.0, 5.0, 5.0], n_resamples=5, seed=1)


def test_hill_refuses_non_positive_k_half_and_negative_substrate():
    with pytest.raises(ValueError, match="k_half must be positive"):
        hill(10.0, 10.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="cannot be negative"):
        hill(-1.0, 10.0, 25.0, 1.0)


def test_mixed_inhibition_refuses_non_positive_alpha():
    with pytest.raises(ValueError, match="alpha must be positive"):
        mixed_inhibition(10.0, 10.0, 25.0, 4.0, 8.0, alpha=0.0)


def test_calibration_cannot_be_inverted_without_a_slope():
    res = CalibrationResult(
        slope=0.0,
        intercept=1.0,
        stderr_slope=0.1,
        stderr_intercept=0.1,
        r_squared=0.0,
        residual_sd=0.1,
        n=4,
    )
    with pytest.raises(ValueError, match="zero slope"):
        res.concentration(2.5)


def test_replicate_sd_must_be_positive_and_the_right_length():
    x = [1.0, 2.0, 3.0, 4.0]
    y = [1.0, 2.0, 3.0, 4.0]
    with pytest.raises(ValueError, match="positive and match"):
        fit_calibration(x, y, replicate_sd=[0.1, 0.1])
    with pytest.raises(ValueError, match="positive and match"):
        fit_calibration(x, y, replicate_sd=[0.1, 0.0, 0.1, 0.1])


def test_report_renders_a_non_numeric_parameter_without_crashing():
    from kinetics.fitting import FitResult

    odd = FitResult(
        model="m", params={"label": "not-a-number"}, n=5, k=1, rss=1.0, extra={"tss": 2.0}
    )
    text = results_to_markdown([odd])
    assert "not-a-number" in text


def test_reports_deduplicate_repeated_notes():
    from kinetics.fitting import FitResult

    same = "identical note"
    fits = [
        FitResult(model="a", params={"x": 1.0}, n=5, k=1, rss=1.0, notes=[same]),
        FitResult(model="b", params={"x": 2.0}, n=5, k=1, rss=1.0, notes=[same]),
    ]
    assert results_to_markdown(fits).count(same) == 1


def test_fit_of_a_single_parameter_model_is_possible():
    """fit_nonlinear is general: the tests use it for MM, the CLI for Hill."""
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = 3.0 * x
    res = fit_nonlinear(x, y, lambda xx, a: a * xx, (1.0,), ("a",), name="straight")
    assert res.params["a"] == pytest.approx(3.0, rel=1e-6)
    assert res.model == "straight"


def test_michaelis_menten_fit_rejects_mismatched_lengths(mm_data):
    S, v = mm_data()
    with pytest.raises(ValueError, match="same length"):
        fit_michaelis_menten(S, v[:-1])


def test_a_constant_response_is_refused_rather_than_fitted():
    """No optimiser can identify a model from a flat response, but some will pretend."""
    with pytest.raises(ValueError, match="response is constant"):
        fit_michaelis_menten([1.0, 2.0, 3.0, 4.0], [5.0, 5.0, 5.0, 5.0])


def test_bootstrap_on_a_flat_response_refuses_at_any_resample_count():
    """Regression test: with few resamples it refused, with many it returned (10, 10)."""
    for n_resamples in (5, 100):
        with pytest.raises(ValueError, match="every bootstrap resample failed"):
            bootstrap_ci(
                [1.0, 2.0, 3.0, 4.0], [5.0, 5.0, 5.0, 5.0], n_resamples=n_resamples, seed=1
            )


def test_an_unestimable_covariance_becomes_nan_and_a_note():
    """SciPy's OptimizeWarning used to be printed and then lost."""
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = np.array([2.0, 3.0, 2.5, 3.5, 3.0])
    result = fit_nonlinear(
        x, y, lambda t, a, b: a * 0 + t * 0 + b, (1.0, 1.0), ("a", "b"), name="degenerate"
    )
    assert math.isnan(result.stderr["a"]) and math.isnan(result.stderr["b"])
    assert any("could not estimate the covariance" in note for note in result.notes)


def test_non_positive_rates_are_fitted_but_flagged():
    """The fit does not refuse a negative rate; it reports it and leaves the call to you.

    Refusing was the first design here and it was wrong: at the bottom of a substrate
    series a negative rate is the detector floor, and dropping such points biases Vmax
    upwards. This package's own cooperative_hill dataset carries one, so a hard guard
    would have made the package unable to fit its own data.
    """
    S = [1.0, 2.5, 5.0, 10.0, 20.0, 40.0]
    v = [1.05, 1.92, 3.4, -0.4, 6.61, 7.98]
    result = fit_michaelis_menten(S, v)
    assert result.params["vmax"] > 0
    flagged = [note for note in result.notes if "are not positive" in note]
    assert flagged and "lowest -0.4" in flagged[0] and "1 of 6" in flagged[0]
    # A point exactly at zero is counted too, and a clean series carries no such note.
    zero_v = [1.05, 1.92, 3.4, 0.0, 6.61, 7.98]
    assert any("are not positive" in note for note in fit_michaelis_menten(S, zero_v).notes)
    clean = [1.05, 1.92, 3.4, 5.2, 6.61, 7.98]
    assert not any("are not positive" in note for note in fit_michaelis_menten(S, clean).notes)
