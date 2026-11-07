"""Nonlinear fitting, information criteria and the FitResult contract."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from kinetics.fitting import (
    FitResult,
    aic,
    aicc,
    bootstrap_ci,
    fit_michaelis_menten,
    fit_nonlinear,
)
from kinetics.models import hill, michaelis_menten


def test_aic_matches_its_definition():
    rss, n, k = 4.0, 10, 2
    assert aic(rss, n, k) == pytest.approx(n * math.log(rss / n) + 2 * k)


def test_aicc_adds_the_small_sample_penalty():
    rss, n, k = 4.0, 10, 2
    expected = aic(rss, n, k) + 2 * k * (k + 1) / (n - k - 1)
    assert aicc(rss, n, k) == pytest.approx(expected)


def test_aicc_is_undefined_when_no_degrees_of_freedom_remain():
    with pytest.raises(ValueError, match="AICc undefined"):
        aicc(1.0, 3, 2)


def test_perfect_fit_gives_minus_infinity_instead_of_a_number():
    assert aic(0.0, 10, 2) == float("-inf")


def test_aic_refuses_negative_rss():
    with pytest.raises(ValueError, match="rss cannot be negative"):
        aic(-1.0, 10, 2)


def test_more_points_costs_two_aic_units_per_parameter():
    assert aic(4.0, 10, 3) - aic(4.0, 10, 2) == pytest.approx(2.0)


def test_nonlinear_fit_returns_the_parameters_it_was_given(mm_data):
    """The only honest way to check a fitter: give it the answer."""
    S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    v = michaelis_menten(S, 10.0, 25.0)
    res = fit_michaelis_menten(S, v)
    assert res.params["vmax"] == pytest.approx(10.0, rel=1e-6)
    assert res.params["km"] == pytest.approx(25.0, rel=1e-6)
    assert res.rss == pytest.approx(0.0, abs=1e-12)


def test_fit_records_its_method_and_shape(mm_data):
    S, v = mm_data()
    res = fit_michaelis_menten(S, v)
    assert res.method == "nonlinear"
    assert (res.n, res.k) == (len(S), 2)
    assert res.model == "michaelis_menten"


def test_estimates_are_near_the_truth_under_moderate_noise(mm_data):
    S, v = mm_data(noise=0.20, seed=7)
    res = fit_michaelis_menten(S, v)
    assert res.params["vmax"] == pytest.approx(10.0, rel=0.10)
    assert res.params["km"] == pytest.approx(25.0, rel=0.20)


def test_weighting_changes_the_answer_when_one_point_is_unreliable(mm_data):
    """Down-weighting a wild point must pull the fit towards the others."""
    S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    v = michaelis_menten(S, 10.0, 25.0)
    v[0] += 4.0  # the lowest-[S] point reads far too high
    unweighted = fit_michaelis_menten(S, v)
    sigma = np.ones_like(S)
    sigma[0] = 3.0  # ... and it is the least trustworthy point
    weighted = fit_michaelis_menten(S, v, sigma=sigma)
    assert abs(weighted.params["vmax"] - 10.0) <= abs(unweighted.params["vmax"] - 10.0)


def test_sigma_must_be_positive_and_match(mm_data):
    S, v = mm_data()
    with pytest.raises(ValueError):
        fit_michaelis_menten(S, v, sigma=[1.0] * (len(S) - 1))


def test_nan_in_the_data_is_an_error_not_a_silent_drop(mm_data):
    """A fit that quietly discards rows cannot tell you the row was bad."""
    S, v = mm_data()
    S_bad = np.array(S, dtype=float)
    v_bad = np.array(v, dtype=float)
    v_bad[2] = np.nan
    with pytest.raises(ValueError):
        fit_michaelis_menten(S_bad, v_bad)


def test_fit_nonlinear_counts_starting_values(mm_data):
    S, v = mm_data()
    with pytest.raises(ValueError, match="each parameter needs exactly one start"):
        fit_nonlinear(S, v, hill, (10.0,), ("vmax", "k_half", "n"), name="hill")


def test_fit_result_serialises_to_json(mm_data):
    S, v = mm_data()
    payload = json.dumps(fit_michaelis_menten(S, v).as_dict())
    back = json.loads(payload)
    assert back["model"] == "michaelis_menten"
    assert isinstance(back["aicc"], float)
    assert set(back["params"]) == {"vmax", "km"}


def test_r_squared_is_nan_when_the_response_has_no_variance():
    fit = FitResult(model="m", params={"a": 1.0}, n=5, k=1, rss=1.0, extra={"tss": 0.0})
    assert math.isnan(fit.r_squared)


def test_aicc_or_none_swallows_only_the_undefined_case():
    small = FitResult(model="m", params={"a": 1.0}, n=3, k=2, rss=1.0, extra={"tss": 2.0})
    assert small.aicc_or_none() is None
    with pytest.raises(ValueError):
        _ = small.aicc
    fine = FitResult(model="m", params={"a": 1.0}, n=10, k=2, rss=1.0, extra={"tss": 2.0})
    assert fine.aicc_or_none() == pytest.approx(aicc(1.0, 10, 2))


def test_hill_prefers_the_right_model_and_pays_for_the_wrong_one():
    """AICc penalises the extra parameter, but only when it buys nothing.

    Counted over 20 seeds, because a single seed proves nothing about a model
    selection rule. Measured on this machine: Hill wins 2 of 20 seeds on
    Michaelis-Menten data and 20 of 20 on genuinely sigmoid (n = 2) data. The
    bounds below are looser than the measurement so the test is not re-measuring
    floating point, but tight enough to catch a broken penalty term.
    """
    S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])

    def hill_wins(shape: float, noise: float) -> int:
        wins = 0
        for seed in range(20):
            rng = np.random.default_rng(seed)
            v = 10.0 * S**shape / (25.0**shape + S**shape)
            v = v + noise * rng.standard_normal(S.size)
            mm = fit_michaelis_menten(S, v)
            hi = fit_nonlinear(S, v, hill, (10.0, 25.0, 1.0), ("vmax", "k_half", "n"), name="hill")
            wins += hi.aicc < mm.aicc
        return wins

    assert hill_wins(1.0, 0.10) <= 5
    assert hill_wins(2.0, 0.10) >= 18


def test_bootstrap_interval_contains_the_point_estimate(mm_data):
    S, v = mm_data(noise=0.3, seed=2)
    point = fit_michaelis_menten(S, v).params
    ci = bootstrap_ci(S, v, n_resamples=300, seed=1)
    for name in ("vmax", "km"):
        assert ci[name][0] <= point[name] <= ci[name][1]


def test_bootstrap_is_reproducible_from_its_seed(mm_data):
    S, v = mm_data(noise=0.3, seed=2)
    assert bootstrap_ci(S, v, n_resamples=200, seed=42) == bootstrap_ci(
        S, v, n_resamples=200, seed=42
    )


def test_bootstrap_changes_with_the_seed_but_not_much(mm_data):
    S, v = mm_data(noise=0.3, seed=2)
    a = bootstrap_ci(S, v, n_resamples=200, seed=1)
    b = bootstrap_ci(S, v, n_resamples=200, seed=2)
    assert a != b
    assert a["vmax"][0] == pytest.approx(b["vmax"][0], rel=0.1)


def test_bootstrap_refuses_too_few_points():
    with pytest.raises(ValueError, match="fewer than four points"):
        bootstrap_ci([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], n_resamples=10, seed=1)


def test_bootstrap_confidence_level_widens_the_interval(mm_data):
    S, v = mm_data(noise=0.4, seed=9)
    narrow = bootstrap_ci(S, v, n_resamples=400, seed=3, confidence=0.50)
    wide = bootstrap_ci(S, v, n_resamples=400, seed=3, confidence=0.99)
    assert (wide["vmax"][1] - wide["vmax"][0]) > (narrow["vmax"][1] - narrow["vmax"][0])
