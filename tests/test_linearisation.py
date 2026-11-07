"""The bias of the linearisations, measured rather than asserted.

Every claim in this file is a number produced by the loops in it. The reference
"truth" is Vmax = 10.0 and Km = 25.0 µM, and the noise is constant-variance
Gaussian at σ = 0.1, which is 1 % of Vmax — a good assay.

Measured over 300 seeds (the numbers the assertions are written against):

=================  ============  ==========  ==================
method             med |ΔVmax|  med |ΔKm|  within 20 % Vmax
=================  ============  ==========  ==================
nonlinear                 0.121        0.840               1.00
hanes                     0.206        1.695               1.00
eadie_hofstee             0.530        2.812               0.97
lineweaver_burk           3.980       13.130               0.22
=================  ============  ==========  ==================

The Lineweaver-Burk estimate is wrong by more than the parameter itself in the
tails (IQR of Vmax runs from 3.0 to 13.9), which is not a mistake in the code: it
is what weighting the least reliable points most heavily does.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from kinetics.fitting import (
    compare_linearisations,
    fit_eadie_hofstee,
    fit_hanes,
    fit_lineweaver_burk,
)
from kinetics.models import michaelis_menten

S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
VMAX_TRUE, KM_TRUE = 10.0, 25.0


def _estimates(noise: float = 0.1, seeds: int = 300) -> dict[str, list[tuple[float, float]]]:
    """Run every method over many seeded datasets and collect (Vmax, Km)."""
    out: dict[str, list[tuple[float, float]]] = {}
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        v = michaelis_menten(S, VMAX_TRUE, KM_TRUE) + noise * rng.standard_normal(S.size)
        if np.any(v <= 0):
            # Hanes-Woolf divides by v; a non-positive simulated rate is not a
            # transformation problem, so it is excluded here and tested separately.
            continue
        for result in compare_linearisations(S, v):
            out.setdefault(result.model, []).append((result.params["vmax"], result.params["km"]))
    return out


@pytest.fixture(scope="module")
def estimates():
    return _estimates()


def test_every_method_recovers_exact_parameters(estimates):
    """Sanity check on the estimators themselves: noiseless data has one answer."""
    v = michaelis_menten(S, VMAX_TRUE, KM_TRUE)
    for result in compare_linearisations(S, v):
        assert result.params["vmax"] == pytest.approx(VMAX_TRUE, rel=1e-6), result.model
        assert result.params["km"] == pytest.approx(KM_TRUE, rel=1e-3), result.model


def test_order_is_fixed_so_reports_do_not_reshuffle(estimates):
    v = michaelis_menten(S, VMAX_TRUE, KM_TRUE) + 0.1
    assert [r.model for r in compare_linearisations(S, v)] == [
        "michaelis_menten",
        "hanes",
        "eadie_hofstee",
        "lineweaver_burk",
    ]


def test_nonlinear_fit_is_the_least_biased_estimator(estimates):
    def median_error(model: str, index: int, truth: float) -> float:
        values = np.array([pair[index] for pair in estimates[model]])
        return float(np.median(np.abs(values - truth)))

    mm = median_error("michaelis_menten", 0, VMAX_TRUE)
    hanes = median_error("hanes", 0, VMAX_TRUE)
    eh = median_error("eadie_hofstee", 0, VMAX_TRUE)
    lb = median_error("lineweaver_burk", 0, VMAX_TRUE)
    assert mm < hanes < eh < lb
    assert lb > 5 * mm


def test_lineweaver_burk_is_unreliable_not_merely_less_precise(estimates):
    """Fewer than half of its Vmax estimates land within 20 % of the truth."""
    values = np.array([pair[0] for pair in estimates["lineweaver_burk"]])
    within = np.mean(np.abs(values - VMAX_TRUE) / VMAX_TRUE < 0.2)
    assert within < 0.5

    mm_values = np.array([pair[0] for pair in estimates["michaelis_menten"]])
    assert np.all(np.abs(mm_values - VMAX_TRUE) / VMAX_TRUE < 0.2)


def test_lineweaver_burk_fails_catastrophically_and_more_often_with_more_noise():
    """Fraction of Vmax estimates wrong by more than 100 % of the true value.

    Measured over 200 seeds: Michaelis-Menten 0.00 at every noise level; Eadie-Hofstee
    0.00; Hanes-Woolf 0.00 / 0.00 / 0.03; Lineweaver-Burk 0.12 / 0.30 / 0.42 for
    σ = 0.05 / 0.10 / 0.20. The nonlinear fit does not merely do better on average —
    it does not produce these failures at all.
    """

    def failure_fraction(noise: float, model: str, seeds: int = 200) -> float:
        failures = 0
        total = 0
        for seed in range(seeds):
            rng = np.random.default_rng(seed)
            v = michaelis_menten(S, VMAX_TRUE, KM_TRUE) + noise * rng.standard_normal(S.size)
            if np.any(v <= 0):
                continue
            for result in compare_linearisations(S, v):
                if result.model == model:
                    total += 1
                    failures += abs(result.params["vmax"] - VMAX_TRUE) / VMAX_TRUE > 1.0
        return failures / total

    assert failure_fraction(0.05, "lineweaver_burk") < failure_fraction(0.20, "lineweaver_burk")
    assert failure_fraction(0.10, "lineweaver_burk") > 0.15
    assert failure_fraction(0.20, "michaelis_menten") == 0.0


def test_linearised_results_carry_their_bias_note():
    v = michaelis_menten(S, VMAX_TRUE, KM_TRUE) + 0.1
    for result in compare_linearisations(S, v):
        assert result.method == "linearised" if result.model != "michaelis_menten" else True
        if result.model != "michaelis_menten":
            assert any("reciprocal" in n or "v/S" in n or "S/v" in n for n in result.notes), (
                result.model
            )


def test_proportional_data_has_no_finite_vmax_and_the_fit_says_so():
    """1/v against 1/S through the origin: no intercept, so no Vmax."""
    result = fit_lineweaver_burk(S, 2.0 * S)
    assert math.isinf(result.params["vmax"])
    assert math.isnan(result.params["km"])
    assert any("intercept is zero" in note for note in result.notes)
    assert "vmax" in result.params  # keys kept: no KeyError downstream


def test_non_positive_rates_are_refused_by_the_transformations():
    with pytest.raises(ValueError, match="positive rates"):
        fit_lineweaver_burk([1.0, 2.0, 5.0], [1.0, 0.0, -1.0])
    with pytest.raises(ValueError, match="all rates must be positive"):
        fit_hanes([1.0, 2.0, 5.0], [1.0, 0.0, 3.0])


def test_identical_x_values_are_refused_not_divided_by():
    """A zero spread in x is an unidentifiable slope: say so instead of returning inf."""
    with pytest.raises(ValueError, match="every x value is identical"):
        fit_hanes([5.0, 5.0, 5.0, 5.0], [1.0, 2.0, 3.0, 4.0])


def test_eadie_hofstee_puts_v_on_both_axes_and_its_note_says_why_that_matters():
    v = michaelis_menten(S, VMAX_TRUE, KM_TRUE) + 0.1
    result = fit_eadie_hofstee(S, v)
    assert any("both axes" in note for note in result.notes)


def test_hanes_reports_the_least_biased_linearisation_on_these_data(estimates):
    def med(model: str) -> float:
        return float(np.median([abs(p[0] - VMAX_TRUE) for p in estimates[model]]))

    assert med("hanes") < med("eadie_hofstee")
