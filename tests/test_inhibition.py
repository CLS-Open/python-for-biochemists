"""Inhibition typing by global fitting — and the bug that made it lie.

The regression test at the top of this file exists because passing the inhibitor
vector through ``extra_args`` bound it to ``ki`` and ``ki`` to ``I``. The fit still
"converged" for some data, which is worse than crashing: the numbers were wrong.
Any future change to how the (S, I, v) models are called must keep the forward model
and the fitted model in agreement, and ``test_forward_and_fitted_models_agree`` is
what checks it.
"""

from __future__ import annotations

import numpy as np
import pytest

from kinetics.inhibition import compare_inhibition_models, fit_inhibition_series
from kinetics.models import INHIBITION_MODELS, competitive

VMAX, KM, KI = 10.0, 25.0, 8.0


def test_forward_and_fitted_models_agree():
    """Noiseless data from a model must return that model's parameters.

    This is the direct check on argument order: if ``I`` and ``ki`` are swapped
    anywhere in the call chain, the recovered ``ki`` is not 8.0.
    """
    for mechanism in ("competitive", "noncompetitive", "uncompetitive"):
        func = INHIBITION_MODELS[mechanism]["func"]
        rows = [
            {"S": S, "I": I, "v": float(func(S, VMAX, KM, I, KI))}
            for S in (2.0, 5.0, 10.0, 25.0, 100.0)
            for I in (0.0, 4.0, 16.0)
        ]
        res = fit_inhibition_series(rows, mechanism)
        assert res.params["ki"] == pytest.approx(KI, rel=1e-4), mechanism
        assert res.params["vmax"] == pytest.approx(VMAX, rel=1e-4), mechanism
        assert res.params["km"] == pytest.approx(KM, rel=1e-4), mechanism


def test_inhibitor_actually_changes_the_rate_in_the_fitted_data(inhibition_rows):
    """Guard against a fit that is 'right' only because I never entered the model."""
    rows = inhibition_rows(mechanism="competitive", noise=0.0)
    at_zero = [r["v"] for r in rows if r["I"] == 0.0 and r["S"] == 10.0][0]
    at_high = [r["v"] for r in rows if r["I"] == 16.0 and r["S"] == 10.0][0]
    assert at_high < at_zero
    assert at_zero == pytest.approx(competitive(10.0, VMAX, KM, 0.0, KI))


@pytest.mark.parametrize("mechanism", ["competitive", "noncompetitive", "uncompetitive"])
def test_the_mechanism_that_generated_the_data_wins_on_a_clean_series(mechanism):
    func = INHIBITION_MODELS[mechanism]["func"]
    rows = [
        {"S": S, "I": I, "v": float(func(S, VMAX, KM, I, KI))}
        for S in (1.0, 2.0, 5.0, 10.0, 25.0, 100.0)
        for I in (0.0, 4.0, 16.0, 64.0)
    ]
    cmp = compare_inhibition_models(rows)
    assert cmp.best == mechanism
    assert cmp.delta_aicc[mechanism] == 0.0


def test_moderate_noise_keeps_the_right_answer(inhibition_rows):
    rows = inhibition_rows(mechanism="competitive", noise=0.10, seed=5)
    assert compare_inhibition_models(rows).best == "competitive"


def test_results_are_sorted_best_first_and_delta_is_zero_for_the_winner(inhibition_rows):
    cmp = compare_inhibition_models(inhibition_rows(noise=0.2))
    assert cmp.results[0].model == cmp.best
    assert cmp.delta_aicc[cmp.best] == 0.0
    assert [r.aicc for r in cmp.results] == sorted(r.aicc for r in cmp.results)
    assert all(d >= 0 for d in cmp.delta_aicc.values())


def test_a_series_without_inhibitor_is_refused_not_fitted():
    rows = [
        {"S": S, "I": 0.0, "v": 10.0 * S / (25.0 + S)} for S in (1.0, 2.0, 5.0, 10.0, 25.0, 100.0)
    ]
    with pytest.raises(ValueError, match="at least two inhibitor concentrations"):
        compare_inhibition_models(rows)


def test_weak_inhibition_cannot_be_typed_and_the_report_says_so():
    """The case that matters: an experiment that cannot answer the question.

    Simulated data from the *competitive* mechanism with Ki = 500 µM, but the highest
    inhibitor concentration used is 2 µM — so [I]/Ki never exceeds 0.004 and all three
    mechanisms describe the data equally well. Measured on these rows: ΔAICc is 0.00
    (uncompetitive), 0.79 (noncompetitive) and 1.19 (competitive). The nominal winner
    is not the mechanism that generated the data, and every difference is inside the
    "indistinguishable" band, which is exactly what the notes have to say.
    """
    rng = np.random.default_rng(1)
    rows = []
    for S in (1.0, 2.0, 5.0, 10.0, 25.0, 100.0):
        for I in (0.0, 1.0, 2.0):
            v = competitive(S, VMAX, KM, I, 500.0) + 0.05 * rng.standard_normal()
            rows.append({"S": S, "I": I, "v": float(v)})

    cmp = compare_inhibition_models(rows)
    assert cmp.best != "competitive"  # the truth is not the winner
    assert cmp.delta_aicc["competitive"] < 2.0
    ambiguity = [n for n in cmp.notes if "indistinguishable" in n]
    assert ambiguity and "noncompetitive" in ambiguity[0]


def test_strong_inhibition_separates_the_mechanisms(inhibition_rows):
    """The contrast: with [I] up to 16 µM and Ki = 8 µM, the differences are decisive."""
    cmp = compare_inhibition_models(inhibition_rows(mechanism="competitive", noise=0.05))
    assert cmp.best == "competitive"
    assert min(v for k, v in cmp.delta_aicc.items() if k != "competitive") > 2.0


def test_every_comparison_carries_the_cannot_prove_a_mechanism_note(inhibition_rows):
    cmp = compare_inhibition_models(inhibition_rows(noise=0.2))
    assert any("cannot rule out" in note for note in cmp.notes)
    assert all(
        any("not proof" in note or "not a proven" in note or "evidence" in note for note in r.notes)
        for r in cmp.results
    )


def test_unknown_model_name_lists_the_alternatives(inhibition_rows):
    with pytest.raises(ValueError, match="unknown inhibition model"):
        fit_inhibition_series(inhibition_rows(), "suicide")


@pytest.mark.parametrize(
    "rows,message",
    [
        ([{"S": 1, "I": 0, "v": 1}, {"S": 2, "I": 0, "v": 2}], "at least six"),
        ([{"S": s, "I": 0, "v": 1.0} for s in (1, 2, 3, 4, 5, 6)], "two inhibitor concentrations"),
        (
            [{"S": s, "I": i, "v": 1.0} for s in (1, 2) for i in (0, 1, 2)],
            "three substrate concentrations",
        ),
        ([{"S": s, "I": i, "v": -1.0} for s in (1, 2, 3) for i in (0, 1)], "negative rates"),
    ],
)
def test_unusable_series_are_refused_with_a_reason(rows, message):
    with pytest.raises(ValueError, match=message):
        compare_inhibition_models(rows)


def test_malformed_rows_name_the_offending_row():
    rows = [{"S": 1, "I": 0, "v": 1.0}] * 5 + [{"S": 1, "I": 0}]
    with pytest.raises(ValueError, match="row 5"):
        compare_inhibition_models(rows)


def test_comparison_serialises_to_a_flat_dict(inhibition_rows):
    payload = compare_inhibition_models(inhibition_rows(noise=0.2)).as_dict()
    assert set(payload) == {"best", "delta_aicc", "models", "notes"}
    assert isinstance(payload["models"], list) and len(payload["models"]) >= 2
    assert all(isinstance(v, float) for v in payload["delta_aicc"].values())


def test_a_total_convergence_failure_is_an_error_not_an_empty_report(monkeypatch):
    """If every mechanism fails, that is "no answer" (exit 1), not a blank table."""
    import kinetics.inhibition as inhibition_module

    def always_fails(rows, model, *, p0=None):
        raise ValueError("simulated convergence failure")

    monkeypatch.setattr(inhibition_module, "fit_inhibition_series", always_fails)
    rows = [{"S": S, "I": I, "v": 1.0} for S in (1.0, 2.0, 5.0) for I in (0.0, 1.0, 2.0)]
    with pytest.raises(ValueError, match="no mechanism could be fitted"):
        inhibition_module.compare_inhibition_models(rows)
