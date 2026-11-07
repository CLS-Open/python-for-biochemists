"""Inhibition typing by global fitting and AICc comparison.

How inhibition is *often* taught: draw a Lineweaver–Burk plot at several inhibitor
concentrations and look at where the lines meet. How it is done here: fit the whole
(S, I, v) surface with each candidate mechanism and compare the fits with AICc.

The difference matters because the graphical method asks a human to judge
intersections on a reciprocal plot — the least reliable transformation in the
toolkit (see :mod:`kinetics.fitting`) — while the global fit uses every point
simultaneously and reports how much better one mechanism explains the data than
another. Neither approach *proves* a mechanism, and the docstrings say so.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np

from .fitting import FitResult, fit_nonlinear
from .models import INHIBITION_MODELS

__all__ = ["InhibitionComparison", "compare_inhibition_models", "fit_inhibition_series"]


@dataclass
class InhibitionComparison:
    """Result of comparing inhibition mechanisms on one dataset.

    Attributes
    ----------
    results:
        One :class:`~kinetics.fitting.FitResult` per mechanism, sorted by AICc
        ascending (best first).
    best:
        The lowest-AICc mechanism name. A *best* model is not a *correct* model.
    delta_aicc:
        ``AICc − min(AICc)`` per model. Conventional reading: Δ ≤ 2 is
        indistinguishable, 4–7 is considerably less support, > 10 is effectively no
        support.
    notes:
        Caveats that must travel with the numbers.
    """

    results: list[FitResult]
    best: str
    delta_aicc: dict[str, float]
    notes: list[str]

    def as_dict(self) -> dict[str, object]:
        return {
            "best": self.best,
            "delta_aicc": {k: round(v, 4) for k, v in self.delta_aicc.items()},
            "models": [r.as_dict() for r in self.results],
            "notes": list(self.notes),
        }


def _prepare(rows: Iterable[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    S: list[float] = []
    I: list[float] = []
    v: list[float] = []
    for i, row in enumerate(rows):
        try:
            S.append(float(row["S"]))
            I.append(float(row["I"]))
            v.append(float(row["v"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"row {i} is not a usable (S, I, v) record: {row!r} ({exc})") from None
    S_arr, I_arr, v_arr = np.asarray(S), np.asarray(I), np.asarray(v)
    if len(v_arr) < 6:
        raise ValueError(
            "an inhibition series needs at least six (S, I, v) points to fit even the "
            "three-parameter mechanisms; six is a floor, not a recommendation"
        )
    if len(np.unique(I_arr)) < 2:
        raise ValueError("at least two inhibitor concentrations are required to estimate Ki")
    if len(np.unique(S_arr)) < 3:
        raise ValueError("at least three substrate concentrations are required")
    if np.any(v_arr < 0):
        raise ValueError("negative rates cannot be fitted by these mechanisms")
    return S_arr, I_arr, v_arr


def _bind_inhibitor(func, I):
    """Close over ``I`` so the fitted parameters keep their declared names.

    The models have the signature ``f(S, vmax, km, I, ki)`` — the inhibitor sits
    *between* the fitted parameters, which is right for a reader and wrong for a
    fitter that appends extra arguments at the end. Passing ``extra_args=(I,)``
    silently bound the inhibitor vector to ``ki`` and ``ki`` to ``I``; the fit then
    failed with a numpy broadcasting message rather than a useful one. Binding by
    name here removes the whole class of mistake, and ``test_inhibition.py`` checks
    that a fit recovers parameters passed to the forward model.
    """

    def bound(S, vmax, km, ki):
        return func(S, vmax, km, I, ki)

    bound.__name__ = getattr(func, "__name__", "inhibition_model") + "_bound"
    return bound


def fit_inhibition_series(
    rows: Iterable[dict],
    model: str,
    *,
    p0: Sequence[float] | None = None,
) -> FitResult:
    """Fit one inhibition mechanism to a full (S, I, v) series.

    Parameters
    ----------
    rows:
        Iterable of dicts with keys ``S``, ``I``, ``v``. Tidy data, one row per
        measurement — the format ``datasets/sample/inhibition_series.csv`` uses.
    model:
        One of ``"competitive"``, ``"noncompetitive"``, ``"uncompetitive"``.
    p0:
        Optional starting guesses ``(vmax, km, ki)``. The default scales the guesses
        from the data (Vmax ≈ max rate, Km ≈ median S, Ki ≈ median positive I), which
        is documented here rather than hidden in a constant.

    Raises
    ------
    ValueError
        Unknown model name, malformed rows, too few points or unique conditions.

    Examples
    --------
    >>> rows = [{"S": s, "I": i, "v": 10 * s / (25 * (1 + i / 8) + s)}
    ...         for s in (5, 25, 100) for i in (0, 8)]
    >>> res = fit_inhibition_series(rows, "competitive")
    >>> res.k
    3
    """
    if model not in INHIBITION_MODELS:
        raise ValueError(
            f"unknown inhibition model {model!r}; choose from {sorted(INHIBITION_MODELS)}"
        )
    S, I, v = _prepare(rows)
    spec = INHIBITION_MODELS[model]
    if p0 is None:
        p0 = (
            float(np.max(v)) * 1.2,
            float(np.median(S)),
            float(np.median(I[I > 0])) if np.any(I > 0) else 1.0,
        )
    return fit_nonlinear(
        S,
        v,
        _bind_inhibitor(spec["func"], I),
        p0,
        spec["params"],
        name=model,
        notes=[
            f"mechanism: {spec['description']}",
            "a lower AICc than a rival mechanism is evidence, not proof, of that mechanism",
        ],
    )


def compare_inhibition_models(rows: Iterable[dict]) -> InhibitionComparison:
    """Fit all three mechanisms and rank them by AICc.

    The same rows are used for every model, so the comparison is like-for-like. If a
    model fails to converge it is omitted and a note says so — silently dropping it
    would make the remaining models look better than they are.

    Examples
    --------
    >>> rows = [{"S": s, "I": i, "v": 10 * s / (25 * (1 + i / 8) + s)}
    ...         for s in (5, 10, 25, 50, 100) for i in (0, 4, 16)]
    >>> cmp = compare_inhibition_models(rows)
    >>> cmp.best
    'competitive'
    """
    rows = list(rows)
    # Validate once, before the loop. These are properties of the data, identical for
    # every mechanism, so they must reach the caller as a data error rather than being
    # folded into "no mechanism could be fitted" — which is what happened when the
    # per-model ``except ValueError`` swallowed them.
    _prepare(rows)
    results: list[FitResult] = []
    notes: list[str] = []
    for name in sorted(INHIBITION_MODELS):
        try:
            results.append(fit_inhibition_series(rows, name))
        except ValueError as exc:
            notes.append(f"{name} could not be fitted and was excluded: {exc}")
    if not results:
        raise ValueError("no mechanism could be fitted to these data")
    results.sort(key=lambda r: r.aicc)
    best = results[0].aicc
    if math.isinf(best):
        # A mechanism that reproduces every point exactly (RSS = 0) makes AICc -inf and
        # every difference undefined. Reporting "ΔAICc = inf" for all models would read
        # like a decisive result, so say what actually happened.
        delta = {r.model: (0.0 if r is results[0] else float("inf")) for r in results}
        notes.append(
            "at least one mechanism fits these points exactly (RSS = 0), so AICc is "
            "-inf and the differences between mechanisms are undefined: these data "
            "carry no error term. Comparison needs measured or simulated replicates."
        )
        close = []
    else:
        delta = {r.model: (r.aicc - best if np.isfinite(r.aicc) else float("inf")) for r in results}
        close = [m for m, d in delta.items() if 0 < d <= 2.0]
    if close:
        notes.append(
            "indistinguishable models (ΔAICc ≤ 2): "
            + ", ".join(sorted(close))
            + " — report the ambiguity rather than the winner"
        )
    notes.append(
        "AICc compares these data under these mechanisms; it cannot rule out a mechanism "
        "that was not among the candidates"
    )
    return InhibitionComparison(
        results=results, best=results[0].model, delta_aicc=delta, notes=notes
    )
