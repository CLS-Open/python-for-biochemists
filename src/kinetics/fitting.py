"""Fitting, model comparison and uncertainty.

Three design decisions, each of which a learner should be able to argue with:

1. **Nonlinear least squares is the default, and linearisation is available only as
   a diagnostic.** Lineweaver–Burk, Eadie–Hofstee and Hanes–Woolf are fitted here
   on purpose: they are what a generation of students were taught, they distort
   error structure, and the distortion is measurable. ``compare_lineararisations``
   makes it measurable in one call.

2. **Uncertainty is reported two ways.** Asymptotic standard errors come free from
   the fit and are optimistic when the model is nonlinear or *n* is small. The
   bootstrap resamples the data and gives an interval you can defend. When the two
   disagree badly, the asymptotic numbers are the ones to distrust.

3. **Model comparison uses AICc, not R².** R² always rises when you add a
   parameter; AICc penalises it, and the correction matters at the sample sizes
   real enzyme assays have (n = 10–30).

Everything here is deterministic except :func:`bootstrap_ci`, whose seed is a
required argument and is recorded in the result.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit

from .models import michaelis_menten

__all__ = [
    "FitResult",
    "aic",
    "aicc",
    "bootstrap_ci",
    "compare_linearisations",
    "fit_eadie_hofstee",
    "fit_hanes",
    "fit_lineweaver_burk",
    "fit_michaelis_menten",
    "fit_nonlinear",
]


@dataclass
class FitResult:
    """The outcome of one fit, including how it was obtained and how sure it is.

    Parameters
    ----------
    model:
        Human-readable model name, e.g. ``"michaelis_menten"`` or
        ``"lineweaver_burk"``.
    params:
        Fitted parameters, in the model's documented order.
    stderr:
        Asymptotic standard errors (``float('nan')`` where SciPy could not
        estimate them, which happens when the model is over-parameterised).
    method:
        ``"nonlinear"``, ``"linearised"``, or ``"bootstrap"`` — never inferred
        later, always recorded.
    n:
        Number of data points used.
    k:
        Number of fitted parameters.
    rss:
        Residual sum of squares (in the units of the response, squared).
    notes:
        Anything a reader needs in order not to over-read the numbers.
    """

    model: str
    params: dict[str, float]
    stderr: dict[str, float] = field(default_factory=dict)
    method: str = "nonlinear"
    n: int = 0
    k: int = 0
    rss: float = float("nan")
    notes: list[str] = field(default_factory=list)
    extra: dict[str, object] = field(default_factory=dict)

    # -- derived quantities -------------------------------------------------
    @property
    def aic(self) -> float:
        return aic(self.rss, self.n, self.k)

    @property
    def aicc(self) -> float:
        return aicc(self.rss, self.n, self.k)

    @property
    def r_squared(self) -> float:
        """Coefficient of determination, computed from the residuals.

        Returned as ``nan`` when the response has no variance (a constant y gives
        a 0/0), instead of inventing a number.
        """
        tss = float(self.extra.get("tss", float("nan")))
        if not math.isfinite(tss) or tss <= 0:
            return float("nan")
        return 1.0 - self.rss / tss

    def aicc_or_none(self) -> float | None:
        """AICc, or ``None`` where the small-sample correction is undefined.

        ``aicc`` raises when ``n - k - 1 <= 0``, which is the right behaviour for a
        model comparison someone asked for — but a report should not explode because
        one fit in the table was tiny. Reports call this instead.
        """
        try:
            value = float(self.aicc)
        except ValueError:
            return None
        return value if math.isfinite(value) else None

    def _rounded_aicc(self) -> float | None:
        value = self.aicc_or_none()
        return None if value is None else round(value, 4)

    def as_dict(self) -> dict[str, object]:
        return {
            "model": self.model,
            "method": self.method,
            "params": {k: float(v) for k, v in self.params.items()},
            "stderr": {k: float(v) for k, v in self.stderr.items()},
            "n": self.n,
            "k": self.k,
            "rss": float(self.rss),
            "r_squared": None if not math.isfinite(self.r_squared) else round(self.r_squared, 6),
            "aic": round(self.aic, 4) if math.isfinite(self.aic) else None,
            "aicc": self._rounded_aicc(),
            "notes": list(self.notes),
        }


def aic(rss: float, n: int, k: int) -> float:
    """Akaike information criterion for Gaussian residuals: ``n·ln(RSS/n) + 2k``.

    ``RSS`` must be positive; a perfect fit has ``RSS = 0`` and an AIC of ``-inf``,
    which is correct but useless, so it is returned as ``-inf`` rather than a
    number. In practice a perfect fit means the data were generated without noise
    and no model comparison is meaningful.

    >>> round(aic(4.0, 10, 2), 6)
    -5.162907
    >>> round(aic(4.0, 10, 3), 6)          # one more parameter costs 2 units
    -3.162907
    """
    if n <= k:
        raise ValueError(f"need more data points than parameters (n={n}, k={k})")
    if rss < 0:
        raise ValueError("rss cannot be negative")
    if rss == 0:
        return float("-inf")
    return n * math.log(rss / n) + 2 * k


def aicc(rss: float, n: int, k: int) -> float:
    """Small-sample corrected AIC: ``AIC + 2k(k+1)/(n−k−1)``.

    Required when ``n/k < 40``, which is the normal situation in enzymology. The
    correction is undefined for ``n − k − 1 <= 0``; that is an error, not a
    warning, because a comparison computed from an undefined number is worse than
    no comparison.

    >>> round(aicc(4.0, 10, 2), 6)
    -3.448622
    >>> round(aicc(4.0, 10, 2) - aic(4.0, 10, 2), 6)   # the small-sample penalty
    1.714286
    >>> round(aicc(4.0, 4, 2), 6)                      # n - k - 1 == 1, penalty 12
    16.0
    >>> aicc(4.0, 3, 2)
    Traceback (most recent call last):
    ...
    ValueError: AICc undefined: n=3, k=2 leaves no degrees of freedom
    """
    if n - k - 1 <= 0:
        raise ValueError(f"AICc undefined: n={n}, k={k} leaves no degrees of freedom")
    base = aic(rss, n, k)
    if not math.isfinite(base):
        return base
    return base + (2 * k * (k + 1)) / (n - k - 1)


def _rss(y: Sequence[float], y_hat: Sequence[float]) -> float:
    y = np.asarray(y, dtype=float)
    y_hat = np.asarray(y_hat, dtype=float)
    return float(np.sum((y - y_hat) ** 2))


def fit_nonlinear(
    x: Sequence[float],
    y: Sequence[float],
    model: Callable[..., float],
    p0: Iterable[float],
    param_names: Sequence[str],
    *,
    name: str | None = None,
    sigma: Sequence[float] | None = None,
    extra_args: tuple = (),
    notes: Iterable[str] = (),
) -> FitResult:
    """Fit ``y = model(x, *params, *extra_args)`` by nonlinear least squares.

    Parameters
    ----------
    x, y:
        Data arrays of equal length. No filtering is done — an invalid point raises.
    p0:
        Starting guesses. Required, not guessed from the data: a fit that only
        converges from a lucky start is a fragile procedure and should be reported
        as such.
    sigma:
        Optional standard deviations for weighted least squares. Enzyme rates often
        have *proportional* error, in which case ``sigma`` should be proportional to
        ``y`` (see module 5).
    notes:
        Appended to the result so the caller cannot lose the caveats.

    Raises
    ------
    ValueError
        If the arrays differ in length, contain non-finite values, or have fewer
        points than parameters, or if the optimiser fails to converge.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError(f"x and y must have the same length, got {x.shape} and {y.shape}")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("x and y must be finite; filter or explain the missing points first")
    k = len(param_names)
    if len(x) <= k:
        raise ValueError(f"need more points ({len(x)}) than parameters ({k})")
    p0 = tuple(float(value) for value in p0)
    if len(p0) != k:
        raise ValueError(
            f"{k} parameter names {tuple(param_names)} but {len(p0)} starting values "
            f"{p0}; each parameter needs exactly one start"
        )
    if sigma is not None and np.asarray(sigma, dtype=float).shape != x.shape:
        raise ValueError("sigma must have one value per data point")
    if np.ptp(y) == 0.0:
        # A constant response carries no information about any parameter, but the
        # optimiser will happily "converge" on one arbitrary pair of values. Refusing is
        # the only honest answer; a bootstrap over such data otherwise returns an
        # interval of zero width, which looks like certainty.
        raise ValueError(
            f"the response is constant ({float(y[0])!r} at every point), so "
            f"{name or 'the model'} "
            "is not identifiable: there is nothing for the fit to explain"
        )

    def wrapped(xx, *params):
        return model(xx, *params, *extra_args)

    absolute = None if sigma is None else np.asarray(sigma, dtype=float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            popt, pcov = curve_fit(
                wrapped,
                x,
                y,
                p0=list(p0),
                sigma=absolute,
                absolute_sigma=absolute is not None,
                maxfev=20000,
            )
        except (RuntimeError, ValueError) as exc:
            raise ValueError(f"fit did not converge for {name or 'model'}: {exc}") from None
    covariance_warning = any("Covariance" in str(w.message) for w in caught)
    if covariance_warning:
        pcov = None

    perr = np.sqrt(np.abs(np.diag(pcov))) if pcov is not None else np.full(k, np.nan)
    residuals = y - wrapped(x, *popt)
    result = FitResult(
        model=name or getattr(model, "__name__", "model"),
        params={nm: float(v) for nm, v in zip(param_names, popt, strict=True)},
        stderr={nm: float(e) for nm, e in zip(param_names, perr, strict=True)},
        method="weighted_nonlinear" if sigma is not None else "nonlinear",
        n=len(x),
        k=k,
        rss=_rss(y, wrapped(x, *popt)),
        notes=list(notes),
        extra={"tss": float(np.sum((y - y.mean()) ** 2)), "residuals": residuals.tolist()},
    )
    if covariance_warning:
        result.notes.append(
            "the optimiser could not estimate the covariance matrix, so every standard "
            "error in this result is NaN: the parameters may be identified by the data but "
            "their uncertainty is not. Usually the model has more parameters than the data "
            "can separate"
        )
    if np.any(np.isnan(perr)) or np.any(perr > 10 * np.abs(popt)):
        result.notes.append(
            "at least one standard error is huge relative to its parameter: the model "
            "may be over-parameterised for this data, or the data may not identify it"
        )
    return result


def fit_michaelis_menten(
    S: Sequence[float],
    v: Sequence[float],
    *,
    sigma: Sequence[float] | None = None,
    p0: tuple[float, float] | None = None,
) -> FitResult:
    """Fit the Michaelis–Menten parameters by nonlinear least squares.

    The starting guesses (when not supplied) are the maximum observed rate and the
    substrate concentration nearest half of it — crude, monotone, and *documented*,
    which is better than a magic number.

    Examples
    --------
    With exact data the fit returns the parameters it was given, which is the only
    way to check a fitter without trusting it:

    >>> import numpy as np
    >>> S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    >>> exact = 10.0 * S / (25.0 + S)                     # Vmax 10, Km 25
    >>> res = fit_michaelis_menten(S, exact)
    >>> round(res.params["vmax"], 6), round(res.params["km"], 6)
    (10.0, 25.0)

    Noise moves both estimates, and the standard error says by how much *if the
    error model is right*:

    >>> rng = np.random.default_rng(7)
    >>> noisy = exact + 0.20 * rng.standard_normal(S.size)
    >>> res = fit_michaelis_menten(S, noisy)
    >>> round(res.params["vmax"], 2), round(res.params["km"], 2)
    (10.1, 26.93)
    """
    S_arr = np.asarray(S, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    # Non-positive rates are flagged, not refused. At the bottom of a substrate series a
    # negative rate is usually the detector floor, and dropping those points biases Vmax
    # upwards; refusing them makes the tool unusable on real assays (this package's own
    # cooperative_hill dataset contains one, at [S] = 0.5 uM). A negative rate far above the
    # floor is a different thing -- a data error -- and only a human can tell the two apart,
    # so the fit carries the count and the worst value instead of deciding for the reader.
    notes: list[str] = ["Vmax is an extrapolation to infinite [S] and is never measured directly"]
    if np.any(v_arr <= 0):
        n_bad = int(np.sum(v_arr <= 0))
        notes.append(
            f"{n_bad} of {len(v_arr)} rates are not positive (lowest {float(np.min(v_arr)):g}). "
            "Near the bottom of a substrate series this is the detector floor and the points "
            "are kept; a non-positive rate well above the floor is a data error, and no fit "
            "can tell you which one you have"
        )
    if p0 is None:
        vmax0 = float(np.max(v_arr)) if len(v_arr) else 1.0
        half = vmax0 / 2.0
        km0 = float(S_arr[np.argmin(np.abs(v_arr - half))]) if len(S_arr) else 1.0
        p0 = (max(vmax0, 1e-9), max(km0, 1e-9))
    return fit_nonlinear(
        S_arr,
        v_arr,
        michaelis_menten,
        p0,
        ("vmax", "km"),
        name="michaelis_menten",
        sigma=sigma,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Linearisations: included to be measured, not to be recommended
# ---------------------------------------------------------------------------


def _linear_fit(x, y, model_name: str, notes: list[str]) -> FitResult:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        raise ValueError("a linear fit of fewer than three points is not identifiable")
    sxx = float(np.sum((x - x.mean()) ** 2))
    if sxx == 0.0:
        raise ValueError(
            "every x value is identical, so the slope is not identifiable; this is a "
            "blocked experiment, not a badly conditioned one. Check the transformation "
            "(Eadie-Hofstee divides by S, Hanes-Woolf divides by v)"
        )
    slope, intercept = np.polyfit(x, y, 1)
    y_hat = slope * x + intercept
    residuals = y - y_hat
    var = float(np.sum(residuals**2) / (len(x) - 2))
    # "Zero" residual variance never happens in floating point: a perfect straight line
    # through four points leaves about 5e-31. Compare it with the scale of the response
    # instead, and say what a meaningless standard error means rather than printing
    # 7e-16 as though it were precision.
    noise_floor = (1e-10 * float(np.max(np.abs(y)))) ** 2 if len(y) else 0.0
    if var <= noise_floor:
        se_slope = float("nan")
        notes = [
            *notes,
            (
                "the residuals are at the level of floating-point noise, so the "
                "slope's standard error cannot be estimated from these data: they "
                "were probably generated without an error term"
            ),
        ]
    else:
        se_slope = math.sqrt(var / sxx)
    return FitResult(
        model=model_name,
        params={"slope": float(slope), "intercept": float(intercept)},
        stderr={"slope": se_slope, "intercept": float("nan")},
        method="linearised",
        n=len(x),
        k=2,
        rss=float(np.sum(residuals**2)),
        notes=notes,
        extra={"tss": float(np.sum((y - y.mean()) ** 2)), "residuals": residuals.tolist()},
    )


_LB_NOTE = (
    "Lineweaver-Burk plots 1/v against 1/S. It linearises the hyperbola but inverts "
    "the error structure: low-[S] points, which have the least reliable rates, get the "
    "largest weight, and the reciprocal inflates their errors quadratically."
)
_EH_NOTE = (
    "Eadie-Hofstee plots v against v/S. It avoids double reciprocals but puts the "
    "measured rate on both axes, so errors in v are correlated with the x variable."
)
_HW_NOTE = (
    "Hanes-Woolf plots S/v against S, the least damaging of the three linearisations, "
    "because the x variable is the error-free one."
)


def fit_lineweaver_burk(S: Sequence[float], v: Sequence[float]) -> FitResult:
    """Linearised Lineweaver–Burk fit. Slope = Km/Vmax, intercept = 1/Vmax.

    Kept because it is still taught and still published; the result carries a note
    explaining its bias, and :func:`compare_linearisations` measures that bias.

    Examples
    --------
    >>> import numpy as np
    >>> S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    >>> rng = np.random.default_rng(7)
    >>> v = 10.0 * S / (25.0 + S) + 0.20 * rng.standard_normal(S.size)   # Vmax 10, Km 25
    >>> res = fit_lineweaver_burk(S, v)
    >>> round(res.params["vmax"], 2), round(res.params["km"], 2)
    (9.65, 23.73)

    Against a truth of 10 and 25. The nonlinear fit on the same data returns 10.1 and
    26.93; the difference is the bias the transformation introduces, measured rather
    than asserted.

    Data that are exactly proportional have no finite Vmax, and the fit says so
    instead of returning a nonsense number:

    >>> prop = fit_lineweaver_burk(S, 2.0 * S)
    >>> prop.params["vmax"], prop.params["km"]
    (inf, nan)
    """
    S_arr = np.asarray(S, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    if np.any(S_arr <= 0) or np.any(v_arr <= 0):
        raise ValueError(
            "Lineweaver-Burk requires positive [S] and positive rates (it reciprocates both)"
        )
    res = _linear_fit(1 / S_arr, 1 / v_arr, "lineweaver_burk", [_LB_NOTE])
    slope = res.params["slope"]
    intercept = res.params["intercept"]
    y_scale = float(np.max(np.abs(1.0 / v_arr))) if len(v_arr) else 1.0
    if abs(intercept) <= 1e-12 * y_scale:
        res.params["vmax"] = float("inf")
        res.params["km"] = float("nan")
        res.stderr["vmax"] = float("nan")
        res.stderr["km"] = float("nan")
        res.notes = [
            *res.notes,
            (
                "the fitted intercept is zero to machine precision, so "
                "Vmax is unbounded and no Km exists: the data are "
                "proportional (v/S constant), which double-reciprocal "
                "plotting cannot represent. The keys are kept so "
                "downstream code does not meet a KeyError."
            ),
        ]
    else:
        res.params["vmax"] = 1.0 / intercept
        res.params["km"] = slope / intercept
        res.stderr["vmax"] = float("nan")
        res.stderr["km"] = float("nan")
    return res


def fit_eadie_hofstee(S: Sequence[float], v: Sequence[float]) -> FitResult:
    """Eadie–Hofstee fit: ``v = Vmax − Km·(v/S)``.

    The regression of ``v`` on ``v/S`` gives ``Vmax`` as the intercept and ``−Km``
    as the slope.

    Examples
    --------
    >>> import numpy as np
    >>> S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    >>> rng = np.random.default_rng(7)
    >>> v = 10.0 * S / (25.0 + S) + 0.20 * rng.standard_normal(S.size)   # Vmax 10, Km 25
    >>> res = fit_eadie_hofstee(S, v)
    >>> round(res.params["vmax"], 2), round(res.params["km"], 2)
    (9.54, 23.54)
    """
    S_arr = np.asarray(S, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    if np.any(S_arr <= 0):
        raise ValueError(
            "Eadie-Hofstee divides by [S]; all substrate concentrations must be positive"
        )
    res = _linear_fit(v_arr / S_arr, v_arr, "eadie_hofstee", [_EH_NOTE])
    res.params["vmax"] = res.params["intercept"]
    res.params["km"] = -res.params["slope"]
    return res


def fit_hanes(S: Sequence[float], v: Sequence[float]) -> FitResult:
    """Hanes–Woolf fit: ``S/v = S/Vmax + Km/Vmax``.

    Examples
    --------
    >>> import numpy as np
    >>> S = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    >>> rng = np.random.default_rng(7)
    >>> v = 10.0 * S / (25.0 + S) + 0.20 * rng.standard_normal(S.size)   # Vmax 10, Km 25
    >>> res = fit_hanes(S, v)
    >>> round(res.params["vmax"], 2), round(res.params["km"], 2)
    (9.97, 25.48)

    The closest of the three to the nonlinear answer (10.1, 26.93), which is what you
    would expect: ``S`` is the variable on the x-axis and ``S`` is the one without
    measurement error.
    """
    S_arr = np.asarray(S, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    if np.any(v_arr <= 0):
        raise ValueError("Hanes-Woolf divides by v; all rates must be positive")
    res = _linear_fit(S_arr, S_arr / v_arr, "hanes", [_HW_NOTE])
    slope = res.params["slope"]
    intercept = res.params["intercept"]
    if slope != 0:
        res.params["vmax"] = 1.0 / slope
        res.params["km"] = intercept / slope
    return res


def compare_linearisations(S: Sequence[float], v: Sequence[float]) -> list[FitResult]:
    """Fit the same data four ways and return them in a fixed order.

    The order is nonlinear, Hanes, Eadie–Hofstee, Lineweaver–Burk: best method
    first. Returning a list rather than a table keeps the function free of
    formatting, which is :mod:`kinetics.report`'s job.

    Examples
    --------
    >>> S = [1, 2, 5, 10, 20, 50, 100]
    >>> v = [0.95, 1.8, 3.4, 5.0, 6.6, 8.0, 8.7]
    >>> [r.model for r in compare_linearisations(S, v)]
    ['michaelis_menten', 'hanes', 'eadie_hofstee', 'lineweaver_burk']
    """
    return [
        fit_michaelis_menten(S, v),
        fit_hanes(S, v),
        fit_eadie_hofstee(S, v),
        fit_lineweaver_burk(S, v),
    ]


def bootstrap_ci(
    S: Sequence[float],
    v: Sequence[float],
    *,
    n_resamples: int = 2000,
    seed: int = 20261106,
    confidence: float = 0.95,
) -> dict[str, tuple[float, float]]:
    """Percentile bootstrap confidence intervals for ``Vmax`` and ``Km``.

    Resamples *pairs* (S, v) with replacement, refits, and reports percentiles. This
    answers "how much would my estimate move if I repeated the experiment?" in a way
    asymptotic standard errors do not, and it makes no assumption that the
    parameters are normally distributed — which, for a ratio like Km, they are not.

    Parameters
    ----------
    n_resamples:
        Number of resamples. 2000 gives stable 2.5/97.5 percentiles; it costs about
        a second on a laptop.
    seed:
        Required and recorded. A "random" interval that cannot be reproduced is not
        a result.

    Returns
    -------
    dict
        ``{"vmax": (low, high), "km": (low, high)}`` in the input units.

    Examples
    --------
    >>> S = [1, 5, 10, 25, 50, 100]
    >>> v = [1.0, 3.3, 5.0, 6.9, 8.1, 8.6]
    >>> ci = bootstrap_ci(S, v, n_resamples=200, seed=7)
    >>> ci["vmax"][0] < ci["vmax"][1]
    True
    """
    S_arr = np.asarray(S, dtype=float)
    v_arr = np.asarray(v, dtype=float)
    if len(S_arr) < 4:
        raise ValueError(
            "bootstrapping fewer than four points produces an interval that means nothing"
        )
    rng = np.random.default_rng(seed)
    vmaxes: list[float] = []
    kms: list[float] = []
    failures = 0
    for _ in range(n_resamples):
        idx = rng.integers(0, len(S_arr), len(S_arr))
        try:
            res = fit_michaelis_menten(S_arr[idx], v_arr[idx])
        except ValueError:
            failures += 1
            continue
        vmaxes.append(res.params["vmax"])
        kms.append(res.params["km"])
    if not vmaxes:
        raise ValueError(
            "every bootstrap resample failed to fit; the data do not identify the model"
        )
    alpha = (1.0 - confidence) / 2.0
    lo, hi = 100 * alpha, 100 * (1 - alpha)
    return {
        "vmax": (float(np.percentile(vmaxes, lo)), float(np.percentile(vmaxes, hi))),
        "km": (float(np.percentile(kms, lo)), float(np.percentile(kms, hi))),
    }
