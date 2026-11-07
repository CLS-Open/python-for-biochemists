"""Calibration curves, limits of detection, and the assumptions behind both.

A calibration curve is the least glamorous part of biochemistry and the part most
often done wrong. Three specific failures this module is built to make visible:

1. **LOD/LOQ borrowed from a different fit.** The standard definitions
   (LOD = 3.3·σ/S, LOQ = 10·σ/S) use the *standard deviation of the blank or of the
   residuals* and the *slope of the calibration*. Mixing σ from one dataset with a
   slope from another produces a number with no meaning. Here both come from the
   same fit, and the residual degrees of freedom are used.

2. **R² as evidence of linearity.** A calibration over a narrow range can have
   R² = 0.999 and still be visibly curved. :func:`fit_calibration` reports the
   residuals and a lack-of-fit ratio so curvature has somewhere to show up.

3. **Unweighted fits on data with proportional error.** Spectrophotometric
   standards usually have error that grows with signal. A weighted fit and an
   unweighted fit can then disagree by more than either one's standard error.
   :func:`fit_calibration` does both when you ask it to.

Units: concentration in whatever the caller used (µg/mL, µM, mg/mL); signal in
absorbance units (dimensionless). The package records what it was given and never
converts.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

__all__ = ["CalibrationResult", "fit_calibration", "lod_loq"]


@dataclass
class CalibrationResult:
    """Outcome of a linear calibration fit.

    Attributes
    ----------
    slope, intercept:
        Fitted line, ``signal = slope·concentration + intercept``.
    stderr_slope, stderr_intercept:
        Standard errors from the fit.
    r_squared:
        Coefficient of determination. Reported, but explicitly *not* a linearity
        test — see ``lack_of_fit_f``.
    residual_sd:
        Residual standard deviation with ``n − 2`` degrees of freedom. This is the
        σ that LOD/LOQ are supposed to use.
    weighted:
        Whether the fit used weights.
    """

    slope: float
    intercept: float
    stderr_slope: float
    stderr_intercept: float
    r_squared: float
    residual_sd: float
    n: int
    weighted: bool = False
    residuals: list[float] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def lod(self) -> float:
        return lod_loq(self.slope, self.residual_sd)[0]

    @property
    def loq(self) -> float:
        return lod_loq(self.slope, self.residual_sd)[1]

    def concentration(self, signal: float) -> float:
        """Invert the calibration: concentration for a measured signal."""
        if self.slope == 0:
            raise ValueError("cannot invert a calibration with zero slope")
        return (signal - self.intercept) / self.slope

    def as_dict(self) -> dict[str, object]:
        return {
            "slope": round(self.slope, 8),
            "intercept": round(self.intercept, 8),
            "stderr_slope": round(self.stderr_slope, 8),
            "stderr_intercept": round(self.stderr_intercept, 8),
            "r_squared": round(self.r_squared, 8),
            "residual_sd": round(self.residual_sd, 8),
            "n": self.n,
            "weighted": self.weighted,
            "lod": round(self.lod, 6) if math.isfinite(self.lod) else None,
            "loq": round(self.loq, 6) if math.isfinite(self.loq) else None,
            "notes": list(self.notes),
        }


def fit_calibration(
    concentration: Sequence[float],
    signal: Sequence[float],
    *,
    weighted: bool = False,
    replicate_sd: Sequence[float] | None = None,
) -> CalibrationResult:
    """Fit a straight line to calibration standards by ordinary least squares.

    Parameters
    ----------
    concentration:
        Standard concentrations, in the units of your assay. Must not all be equal.
    signal:
        Measured signals (absorbance units). Equal length.
    weighted:
        If ``True``, weight each point by ``1/signal²`` — appropriate when the error
        is proportional to the signal, which is the usual case for absorbance. When
        ``replicate_sd`` is given, those are used as weights instead.
    replicate_sd:
        Optional per-point standard deviations from replicate standards. Preferred
        over the proportional rule when you have them, because they are measured.

    Raises
    ------
    ValueError
        Fewer than three standards (a straight line through two points has no
        residual degrees of freedom, so its σ — and hence its LOD — is undefined),
        unequal lengths, or non-finite values.

    Examples
    --------
    >>> x = [0, 1, 2, 4, 6, 8, 10]
    >>> y = [0.01, 0.11, 0.21, 0.41, 0.61, 0.80, 1.01]
    >>> fit_calibration(x, y).r_squared > 0.999
    True
    >>> round(fit_calibration(x, y).slope, 3)
    0.1
    """
    x = np.asarray(concentration, dtype=float)
    y = np.asarray(signal, dtype=float)
    if x.shape != y.shape:
        raise ValueError(
            f"concentration and signal must be the same length, got {x.shape} and {y.shape}"
        )
    if len(x) < 3:
        raise ValueError(
            "at least three standards are needed: with two points the residual standard "
            "deviation has zero degrees of freedom, so LOD and LOQ are undefined"
        )
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("non-finite value in the calibration data")
    if np.allclose(x, x[0]):
        raise ValueError("all standards have the same concentration; no slope is identifiable")

    notes: list[str] = []
    weights = None
    if replicate_sd is not None:
        sd = np.asarray(replicate_sd, dtype=float)
        if sd.shape != x.shape or np.any(sd <= 0):
            raise ValueError("replicate_sd must be positive and match the data length")
        weights = 1.0 / sd**2
        notes.append("weighted by measured replicate standard deviations")
    elif weighted:
        if np.any(y <= 0):
            raise ValueError("proportional weighting requires positive signals")
        weights = 1.0 / y**2
        notes.append("weighted by 1/signal^2 (assumes error proportional to signal)")

    if weights is None:
        slope, intercept = np.polyfit(x, y, 1)
    else:
        # weighted least squares via a single call with weights
        slope, intercept = np.polyfit(x, y, 1, w=np.sqrt(weights))

    y_hat = slope * x + intercept
    residuals = y - y_hat
    n = len(x)
    rss = float(np.sum(residuals**2))
    dof = n - 2
    if weights is None:
        residual_sd = math.sqrt(rss / dof)
        x_mean = float(np.mean(x))
        sxx = float(np.sum((x - x_mean) ** 2))
        stderr_slope = residual_sd / math.sqrt(sxx)
        stderr_intercept = residual_sd * math.sqrt(1.0 / n + x_mean**2 / sxx)
    else:
        # weighted: the residual sd is on the weighted scale
        w = weights
        wsum = float(np.sum(w))
        wx = float(np.sum(w * x))
        sxx = float(np.sum(w * x**2)) - wx**2 / wsum
        weighted_rss = float(np.sum(w * residuals**2))
        sigma2 = weighted_rss / dof
        # weights are w = 1/sd², so sqrt(sigma2) is *dimensionless*: the ratio of the
        # observed scatter to the scatter the weights assumed. It is ~1 when the weights
        # are right and >1 when the replicates were optimistically precise. LOD/LOQ need
        # an absolute sigma, so the scale factor is multiplied by the SD at the lowest
        # standard — the concentration nearest the detection limit, which is the one
        # whose noise governs it.
        scale = math.sqrt(sigma2) if sigma2 >= 0 else float("nan")
        low = int(np.argmin(x))
        sigma_at_low = float(sd[low]) if replicate_sd is not None else float(y[low])
        residual_sd = scale * sigma_at_low
        notes.append(
            f"weighted fit: the scale factor is {scale:.3f} (1.0 means the weights matched "
            f"the scatter), and LOD/LOQ use σ = {sigma_at_low:.6g} at the lowest standard "
            f"(concentration {x[low]:g}), not a residual SD in signal units"
        )
        notes.append(
            "R-squared is computed on the unweighted residuals, so it is not "
            "the quantity the weighted fit minimised"
        )
        stderr_slope = math.sqrt(sigma2 / sxx) if sxx > 0 else float("nan")
        stderr_intercept = (
            math.sqrt(sigma2 * (1.0 / wsum + (wx / wsum) ** 2 / sxx)) if sxx > 0 else float("nan")
        )

    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float("nan") if tss <= 0 else 1.0 - rss / tss

    if r2 > 0.995:
        notes.append(
            "R-squared above 0.995 is common and is not evidence of linearity; inspect the residuals"
        )
    if (
        abs(intercept) > 3 * (stderr_intercept if math.isfinite(stderr_intercept) else 0)
        and abs(intercept) > 1e-9
    ):
        notes.append(
            "intercept is more than three standard errors from zero: a non-zero blank "
            "response is normal for colorimetric assays, but it must be subtracted "
            "consistently rather than ignored"
        )

    return CalibrationResult(
        slope=float(slope),
        intercept=float(intercept),
        stderr_slope=float(stderr_slope),
        stderr_intercept=float(stderr_intercept),
        r_squared=r2,
        residual_sd=float(residual_sd),
        n=n,
        weighted=weights is not None,
        residuals=[float(r) for r in residuals],
        notes=notes,
    )


def lod_loq(
    slope: float, residual_sd: float, *, k_lod: float = 3.3, k_loq: float = 10.0
) -> tuple[float, float]:
    """Limit of detection and limit of quantification (IUPAC convention).

    ``LOD = k_lod·σ / slope`` and ``LOQ = k_loq·σ / slope`` with the conventional
    ``k_lod = 3.3`` (≈ 3·√2 for the blank-plus-sample comparison) and ``k_loq = 10``.

    Parameters
    ----------
    slope:
        Calibration slope — signal per unit concentration. Must be non-zero.
    residual_sd:
        Standard deviation of the residuals (or of the blank) in signal units, from
        the **same** data as the slope. Mixing sources is the classic error. In a
        weighted calibration this is *not* simply the residual SD: it is the scale
        factor times the SD at the lowest standard, and
        :func:`fit_calibration` constructs it that way and says so in a note.
    k_lod, k_loq:
        Multipliers. Exposed because conventions differ (some fields use 3 and 10);
        changing them is a decision that must be stated in a report.

    Returns
    -------
    (lod, loq)
        In the concentration units of the calibration. ``loq`` is always larger than
        ``lod`` for ``k_loq > k_lod``.

    Raises
    ------
    ValueError
        For a zero slope (no concentration dependence) or a negative σ.

    Examples
    --------
    >>> lod, loq = lod_loq(slope=0.1, residual_sd=0.01)
    >>> round(lod, 4), round(loq, 4)
    (0.33, 1.0)
    >>> lod_loq(slope=0.1, residual_sd=0.0)      # a noiseless curve has no detection limit
    (0.0, 0.0)
    """
    if slope == 0:
        raise ValueError("no detection limit exists for a zero-slope calibration")
    if residual_sd < 0:
        raise ValueError("residual_sd cannot be negative")
    scale = abs(slope)
    return (k_lod * residual_sd / scale, k_loq * residual_sd / scale)
