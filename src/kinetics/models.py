"""Rate laws, written out in full.

Each function takes numbers in, returns numbers out, and does nothing clever. That
is deliberate: a rate law you cannot read is a rate law you cannot check, and every
model below is something a learner should be able to derive on paper.

Units
-----
================  ==============================================
quantity          unit used throughout this course
================  ==============================================
``S``             substrate concentration, µM
``v``             initial rate, µmol min⁻¹ (mg protein)⁻¹
``vmax``          same as ``v``
``km``, ``ki``    µM
``I``             inhibitor concentration, µM
``n`` (Hill)      dimensionless
================  ==============================================

The package does not convert units. If you feed it mM you will get a Km in mM and
a Vmax in whatever your rate was in — which is fine, provided you say so. Module 1
is about why that sentence is the whole module.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "INHIBITION_MODELS",
    "competitive",
    "hill",
    "michaelis_menten",
    "mixed_inhibition",
    "noncompetitive",
    "uncompetitive",
]


def michaelis_menten(S, vmax, km):
    """Michaelis–Menten rate law: ``v = Vmax·S / (Km + S)``.

    Parameters
    ----------
    S : float or array
        Substrate concentration (µM). Must be non-negative.
    vmax : float
        Maximum rate (same units as the returned rate).
    km : float
        Michaelis constant (µM), strictly positive.

    Returns
    -------
    float or ndarray

    Raises
    ------
    ValueError
        If ``km <= 0`` (the model has no finite Km at zero), or if any ``S < 0``.

    Notes
    -----
    The curve is *rectangular hyperbolic*: it approaches Vmax asymptotically and
    therefore **Vmax can never be measured directly** — only estimated by
    extrapolation, which is why the estimate's uncertainty is larger than a fit's
    goodness-of-fit suggests.

    Examples
    --------
    >>> round(michaelis_menten(25.0, 10.0, 25.0), 6)
    5.0
    >>> round(michaelis_menten(1000.0, 10.0, 1.0), 4)
    9.99
    >>> michaelis_menten(0.0, 10.0, 25.0)
    0.0
    """
    S = np.asarray(S, dtype=float)
    if km <= 0:
        raise ValueError(f"km must be positive, got {km!r}")
    if np.any(S < 0):
        raise ValueError("substrate concentration cannot be negative")
    return _out(vmax * S / (km + S))


def hill(S, vmax, k_half, n=1.0):
    """Hill equation: ``v = Vmax·Sⁿ / (K½ⁿ + Sⁿ)``.

    ``n = 1`` reduces to Michaelis–Menten. A fit with ``n > 1`` is often reported as
    "positive cooperativity", which is an *interpretation*, not a measurement: the
    same curve is produced by a monomer with a bad assay. Module 5 says so.

    Examples
    --------
    >>> round(hill(25.0, 10.0, 25.0, 1.0), 6)
    5.0
    >>> round(hill(50.0, 10.0, 25.0, 2.0), 3)   # steeper than MM at S = 2·K
    8.0
    """
    S = np.asarray(S, dtype=float)
    if k_half <= 0:
        raise ValueError(f"k_half must be positive, got {k_half!r}")
    if np.any(S < 0):
        raise ValueError("substrate concentration cannot be negative")
    Sn = np.power(S, n)
    return _out(vmax * Sn / (np.power(k_half, n) + Sn))


def _out(value):
    """Return a plain ``float`` for scalar input, an array otherwise.

    Without this, ``michaelis_menten(25.0, 10.0, 25.0)`` returns ``np.float64``,
    which reprs as ``np.float64(5.0)`` in numpy 2 and refuses ``json.dumps``. The
    models are also used as teaching text: the number printed in a docstring should
    be the number a student types at the prompt.
    """
    arr = np.asarray(value, dtype=float)
    return float(arr) if arr.ndim == 0 else arr


def _check_inhibitor(I, ki):
    I = np.asarray(I, dtype=float)
    if np.any(I < 0):
        raise ValueError("inhibitor concentration cannot be negative")
    if ki <= 0:
        raise ValueError(f"ki must be positive, got {ki!r}")
    return I


def competitive(S, vmax, km, I, ki):
    """Competitive inhibition: apparent Km increases, Vmax unchanged.

    ``v = Vmax·S / (Km·(1 + I/Ki) + S)``

    Mechanically: inhibitor binds free enzyme only, so saturating substrate
    outcompetes it and Vmax is recovered at infinite S.

    Examples
    --------
    >>> round(competitive(25.0, 10.0, 25.0, 0.0, 8.0), 6)     # no inhibitor
    5.0
    >>> round(competitive(25.0, 10.0, 25.0, 8.0, 8.0), 4)     # I = Ki -> Km doubles
    3.3333
    """
    S = np.asarray(S, dtype=float)
    I = _check_inhibitor(I, ki)
    return _out(vmax * S / (km * (1.0 + I / ki) + S))


def noncompetitive(S, vmax, km, I, ki):
    """Non-competitive inhibition: Vmax decreases, Km unchanged.

    ``v = Vmax·S / ((Km + S)·(1 + I/Ki))``

    Mechanically: inhibitor binds equally well to free enzyme and to the ES complex,
    so it removes a fraction of the enzyme from the reaction entirely.

    Examples
    --------
    >>> round(noncompetitive(25.0, 10.0, 25.0, 0.0, 12.0), 6)
    5.0
    >>> round(noncompetitive(25.0, 10.0, 25.0, 12.0, 12.0), 4)   # I = Ki -> Vmax halves
    2.5
    """
    S = np.asarray(S, dtype=float)
    I = _check_inhibitor(I, ki)
    return _out(vmax * S / ((km + S) * (1.0 + I / ki)))


def uncompetitive(S, vmax, km, I, ki):
    """Uncompetitive inhibition: Vmax and Km both decrease by the same factor.

    ``v = Vmax·S / (Km + S·(1 + I/Ki))``

    Mechanically: inhibitor binds only the ES complex. The Km *appears* to fall —
    which is why an uncompetitive inhibitor is easy to mistake for an activator if
    you only look at Km.

    Examples
    --------
    >>> round(uncompetitive(25.0, 10.0, 25.0, 0.0, 20.0), 6)
    5.0
    >>> round(uncompetitive(25.0, 10.0, 25.0, 20.0, 20.0), 4)   # I = Ki
    3.3333
    >>> round(uncompetitive(1e6, 10.0, 25.0, 20.0, 20.0), 3)     # saturating [S] -> Vmax/2
    5.0
    """
    S = np.asarray(S, dtype=float)
    I = _check_inhibitor(I, ki)
    return _out(vmax * S / (km + S * (1.0 + I / ki)))


def mixed_inhibition(S, vmax, km, I, ki, alpha=1.0):
    """Mixed inhibition: the general case that contains the others.

    ``v = Vmax·S / (Km·(1 + I/Ki) + S·(1 + I/(α·Ki)))``

    The two limiting cases are the ones that matter:

    * ``α = 1`` → non-competitive (equal affinity for E and ES);
    * ``α → ∞`` → competitive (no binding to ES).

    Fitting ``α`` is *not* the same as demonstrating mixed inhibition: with four
    parameters and a handful of concentrations, ``α`` is usually unidentifiable, and
    a confidence interval that spans two orders of magnitude is a warning, not a
    result. Module 6 covers this properly; here the function exists so the
    statement can be tested.
    """
    S = np.asarray(S, dtype=float)
    I = _check_inhibitor(I, ki)
    if alpha <= 0:
        raise ValueError(f"alpha must be positive, got {alpha!r}")
    return _out(vmax * S / (km * (1.0 + I / ki) + S * (1.0 + I / (alpha * ki))))


#: Name → (callable, parameter names, parameter starting guesses).
#: The starting guesses are data-independent because a fit whose starting point is
#: tuned per dataset cannot be automated; convergence problems are diagnosed
#: instead of hidden (see ``fit_nonlinear``).
INHIBITION_MODELS: dict[str, dict] = {
    "competitive": {
        "func": competitive,
        "params": ("vmax", "km", "ki"),
        "p0": (1.0, 1.0, 1.0),
        "description": "binds free enzyme only; apparent Km rises, Vmax unchanged",
    },
    "noncompetitive": {
        "func": noncompetitive,
        "params": ("vmax", "km", "ki"),
        "p0": (1.0, 1.0, 1.0),
        "description": "binds E and ES equally; Vmax falls, Km unchanged",
    },
    "uncompetitive": {
        "func": uncompetitive,
        "params": ("vmax", "km", "ki"),
        "p0": (1.0, 1.0, 1.0),
        "description": "binds ES only; both Vmax and apparent Km fall",
    },
}
