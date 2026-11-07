"""kinetics — enzyme kinetics models, fitting and diagnostics for teaching.

The point of this package is not to fit curves. Fitting curves is three lines of
SciPy. The point is to make the *decisions inside a fit* visible:

* which model, and why (competitive vs non-competitive inhibition are different
  mechanisms, not different curves);
* which noise model (ordinary least squares assumes constant variance; enzyme
  assays usually have proportional error);
* which estimate of uncertainty (asymptotic standard errors from the fit are
  optimistic; bootstrap intervals are honest about small *n*);
* and which transformations are *wrong but popular* (Lineweaver–Burk linearisation
  distorts errors — this package refits it deliberately so learners can see the
  bias rather than be told about it).

Everything is deterministic: no fit depends on a random seed except the bootstrap,
and that seed is explicit and recorded in the result.

Units convention
----------------
Every function documents its units and none of them guesses. The package never
adds a number with units to a number without units, and a fit result records the
units it was given so that a report cannot silently mix µM with mM. See module 1.
"""

from __future__ import annotations

__version__ = "0.1.0"

__all__ = [  # noqa: RUF022 - grouped by module, not alphabetised
    "__version__",
    # models
    "michaelis_menten",
    "hill",
    "competitive",
    "noncompetitive",
    "uncompetitive",
    "mixed_inhibition",
    "INHIBITION_MODELS",
    # fitting
    "FitResult",
    "aicc",
    "fit_nonlinear",
    "fit_michaelis_menten",
    "fit_lineweaver_burk",
    "fit_eadie_hofstee",
    "fit_hanes",
    "compare_linearisations",
    "bootstrap_ci",
    # inhibition
    "fit_inhibition_series",
    "compare_inhibition_models",
    # calibration
    "CalibrationResult",
    "fit_calibration",
    "lod_loq",
    # report
    "results_to_markdown",
    "results_to_json",
    "results_to_html",
]

from .calibration import CalibrationResult, fit_calibration, lod_loq
from .fitting import (
    FitResult,
    aicc,
    bootstrap_ci,
    compare_linearisations,
    fit_eadie_hofstee,
    fit_hanes,
    fit_lineweaver_burk,
    fit_michaelis_menten,
    fit_nonlinear,
)
from .inhibition import compare_inhibition_models, fit_inhibition_series
from .models import (
    INHIBITION_MODELS,
    competitive,
    hill,
    michaelis_menten,
    mixed_inhibition,
    noncompetitive,
    uncompetitive,
)
from .report import results_to_html, results_to_json, results_to_markdown
