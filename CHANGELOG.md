# Changelog

All notable changes to this repository are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Because the documentation quotes numbers computed from the committed datasets, a change to a
dataset or a generator is a change to the documentation and is versioned as one.

## [0.1.0] — 2026-11-07

First complete release. Six modules, five practicals, quiz, instructor pack, two mini-projects
and a capstone, with a tested package behind every number.

### Added

* **`kinetics` package** — `models` (Michaelis–Menten, Hill, the three inhibition mechanisms,
  mixed inhibition as a forward model), `fitting` (nonlinear least squares with optional
  weights, four linearisations, AICc, percentile bootstrap), `inhibition` (mechanism comparison
  by AICc), `calibration` (weighted and unweighted linear calibration, LOD and LOQ),
  `report` (Markdown, JSON, self-contained HTML) and `cli` (`mm`, `compare`, `inhibition`,
  `calibrate`; exit codes 0/1/2).
* **206 tests**, 100 % statement coverage, 18 doctests that also pass under
  `PYTHONHASHSEED=1..4`.
* **83 comparisons** against independent SciPy/NumPy implementations
  (`scripts/validate_against_scipy.py`): 81 deterministic to 5.0e-07 relative, 2 Monte-Carlo
  bootstrap comparisons at 4.7e-02, which is the agreement a correct resampling loop should
  show.
* **Eight simulated datasets**, 115 rows, each with metadata carrying `synthetic: true`, the
  generating seed, the error model and a SHA-256 hash (`scripts/check_data_metadata.py`
  verifies all eight).
* **Four figures** whose every printed number is recorded in `figures/figure_values.json`;
  `scripts/make_figures.py --check` rebuilds them into a temporary directory and fails if the
  PNGs or the JSON have drifted.
* **Executed notebook** `01_fitting_decisions.ipynb` (8 code cells), re-executed in CI from both
  the repository root and `/tmp`.
* **Documentation**: six modules, five practicals with additive rubrics, a 16-question quiz with
  a full answer key, an instructor pack, validation, limitations, references and a roadmap.
* **22 DOI-verified references** (`scripts/check_references.py`), each checked against the title
  Crossref returns — a DOI that resolves to the wrong paper is treated as a defect.
* **CI**: tests on CPython 3.10–3.13, lint, doctests per hash seed, notebook execution, figure
  reproduction, dataset hashes, documentation links, strict `mkdocs` build, SciPy
  cross-validation and weekly DOI resolution.

### Design decisions worth knowing about

* **Non-positive rates are flagged, not refused.** An early version raised on any `v <= 0`; it
  made the package unable to fit its own `cooperative_hill` dataset, whose lowest-rate point
  sits below the detector floor. The fit now carries the count and the worst value in a note,
  because only a human can tell a detector floor from a data error.
* **A constant response raises.** `bootstrap_ci` on a constant `y` used to return zero-width
  intervals, which is a confident answer to an impossible question.
* **A linearised fit never gets a comparable AICc.** The report prints "not comparable" rather
  than a number, because the response scale differs.
* **`bootstrap_ci` requires a seed**, and records it in the result.

### Known limitations

See [`docs/limitations.md`](docs/limitations.md). The short version: no global fitting, no
substrate inhibition, no robust losses, no errors-in-x, and every dataset is simulated.
