# Validation — exactly what has been checked, and against what

This page exists because "tested" is a claim, and a claim without a command is a slogan.
Everything below was run on this repository at the commit it ships from. Where a number
appears, it is the number the command printed.

## The short version

| Check | Command | Result |
|---|---|---|
| Unit tests | `python -m pytest -q` | **206 passed** |
| Statement coverage | `python -m pytest --cov=kinetics --cov-report=term-missing` | **100 %** (689 statements, 0 missed) |
| Doctests | `python -m pytest --doctest-modules src/kinetics` | **18 passed** |
| Doctests are order-independent | same, under `PYTHONHASHSEED=1..4` | **18 passed × 4 seeds** |
| Independent numerical check | `python scripts/validate_against_scipy.py` | **83 comparisons**, see below |
| Dataset integrity | `python scripts/check_data_metadata.py` | **8/8 datasets**, 115 rows, all SHA-256 verified |
| Figures reproduce | `python scripts/make_figures.py --check` | **4/4** PNGs and `figure_values.json` unchanged |
| Notebook executes | `python scripts/check_notebooks.py` | **1/1**, 8 code cells, no error |
| Citations resolve | `python scripts/check_references.py` | **22/22** DOIs, titles match Crossref |
| Documentation builds | `mkdocs build --strict` | **0 warnings** |
| Documentation links | `python scripts/check_docs_links.py` | **0 broken** internal links |
| Package builds and installs | `pip install -e .` in a clean environment | entry point `kinetics` works |

## What "100 % coverage" does and does not mean

Every statement in `src/kinetics/` is executed by the test suite. That is a floor, not a
guarantee: coverage cannot tell you the numbers are right, only that no line is dead code
nobody has ever run. The checks that address correctness are the next two sections.

The suite is deliberately adversarial. `tests/test_edges.py` holds one test per guard,
written so that every user-visible error message in the package is asserted somewhere. If a
guard is removed, a message changes, or a `nan` stops being explained, that file fails.

## Independent numerical validation

`scripts/validate_against_scipy.py` recomputes the package's results with separate
implementations — `scipy.optimize.curve_fit` directly, closed-form ordinary least squares for
the linearisations, the AICc definition written out by hand, an independent resampling loop
for the bootstrap, and closed-form calibration identities.

```
OK   Michaelis-Menten fit vs curve_fit           50 comparisons
OK   linearisations vs closed-form OLS            9 comparisons
OK   AICc vs definition                           8 comparisons
OK   bootstrap CI vs independent resampling       2 comparisons
OK   calibration and LOD vs closed form           7 comparisons
OK   inhibition fit vs curve_fit                  3 comparisons
OK   Hill(n=1) == Michaelis-Menten                4 comparisons

83 comparisons between kinetics and independent SciPy/NumPy implementations
  81 deterministic comparisons: worst relative difference 5.00e-07
  2 Monte-Carlo comparisons (bootstrap): worst relative difference 4.67e-02, expected to be
    percent-level because the random streams differ
```

The two Monte-Carlo comparisons are reported separately on purpose. A bootstrap that agreed
to 1e-7 with an independent implementation would mean one of the two was not resampling —
the agreement should be percent-level, and it is.

## The claims the documentation makes are tests

Several statements in the modules would be marketing if they were not executable. Each one
below is asserted, so the documentation cannot drift away from the code:

| Claim | Where it is asserted |
|---|---|
| Lineweaver–Burk fails catastrophically in ~30 % of assays at σ = 0.10, the nonlinear fit in none | `tests/test_linearisation.py` (300 seeded assays) |
| Hanes–Woolf is the least damaging linearisation | `tests/test_linearisation.py` |
| A linearised fit's information criterion is **not** comparable with a nonlinear one | `tests/test_report.py` (asserts the cell reads "not comparable") |
| A constant response raises rather than returning a zero-width interval | `tests/test_edges.py` |
| An over-parameterised fit reports `nan` standard errors and says why | `tests/test_edges.py` |
| Non-positive rates are flagged in a note, not silently dropped | `tests/test_edges.py` |
| `cooperative_hill` really is sigmoid (Hill n > Michaelis–Menten) | `tests/test_datasets.py` |
| `inhibition_weak` really is indistinguishable between mechanisms | `tests/test_inhibition.py` |
| The QC dataset still contains exactly its documented defects | `tests/test_datasets.py` + `scripts/check_data_metadata.py` |

## Determinism

* `scripts/make_datasets.py` regenerates all eight CSVs; `check_data_metadata.py` compares
  them against the SHA-256 recorded in each YAML file. A silent change to a generator — or to
  the digit of a data file — fails the check.
* `scripts/make_figures.py --check` rebuilds every figure and compares the PNG bytes and the
  `figure_values.json` payload.
* `bootstrap_ci` takes a **required** seed and records it in the result, so no interval in
  this repository is unreproducible.
* The notebook is re-executed by `check_notebooks.py`; committed outputs are replaced, never
  trusted.

## What has not been validated

Stated plainly, because this list is the one that matters:

* **No comparison against wet-lab data.** Every dataset here is simulated. The code is
  validated against independent implementations, not against a spectrophotometer.
* **No validation of the model, only of the fitter.** Agreement with `curve_fit` proves the
  arithmetic; it says nothing about whether Michaelis–Menten is the right model for your
  enzyme.
* **No guarantee of numerical stability outside the tested range.** Substrate concentrations
  spanning more than about six orders of magnitude, or rates within 1e-12 of zero, are outside
  what the tests exercise.
* **The interpreter matrix is 3.10–3.13.** Locally everything above was verified on CPython
  3.13.16; the older three versions are covered by the CI matrix, not by this page, so a
  version-specific failure would surface there rather than here. No wheels are built for other
  interpreters.
* **The reference checker needs the network.** It is skipped, not faked, when offline
  (`--offline` reports internal consistency only).

## How to run all of it

```bash
git clone https://github.com/cls-open/python-for-biochemists
cd python-for-biochemists
pip install -e ".[dev]"
python -m pytest -q
python -m pytest -q --doctest-modules src/kinetics
python scripts/validate_against_scipy.py
python scripts/check_data_metadata.py
python scripts/make_figures.py --check
python scripts/check_notebooks.py
python scripts/check_references.py          # needs network
python scripts/check_docs_links.py
mkdocs build --strict                       # needs mkdocs-material
```

If any of these prints something different from the tables above, the honest conclusion is
that this page is wrong — please open an issue with the command and the output. Documentation
errors are treated as defects in this repository, because a wrong number teaches a wrong
habit.
