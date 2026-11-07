# Installation

Two routes: use the package as a tool (five minutes), or work on it as a project
(fifteen minutes). Both end with the same check.

## Requirements

| Component | Minimum | Tested on |
|---|---|---|
| Python | 3.10 | 3.13.16 |
| numpy | 1.24 | 2.3.5 |
| scipy | 1.10 | 1.17.1 |
| matplotlib (figures, optional) | 3.8 | 3.10 |
| pytest (tests, optional) | 8 | 9.0 |

No compiler is needed: `numpy` and `scipy` ship wheels for every supported platform.

## Route A — install and run the CLI

```bash
git clone https://github.com/cls-open/python-for-biochemists.git
cd python-for-biochemists
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e .
kinetics --version
```

Expected output:

```
kinetics 0.1.0
```

That one line verifies the package, the entry point and the dependency chain at once. If
it prints a traceback, the usual causes are a Python older than 3.10 or an `numpy` too
old for the SciPy build — read the error rather than upgrading everything.

## Route B — the full development environment

```bash
pip install -e ".[dev]"
python -m pytest -q
python scripts/check_data_metadata.py
```

The test suite is fast (about ten seconds on a laptop) and needs no network. The reference
checker does need network access; it is the only script that does:

```bash
python scripts/check_references.py
```

## Your first fit, in one command

```bash
kinetics mm datasets/sample/michaelis_menten_mean.csv
```

Expected output (abridged):

```
# Michaelis-Menten fit

| field | value |
|---|---|
| dataset | datasets/sample/michaelis_menten_mean.csv |
| points | 9 |
| weighting | 1/sd² from the sd column |

| model | method | n | k | R² | AICc |
|---|---|---:|---:|---:|---:|
| michaelis_menten | weighted_nonlinear | 9 | 2 | 0.99873 | -34.25 |
```

The closing digits may differ if you have a different SciPy version; the parameters will
not drift by more than the printed standard error. If a number here disagrees with the
one your own run prints, that is a documentation bug — report it.

## Exit codes

The CLI is designed to be used in a shell script, so its exit codes mean something:

| Code | Meaning | What to do |
|---|---|---|
| 0 | an answer was produced | read it |
| 1 | the input was readable but no answer can be given from it (too few points, no convergence) | look at stderr; the reason is printed there |
| 2 | the input could not be read or has the wrong columns | fix the file |

```bash
kinetics mm datasets/sample/malformed_example.csv; echo "exit=$?"
# error: datasets/sample/malformed_example.csv line 3: v='3.4x' is not a number
# exit=2
```

## Accepted column names

Concentrations and rates get written many ways. Rather than make users edit data files,
the reader accepts aliases (case-insensitive):

| Meaning | Accepted names |
|---|---|
| Substrate concentration | `S`, `[S]`, `substrate`, `s_um` |
| Rate | `v`, `V`, `rate`, `v0`, `velocity` |
| Inhibitor | `I`, `[I]`, `inhibitor` |
| Standard concentration | `concentration`, `conc`, `c` |
| Signal | `signal`, `response`, `absorbance`, `a` |
| Standard deviation | `sd`, `replicate_sd` |

Units are **not** converted and never inferred: whatever unit your file uses is the unit of
the answer. See [Module 1](modules/01-units-rates-and-what-an-assay-measures.md).

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `kinetics: command not found` | the venv is not active, or `pip install -e .` was run under a different interpreter | `python -c "import kinetics; print(kinetics.__file__)"` shows which one |
| `ValueError: fit did not converge for michaelis_menten` | fewer than four points, or data that do not increase with [S] | plot the data; the fit needs a monotone substrate series |
| `ValueError: need more data points than parameters` | two points for a two-parameter model | collect more data — this is not a software problem |
| `InputError: ... missing column(s)` | header name not recognised | check the alias table above |
| Figures look different from the ones in `docs/` | different matplotlib version, not different data | the numbers in `figures/figure_values.json` are what matter |

## Uninstalling

```bash
pip uninstall kinetics
```

The package writes nothing outside its own installation directory. It never modifies your
data files. `--out` is the only way it writes anything, and it will overwrite an existing
file without warning — that is a deliberate choice against the alternative of `--out`
failing silently on a path you had already created.
