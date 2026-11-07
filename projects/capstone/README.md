# Capstone (20–25 hours)

**Question.** Given a substrate series, an inhibitor series and a calibration curve, produce a
single report that someone else could audit — and be explicit about the parts the data cannot
support.

The capstone is deliberately not "fit four datasets". It is one analysis, done end to end, with
the decisions recorded as they are made. Most of the marks are for the reasoning, not the
numbers.

---

## What you start with

| File | What it contains |
|---|---|
| `datasets/sample/michaelis_menten.csv` | 9 substrate concentrations × 3 replicates, known truth |
| `datasets/sample/inhibition_series.csv` | 5 [S] × 4 [I], competitive inhibitor, Ki = 8 µM |
| `datasets/sample/calibration_standards.csv` | 7 standards × 3 replicates, no `sd` column |
| `datasets/sample/calibration_mean.csv` | the same standards as means, with replicate SD |

You may use any of the other committed datasets, but you must say which you used and why.

## The five deliverables

### 1. A analysis script (or notebook) that runs end to end

One entry point, no manual steps, no numbers typed by hand. It must write:

* a machine-readable results file (JSON) with every number the report quotes;
* a figure.

`scripts/make_figures.py` in this repository is a worked example of the pattern: the figure and
the JSON are produced by the same run, so they cannot disagree.

### 2. A methods section (300–400 words)

It must state, for each analysis: the model as a formula, the fitting method, the weighting and
why, the number of points and replicates, and how uncertainty was obtained (asymptotic or
bootstrap, and with which seed). Vague phrases — "standard nonlinear regression", "appropriate
weights" — lose marks.

### 3. A results section with the honest numbers

For each fit: the estimate, its standard error **or** confidence interval (state which), and the
comparison statistics that are actually comparable. Where a number is not identified, say so
instead of printing it.

### 4. One paragraph the data do not support, and you refuse to write

Choose at least one claim implied by the data that you decline to make — a mechanism, an LOD
without its σ, a Vmax quoted to four figures — and explain the refusal in three or four
sentences. This is a required deliverable and it is worth a fifth of the marks.

### 5. A reproducibility note

The exact commands, the Python and package versions, the seeds, and how long it takes. A reader
with the repository should be able to reproduce every number in your report in one command.

## What "good" looks like

A strong submission reads like a paper's methods and results written by someone who does not
trust their own numbers yet. It says *"the competitive model is better supported by ΔAICc =
42.55; the noncompetitive fit is retained in the table because the ranking alone does not
establish the mechanism"*. A weak submission reports the same numbers as conclusions.

Two specific habits are marked:

* **Comparability.** AICc is reported only between fits to the same response with the same
  weighting. If you put a linearised fit's AICc next to a nonlinear one, you lose marks — the
  package prints "not comparable" where this applies, and the report module does it for a
  reason.
* **The σ definition.** Every detection limit comes with the σ that produced it, in signal
  units, and the sentence that says where the σ came from.

## Marking

| Criterion | Points | What earns full credit |
|---|---:|---|
| End-to-end script; every number in the report generated, none typed | 20 | One command reproduces the report and the JSON matches it |
| Methods section: model, method, weighting, n, uncertainty, seed | 20 | A reader could re-run it without asking a question |
| Results: estimates with defensible uncertainty; comparability respected | 20 | No AICc across different responses or weightings |
| The refusal paragraph | 20 | A specific claim, correctly identified as unsupported, with the evidence |
| Reproducibility note | 10 | Versions, seeds, runtime, and the one command |
| Code quality: guards, no silent filtering, error messages a reader can act on | 10 | Matches the standard of `src/kinetics/` |
| **Total** | **100** |

A report in which any number cannot be traced to the submitted script is capped at 50.

## Time budget

| Phase | Hours |
|---|---:|
| Read the modules you skipped; set up the data | 3 |
| Fitting and model comparison, with wrong turns | 8 |
| Calibration and LOD, including the weighting decision | 4 |
| Figures and the JSON writer | 3 |
| Writing, and the refusal paragraph (usually the hardest) | 5 |
| Reproducibility check: delete the outputs, re-run, compare | 2 |

The last row is not optional. Deleting your outputs and regenerating them is the only way to
know the report is a product of the code rather than of your afternoon.
