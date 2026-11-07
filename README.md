# Python for Biochemical Data and Enzyme Kinetics — CLS-Open #2

**CLS-Open repository 2 of 15** · Level: **intermediate** · Status: **v0.1.0 — complete**

> Most kinetics teaching stops at "plot it and read the slope". This course starts where that
> fails: the moment you have to decide *which* model, *which* weighting, *which* uncertainty,
> and what to write down when the data cannot answer the question.

The first CLS-Open repository with real numerical dependencies (`numpy`, `scipy`), and the
first where a wrong answer looks exactly like a right one.

## What you get in v0.1.0

| Component | Status |
|---|---|
| **Module 1** — Units, rates and what an assay measures | ✅ complete |
| **Module 2** — The Michaelis–Menten model and the two numbers people misquote | ✅ complete |
| **Module 3** — Linearisations and what the transformation costs | ✅ complete |
| **Module 4** — Uncertainty, standard errors and the bootstrap | ✅ complete |
| **Module 5** — Weighting, error models and the shape of noise | ✅ complete |
| **Module 6** — Choosing between models without fooling yourself | ✅ complete |
| **`kinetics`** — models, fitting, linearisations, inhibition, calibration, reporting, CLI | ✅ **206 tests, 100 % statement coverage**, 83 comparisons against independent SciPy implementations |
| 5 practicals with rubrics | ✅ complete |
| Quiz (Modules 1–6) with full answer key | ✅ complete |
| Instructor pack, two mini-projects, capstone specification | ✅ complete |
| Executed notebook (`01_fitting_decisions.ipynb`) | ✅ 8 code cells, outputs stored, re-executed by CI |
| 8 simulated datasets (115 rows) with SHA-256-verified metadata | ✅ complete |
| 4 figures, every displayed number recorded in `figures/figure_values.json` | ✅ reproduced by `make_figures.py --check` |
| Every citation DOI-verified against Crossref | ✅ 22/22 titles matched |

## Install

```bash
git clone https://github.com/cls-open/python-for-biochemists
cd python-for-biochemists
python3 -m venv .venv && source .venv/bin/activate    # Windows: use WSL2
python -m pip install -e ".[dev]"
kinetics --version        # kinetics 0.1.0
python -m pytest -q       # expect: 206 passed
```

`kinetics` needs only `numpy` and `scipy` at runtime.

## Run in five minutes

```bash
kinetics mm datasets/sample/michaelis_menten_mean.csv            # weighted fit + AICc
kinetics compare datasets/sample/michaelis_menten_mean.csv       # four methods side by side
kinetics inhibition datasets/sample/inhibition_series.csv        # mechanism by AICc
kinetics calibrate datasets/sample/calibration_mean.csv          # slope, LOD, LOQ
kinetics mm datasets/sample/malformed_example.csv; echo $?       # 2: the file is broken, loudly
jupyter lab notebooks/01_fitting_decisions.ipynb                 # executed walkthrough
python scripts/make_figures.py --check                           # figures match the data
```

Every subcommand speaks Markdown, JSON or self-contained HTML (`--format`, `--out`), and the
JSON is the same numbers as the Markdown — asserted by tests.

## The claims are measurements, not assertions

This repository's whole reason for existing:

* **The Lineweaver–Burk bias is demonstrated.** On `michaelis_menten_mean.csv`, where the
  nonlinear fit returns Vmax 9.98381 (−0.16 % from a simulated truth of 10.0), the double
  reciprocal returns **11.9869 (+19.9 %)**, and Km 31.9047 (+27.6 %). Hanes–Woolf and
  Eadie–Hofstee both land within 0.7 % — the damage is specific to the double reciprocal, not
  to linearisation in general.
* **The distribution, not the anecdote.** Over 400 simulated assays only **57 %** of
  Lineweaver–Burk Vmax estimates land within 20 % of truth, against 100 % for the nonlinear fit.
  Over 200 assays, exactly **5** recover Vmax within 1 % — which is why "my plot looked fine"
  is not evidence.
* **A weighting decision changes an LOD by a factor of 58.** The same seven calibration
  standards give LOD **0.109296** weighted by replicate SD and **6.357793** unweighted. Both
  numbers are "the LOD"; the σ that produced them is the result.
* **The negative result is included.** `datasets/sample/inhibition_weak.csv` is a real
  competitive inhibitor (Ki = 500 µM) measured with [I] that never exceeds 0.4 % of Ki. The fit
  returns "uncompetitive" as the best model with ΔAICc 0.48 — and the documentation teaches you
  to write "indistinguishable on these data" instead.
* **Nothing is silently fixed.** No row filtering, no unit conversion, no fallback optimiser.
  A negative rate is fitted but flagged in a note with the count and the worst value; a
  constant response raises; a linearised fit's AICc reads "not comparable" rather than printing
  a number that invites the wrong comparison.

## What is not yet present

Nothing here is an outline dressed as a course: if a module is listed above, it is written and
tested, and if it is not listed, it does not exist yet. The things this version deliberately does
**not** include are:

* **global fitting** (one Vmax/Km shared across datasets or mutants) — specified, not yet
  implemented; today it is a per-dataset fit, which is the wrong tool for that question;
* **substrate inhibition** as a first-class model;
* **robust regression** and errors-in-x — discussed in the limitations, not implemented;
* **any real lab data.** All eight datasets are synthetic, with seeds and error models published;
* **workflows beyond steady-state initial rates** — no progress curves, no stopped-flow.

[`docs/roadmap.md`](docs/roadmap.md) lists what is wanted in what order, and
[`docs/limitations.md`](docs/limitations.md) lists what is refused on purpose and why.

## The data are simulated, and say so

Every CSV carries `synthetic: true`, its generator seed, its error model
(`proportional_plus_absolute`, with the parameters) and a SHA-256 of the file, in
[`datasets/metadata/`](datasets/metadata). `scripts/make_datasets.py` regenerates all eight
files deterministically; `scripts/check_data_metadata.py` verifies the hashes. Nothing here is
passed off as wet-lab data, and the truth values are published in the metadata precisely so you
can check a fit against them.

## Repository layout

```
src/kinetics/      models, fitting, linearisations, inhibition, calibration, report, cli
tests/             206 tests; test_edges.py holds one test per user-visible guard
scripts/           make_datasets, make_figures, validate_against_scipy, check_* (data,
                   docs links, references, notebooks)
datasets/          sample CSVs + metadata YAML + index.json (hashes)
figures/           4 PNGs + figure_values.json (every number the modules quote)
notebooks/         01_fitting_decisions.ipynb (executed, re-run in CI)
docs/              6 modules, 5 practicals, quiz, instructor pack, validation, limitations,
                   references, roadmap
projects/          two mini-projects and a capstone
```

## Exit codes mean something

| Code | Meaning |
|---:|---|
| `0` | fitted and reported |
| `1` | the data were read but the fit refused: constant response, unidentifiable model, no inhibitor present |
| `2` | the file itself is unusable: missing column, non-numeric cell (the message names the line) |

## Documentation

Full course, figures and validation: **<https://cls-open.github.io/python-for-biochemists/>**

* [Validation](docs/validation.md) — exactly what has been checked, and against what
* [Limitations](docs/limitations.md) — what this cannot do, and why the omissions are deliberate
* [References](docs/references.md) — every citation, DOI-verified

## Licence

Content CC BY 4.0, code MIT — see [`LICENSES/`](LICENSES/README.md). Part of the
[CLS-Open](https://github.com/CLS-Open) initiative of fifteen interlinked courses in
computational life science.
