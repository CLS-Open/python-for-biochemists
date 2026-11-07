# Contributing

This repository is part of the CLS-Open initiative: fifteen interlinked courses that share one
standard, enforced by [`cls-toolkit`](https://github.com/cls-open/cls-toolkit). The standard is
narrow on purpose, and the two rules below cover most of it.

## Rule 1 — every number in the documentation comes from executed code

If you change a number in a module, a practical or a README, include the command that produced
it. If a claim is new, add the test or script that asserts it. The repository already contains
the machinery:

| Change | What must pass |
|---|---|
| a claim about the Lineweaver–Burk bias | `python -m pytest tests/test_linearisation.py -q` |
| a change to a dataset or generator | `python scripts/check_data_metadata.py` |
| a change to a figure or the numbers it prints | `python scripts/make_figures.py --check` |
| a change to a notebook's code | `python scripts/check_notebooks.py` |
| a new or changed citation | `python scripts/check_references.py` (needs network) |
| any documentation edit | `python scripts/check_docs_links.py` and `mkdocs build --strict` |

Run the whole set before opening a pull request:

```bash
python -m pytest -q
python -m pytest -q --doctest-modules src/kinetics
ruff check . && ruff format --check .
python scripts/validate_against_scipy.py
python scripts/check_data_metadata.py
python scripts/make_figures.py --check
python scripts/check_notebooks.py
python scripts/check_docs_links.py
mkdocs build --strict
```

## Rule 2 — nothing is silently fixed

No silent row filtering, no silent unit conversion, no silent fallback optimiser, no injectable
defaults that hide a decision. A guard that refuses or a note that warns is preferred to a
correction the user never sees. If you add validation, add a test in `tests/test_edges.py` —
that file exists so every user-visible message is asserted somewhere.

## Datasets

All data in this repository are synthetic and must stay that way unless a pull request also adds
the provenance a real dataset needs (source, retrieval date, licence, and `synthetic: false` in
the metadata). New datasets must be produced by `scripts/make_datasets.py` — deterministic,
seeded — and must carry `true_parameters` and an `error_model` in their metadata.

## Citations

A DOI that resolves to the wrong paper is worse than a dead link, because it looks correct.
`check_references.py` compares the title printed in `docs/references.md` with the title Crossref
returns for that DOI. If the check fails, the reference is wrong even if the DOI is valid.

## Style

* British spelling in prose (behaviour, modelling, linearisation), `µ` for micro, `·` for
  multiplication in units.
* Line length 100, `ruff` for lint and format; the configuration is in `pyproject.toml` and the
  ignores are documented there.
* Prose should say what is **not** covered. Modules end with a "what this module did not cover"
  list; keep it honest when you add material.
* Commit messages: an imperative subject line, then why — not what. The diff already says what.

## Reporting a problem

Documentation errors are treated as defects, not as cosmetic issues: a wrong number teaches a
wrong habit. Open an issue with the command you ran and the output you got. If you believe a
guard is wrong — a fit that refuses something it should accept — say so; the design decisions in
[`docs/limitations.md`](docs/limitations.md) are recorded precisely so they can be argued with.
