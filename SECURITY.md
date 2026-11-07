# Security policy

## What this repository is

`kinetics` is a teaching package for fitting enzyme kinetics data. It reads CSV files you give
it and does arithmetic. It does not open network connections, does not execute data, and has no
authentication surface. Its dependencies are `numpy` and `scipy`, plus `matplotlib`, `pyyaml`,
`nbformat` and `nbclient` for development and validation scripts.

The realistic risks are therefore not remote code execution but:

* **Malformed input** — a CSV with a huge number of columns, absurd values, or a non-numeric
  cell. These raise a `ValueError`/`InputError` with the offending line, or exit with code 2;
  they are not silently coerced.
* **Resource exhaustion** — `bootstrap_ci` with a very large `n_resamples` will run for a long
  time. It is a plain Python loop by design (see Module 4 §4.3); there is no hidden parallelism
  or thread pool to escape.
* **Supply chain** — a compromised release of `numpy` or `scipy` would affect any result here.
  Pin versions in your own environment if you are reproducing published numbers.

## Supported versions

| Version | Supported |
|---|---|
| 0.1.x | ✅ security fixes |
| < 0.1 | ❌ |

## Reporting a vulnerability

Open a **private security advisory** on the GitHub repository
(<https://github.com/cls-open/python-for-biochemists/security/advisories/new>) rather than a
public issue. Include the version (`kinetics --version`), a minimal input file that reproduces
the problem, and what you expected to happen.

Expect an acknowledgement within 7 days. Because this is an educational project maintained
voluntarily, please allow 30 days before public disclosure.

## Out of scope

* The correctness of a scientific conclusion drawn from your own data. That is what
  [`docs/limitations.md`](docs/limitations.md) is for, and it is a discussion, not a
  vulnerability.
* Anything requiring an attacker who can already modify the files under `datasets/` — at that
  point the data are the attack, and `check_data_metadata.py` is the defence: it fails if a
  committed CSV no longer matches the SHA-256 recorded in its metadata.

## A note on what "secure" means for a teaching tool

The most likely way this repository could harm a user is by teaching a habit that produces a
wrong number that looks right. That is why the guards are loud, the `nan` values are explained
rather than hidden, and the limitations page is longer than the feature list. If you find a case
where the package is quietly wrong, that is a security report as far as this project is
concerned — please send it.
