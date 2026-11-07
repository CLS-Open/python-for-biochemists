"""Command-line interface.

Design rules, all of them deliberate:

* **stdout is the result, stderr is for humans.** ``--format json`` can be piped
  into ``jq`` without a warning line corrupting the document.
* **Exit codes mean something.** ``0`` an answer was produced, ``1`` the input was
  readable but no answer can be given from it (a fit that did not converge, too few
  points, an unidentifiable model), ``2`` the input could not be read or does not
  have the columns this subcommand needs. A ``1`` always comes with a report on
  stdout explaining what failed; a ``2`` comes with a message on stderr.
* **Every number in the output was computed by this run.** Nothing is copied from a
  tutorial or a remembered value.

Typical use::

    kinetics mm datasets/sample/michaelis_menten.csv
    kinetics compare datasets/sample/michaelis_menten.csv --format json
    kinetics inhibition datasets/sample/inhibition_series.csv
    kinetics calibrate datasets/sample/calibration_standards.csv

Input files are CSV with a header row and one measurement per row. Column names
are lower-cased and stripped; ``[S]``, ``S_uM`` and ``substrate`` are all accepted
for substrate concentration, because the alternative is users editing data files to
suit a parser.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from pathlib import Path

from . import __version__
from .calibration import fit_calibration
from .fitting import bootstrap_ci, compare_linearisations, fit_michaelis_menten, fit_nonlinear
from .inhibition import compare_inhibition_models
from .models import hill
from .report import results_to_html, results_to_json, results_to_markdown

EXIT_OK = 0
EXIT_NO_ANSWER = 1
EXIT_UNREADABLE = 2

_ALIASES = {
    "s": "S",
    "substrate": "S",
    "s_um": "S",
    "[s]": "S",
    "s_µm": "S",
    "v": "v",
    "rate": "v",
    "v0": "v",
    "velocity": "v",
    "v_umol_min_mg": "v",
    "i": "I",
    "inhibitor": "I",
    "[i]": "I",
    "concentration": "concentration",
    "conc": "concentration",
    "c": "concentration",
    "signal": "signal",
    "response": "signal",
    "absorbance": "signal",
    "a": "signal",
    "sd": "sd",
    "replicate_sd": "sd",
}


class InputError(Exception):
    """The file cannot be used by this subcommand at all (exit code 2)."""


def _norm(name: str) -> str:
    key = name.strip().lower()
    return _ALIASES.get(key, key)


def read_csv(path: str) -> list[dict[str, float]]:
    """Read a CSV of measurements into a list of dicts with numeric values.

    Raises
    ------
    InputError
        File missing, empty, without a header, or containing a non-numeric cell.
        The message names the row number so the user can look at it.
    """
    p = Path(path)
    if not p.exists():
        raise InputError(f"no such file: {path}")
    try:
        text = p.read_text(encoding="utf-8-sig")
    except OSError as exc:  # permissions, a directory given as a path, ...
        raise InputError(f"cannot read {path}: {exc}") from None
    rows = list(csv.reader(io.StringIO(text)))
    rows = [r for r in rows if any(cell.strip() for cell in r)]
    if not rows:
        raise InputError(f"{path} is empty")
    if len(rows) < 2:
        raise InputError(f"{path} has a header but no measurements")
    header = [_norm(c) for c in rows[0]]
    out: list[dict[str, float]] = []
    for lineno, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise InputError(
                f"{path} line {lineno}: {len(row)} fields but the header has {len(header)}"
            )
        record: dict[str, float] = {}
        for key, cell in zip(header, row, strict=True):
            cell = cell.strip()
            if cell == "":
                continue
            try:
                record[key] = float(cell)
            except ValueError:
                raise InputError(f"{path} line {lineno}: {key}={cell!r} is not a number") from None
        out.append(record)
    return out


def _require(rows: list[dict[str, float]], columns: tuple[str, ...], path: str) -> None:
    present = set().union(*(r.keys() for r in rows))
    missing = [c for c in columns if c not in present]
    if missing:
        raise InputError(
            f"{path} is missing column(s) {missing}; found {sorted(present)}. "
            f"This subcommand needs {list(columns)}; see docs/installation.md for the "
            "accepted column names."
        )


def _emit(text: str, out: str | None) -> None:
    if out:
        Path(out).write_text(text, encoding="utf-8")
        print(f"wrote {out}", file=sys.stderr)
    else:
        sys.stdout.write(text)


def _render(results, fmt: str, title: str, context: dict) -> str:
    if fmt == "json":
        return results_to_json(results, context=context) + "\n"
    if fmt == "html":
        return results_to_html(results, title=title, context=context)
    return results_to_markdown(results, title=title, context=context)


# ---------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------


def cmd_mm(args) -> int:
    rows = read_csv(args.dataset)
    _require(rows, ("S", "v"), args.dataset)
    S = [r["S"] for r in rows]
    v = [r["v"] for r in rows]
    sigma = [r["sd"] for r in rows] if all("sd" in r for r in rows) else None

    results = [fit_michaelis_menten(S, v, sigma=sigma)]
    context = {
        "dataset": args.dataset,
        "points": str(len(S)),
        "weighting": "1/sd² from the sd column" if sigma else "unweighted (ordinary least squares)",
    }
    if args.bootstrap:
        ci = bootstrap_ci(S, v, n_resamples=args.bootstrap, seed=args.seed)
        context["bootstrap"] = f"{args.bootstrap} resamples, seed {args.seed}"
        for name, (lo, hi) in ci.items():
            context[f"{name} 95% CI"] = f"{lo:.4g} to {hi:.4g}"
    if args.hill:
        results.append(
            fit_nonlinear(
                S,
                v,
                hill,
                (max(v), sorted(S)[len(S) // 2], 1.0),
                ("vmax", "k_half", "n"),
                name="hill",
                sigma=sigma,
                notes=[
                    (
                        "a Hill coefficient is a curve-shape parameter until a mechanism "
                        "is demonstrated; module 5 explains why, and why a fit to n = 1.0 "
                        "is not evidence against cooperativity"
                    )
                ],
            )
        )
    _emit(_render(results, args.format, "Michaelis-Menten fit", context), args.out)
    return EXIT_OK


def cmd_compare(args) -> int:
    rows = read_csv(args.dataset)
    _require(rows, ("S", "v"), args.dataset)
    S = [r["S"] for r in rows]
    v = [r["v"] for r in rows]

    # compare_linearisations already returns the nonlinear Michaelis-Menten fit
    # alongside the three linearisations, so it is not fitted twice.
    results = compare_linearisations(S, v)
    results.append(
        fit_nonlinear(
            S,
            v,
            hill,
            (max(v), sorted(S)[len(S) // 2], 1.0),
            ("vmax", "k_half", "n"),
            name="hill",
            notes=[
                ("the Hill coefficient is a shape parameter, not evidence of subunit cooperativity")
            ],
        )
    )

    # Rank by AICc where it is defined. A linearised fit is logged on a different
    # residual scale (1/v or S/v), so its AICc is NOT comparable to the nonlinear
    # fits' and it is excluded from the ranking rather than silently compared.
    nonlinear = [r for r in results if r.method != "linearised"]
    best = min(nonlinear, key=lambda r: r.aicc)
    deltas = {r.model: r.aicc - best.aicc for r in nonlinear}
    for r in nonlinear:
        r.notes = [
            *r.notes,
            f"ΔAICc = {deltas[r.model]:.2f} against the best model ({best.model})"
            + (
                "; Δ ≤ 2 means the data do not separate these two" if deltas[r.model] <= 2.0 else ""
            ),
        ]
    context = {
        "dataset": args.dataset,
        "points": str(len(S)),
        "AICc ranking": ", ".join(
            f"{m} (Δ={d:.2f})" for m, d in sorted(deltas.items(), key=lambda kv: kv[1])
        ),
        "excluded from ranking": (
            "linearised fits: different residual scale, so AICc is not defined comparably"
        ),
    }
    order = {"michaelis_menten": 0, "hill": 1, "hanes": 2, "eadie_hofstee": 3, "lineweaver_burk": 4}
    results = sorted(results, key=lambda r: (order.get(r.model, 9), r.method))
    _emit(_render(results, args.format, "Model comparison", context), args.out)
    return EXIT_OK


def cmd_inhibition(args) -> int:
    rows = read_csv(args.dataset)
    _require(rows, ("S", "I", "v"), args.dataset)
    comparison = compare_inhibition_models(rows)
    if args.format == "json":
        text = json.dumps({"dataset": args.dataset, **comparison.as_dict()}, indent=2) + "\n"
    else:
        results = comparison.results
        for r in results:
            r.notes = [*r.notes, f"ΔAICc = {comparison.delta_aicc[r.model]:.2f}"]
        text = _render(
            results,
            args.format,
            "Inhibition mechanism comparison",
            {
                "dataset": args.dataset,
                "best by AICc": comparison.best,
                "caveat": "a best-fitting mechanism is not a proven mechanism",
            },
        )
    _emit(text, args.out)
    return EXIT_OK


def cmd_calibrate(args) -> int:
    rows = read_csv(args.dataset)
    _require(rows, ("concentration", "signal"), args.dataset)
    conc = [r["concentration"] for r in rows]
    sig = [r["signal"] for r in rows]
    sd = [r["sd"] for r in rows] if all("sd" in r for r in rows) else None
    result = fit_calibration(conc, sig, weighted=bool(sd), replicate_sd=sd)
    if args.format == "json":
        text = json.dumps({"dataset": args.dataset, **result.as_dict()}, indent=2) + "\n"
    elif args.format == "html":
        lines = ["<h2>Calibration</h2><table>"]
        for key, value in result.as_dict().items():
            if key == "notes":
                continue
            lines.append(f"<tr><th>{key}</th><td>{value}</td></tr>")
        lines.append("</table>")
        lines += [f"<p>{note}</p>" for note in result.notes]
        text = "<!DOCTYPE html><meta charset='utf-8'><title>Calibration</title>\n" + "\n".join(
            lines
        )
    else:
        d = result.as_dict()
        lines = [
            "# Calibration",
            "",
            "| quantity | value |",
            "|---|---|",
            f"| slope | {d['slope']:.6g} |",
            f"| intercept | {d['intercept']:.6g} |",
            f"| slope std. error | {d['stderr_slope']:.4g} |",
            f"| intercept std. error | {d['stderr_intercept']:.4g} |",
            f"| R² | {d['r_squared']:.6f} |",
            f"| residual SD | {d['residual_sd']:.6g} |",
            f"| standards (n) | {d['n']} |",
            f"| weighted | {d['weighted']} |",
            f"| LOD (3.3σ/S) | {d['lod'] if d['lod'] is not None else 'not finite'} |",
            f"| LOQ (10σ/S) | {d['loq'] if d['loq'] is not None else 'not finite'} |",
            "",
        ]
        if result.notes:
            lines += ["## Notes", ""] + [f"- {note}" for note in result.notes]
        text = "\n".join(lines) + "\n"
    _emit(text, args.out)
    return EXIT_OK


# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="kinetics",
        description="Fit enzyme kinetics and calibration data, and show the working.",
    )
    p.add_argument("--version", action="version", version=f"kinetics {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("dataset", help="CSV file, one measurement per row")
        sp.add_argument(
            "--format",
            choices=("markdown", "json", "html"),
            default="markdown",
            help="output format (default: markdown)",
        )
        sp.add_argument("--out", help="write to this path instead of stdout")

    mm = sub.add_parser("mm", help="fit Vmax and Km (optionally with bootstrap intervals)")
    common(mm)
    mm.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        metavar="N",
        help="percentile bootstrap with N resamples (needs >= 4 points)",
    )
    mm.add_argument(
        "--seed",
        type=int,
        default=20261106,
        help="bootstrap seed, recorded in the output (default: 20261106)",
    )
    mm.add_argument("--hill", action="store_true", help="also fit the Hill equation")
    mm.set_defaults(func=cmd_mm)

    cmp_ = sub.add_parser(
        "compare", help="compare Michaelis-Menten, Hill and the three linearisations"
    )
    common(cmp_)
    cmp_.set_defaults(func=cmd_compare)

    inh = sub.add_parser("inhibition", help="compare competitive / non-competitive / uncompetitive")
    common(inh)
    inh.set_defaults(func=cmd_inhibition)

    cal = sub.add_parser("calibrate", help="linear calibration with LOD and LOQ")
    common(cal)
    cal.set_defaults(func=cmd_calibrate)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_UNREADABLE
    except ValueError as exc:
        # Readable input, no usable answer: say why on stderr, and keep exit 1.
        print(f"no answer: {exc}", file=sys.stderr)
        return EXIT_NO_ANSWER


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
