"""Rendering fit results as Markdown, JSON or a self-contained HTML page.

One data structure, three renderings — so a report cannot disagree with the numbers
it was generated from. The HTML carries inline CSS only: no external fonts, no
scripts, no images fetched over the network, because a report is something you email
to a collaborator who may open it on a locked-down machine.
"""

from __future__ import annotations

import json
import math
from html import escape

from .fitting import FitResult

__all__ = ["results_to_html", "results_to_json", "results_to_markdown"]


def _fmt(value, spec: str = ".4g") -> str:
    """Format a number, or an em dash where there is no number to print.

    ``nan`` is a value, but it is never an *estimate*: printing "nan" in a table invites
    a reader to think a number was computed and came out odd, when what happened is
    that no number exists.
    """
    if value is None:
        return "—"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return format(number, spec) if math.isfinite(number) else "—"


def _aicc_cell(result: FitResult) -> str:
    """AICc for comparable fits, an explicit refusal for transformed ones.

    A linearised fit is computed on a transformed response (1/v, v/S, S/v), so its
    residual sum of squares is in different units from a nonlinear fit's. Its AICc is
    therefore computable but *not comparable*, and printing it in the same column would
    make Lineweaver-Burk look like the best model on almost every dataset — which it is
    not. The CLI excludes these fits from its ranking for the same reason.
    """
    if result.method == "linearised":
        return "not comparable"
    value = result.aicc_or_none()
    return "—" if value is None else f"{value:.2f}"


def _md(text) -> str:
    """Escape a cell for a Markdown table.

    A pipe in a value (a dataset path, a note) silently adds a column and shifts
    every cell to its right. Escaping it keeps the table readable by both humans and
    parsers.
    """
    return str(text).replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ")


def _param_rows(results: list[FitResult]) -> list[tuple[str, dict[str, float], dict[str, float]]]:
    return [(r.model, r.params, r.stderr) for r in results]


def results_to_json(results: list[FitResult], *, context: dict | None = None) -> str:
    """Serialise one or more fits. ``None`` is preserved (it means "not applicable")."""
    payload = {
        "tool": "kinetics",
        "context": context or {},
        "results": [r.as_dict() for r in results],
    }
    return (
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)
        if _all_finite(payload)
        else json.dumps(_clean_nans(payload), indent=2, ensure_ascii=False)
    )


def _all_finite(obj) -> bool:
    if isinstance(obj, float):
        return math.isfinite(obj)
    if isinstance(obj, dict):
        return all(_all_finite(v) for v in obj.values())
    if isinstance(obj, list):
        return all(_all_finite(v) for v in obj)
    return True


def _clean_nans(obj):
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _clean_nans(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean_nans(v) for v in obj]
    return obj


def results_to_markdown(
    results: list[FitResult], *, title: str = "Kinetics fit", context: dict | None = None
) -> str:
    """Render fits as a Markdown report with a model-comparison table.

    Examples
    --------
    >>> from kinetics.fitting import FitResult
    >>> fit = FitResult(model="michaelis_menten", params={"vmax": 10.0, "km": 25.0},
    ...                 stderr={"vmax": 0.4, "km": 1.2}, n=7, k=2, rss=0.5,
    ...                 extra={"tss": 60.0}, notes=["Vmax is an extrapolation"])
    >>> print(results_to_markdown([fit], title="Demo", context={"dataset": "demo.csv"}))
    # Demo
    <BLANKLINE>
    | field | value |
    |---|---|
    | dataset | demo.csv |
    <BLANKLINE>
    ## Model comparison
    <BLANKLINE>
    | model | method | n | k | R² | AICc |
    |---|---|---:|---:|---:|---:|
    | michaelis_menten | nonlinear | 7 | 2 | 0.99167 | -11.47 |
    <BLANKLINE>
    ## Parameters
    <BLANKLINE>
    | model | parameter | estimate | std. error |
    |---|---|---:|---:|
    | michaelis_menten | vmax | 10 | 0.4 |
    | michaelis_menten | km | 25 | 1.2 |
    <BLANKLINE>
    ## Notes (kept with the numbers, not behind them)
    <BLANKLINE>
    - **michaelis_menten**: Vmax is an extrapolation
    <BLANKLINE>
    >>> "TODO" in results_to_markdown([fit])
    False
    """
    lines = [f"# {title}", ""]
    if context:
        lines += ["| field | value |", "|---|---|"]
        lines += [f"| {_md(k)} | {_md(v)} |" for k, v in context.items()]
        lines.append("")
    lines += [
        "## Model comparison",
        "",
        "| model | method | n | k | R² | AICc |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for r in results:
        r2 = "—" if not math.isfinite(r.r_squared) else f"{r.r_squared:.5f}"
        lines.append(f"| {r.model} | {r.method} | {r.n} | {r.k} | {r2} | {_aicc_cell(r)} |")
    lines += [
        "",
        "## Parameters",
        "",
        "| model | parameter | estimate | std. error |",
        "|---|---|---:|---:|",
    ]
    for model, params, errors in _param_rows(results):
        for name, value in params.items():
            lines.append(
                f"| {model} | {name} | {_fmt(value, '.6g')} | {_fmt(errors.get(name), '.3g')} |"
            )
    if any(r.method == "linearised" for r in results):
        lines += [
            "",
            (
                "AICc is reported only for fits on the untransformed response: a linearised fit "
                "is computed on 1/v, v/S or S/v, so its AICc is computable but not comparable "
                "with these."
            ),
        ]
    notes = [(r.model, note) for r in results for note in r.notes]
    if notes:
        lines += ["", "## Notes (kept with the numbers, not behind them)", ""]
        seen = set()
        for model, note in notes:
            if note in seen:
                continue
            seen.add(note)
            lines.append(f"- **{model}**: {note}")
    return "\n".join(lines) + "\n"


_CSS = """
:root { color-scheme: light dark; }
body { font: 16px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
       max-width: 62rem; margin: 0 auto; padding: 2rem 1rem; }
h1 { font-size: 1.5rem; margin-bottom: .2rem; }
h2 { font-size: 1.1rem; margin-top: 2rem; border-bottom: 1px solid #8884; padding-bottom: .2rem; }
table { border-collapse: collapse; width: 100%; font-size: .9rem; margin: .5rem 0; }
th, td { text-align: left; padding: .3rem .5rem; border-bottom: 1px solid #8883; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .9em; }
ul.notes li { margin-bottom: .4rem; }
.muted { color: #6b7280; font-size: .85rem; }
"""


def results_to_html(
    results: list[FitResult], *, title: str = "Kinetics fit", context: dict | None = None
) -> str:
    """Render a self-contained HTML page (inline CSS, no JavaScript, no remote assets)."""
    cmp_rows = []
    for r in results:
        r2 = "—" if not math.isfinite(r.r_squared) else f"{r.r_squared:.5f}"
        aicc = _aicc_cell(r)
        cmp_rows.append(
            "<tr>"
            f"<td><code>{escape(r.model)}</code></td><td>{escape(r.method)}</td>"
            f'<td class="num">{r.n}</td><td class="num">{r.k}</td>'
            f'<td class="num">{r2}</td><td class="num">{aicc}</td>'
            "</tr>"
        )
    param_rows = []
    for model, params, errors in _param_rows(results):
        for name, value in params.items():
            param_rows.append(
                "<tr>"
                f"<td>{escape(model)}</td><td><code>{escape(name)}</code></td>"
                f'<td class="num">{_fmt(value, ".6g")}</td>'
                f'<td class="num">{_fmt(errors.get(name), ".3g")}</td>'
                "</tr>"
            )
    notes = []
    seen = set()
    for r in results:
        for note in r.notes:
            if note in seen:
                continue
            seen.add(note)
            notes.append(f"<li><strong>{escape(r.model)}</strong>: {escape(note)}</li>")
    context_html = ""
    if context:
        context_html = (
            "<table><tbody>"
            + "".join(
                f"<tr><th>{escape(str(k))}</th><td>{escape(str(v))}</td></tr>"
                for k, v in context.items()
            )
            + "</tbody></table>"
        )
    notes_html = f'<h2>Notes</h2><ul class="notes">{"".join(notes)}</ul>' if notes else ""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title>
<style>{_CSS}</style>
</head>
<body>
<h1>{escape(title)}</h1>
<p class="muted">Generated by <code>kinetics</code> 0.1.0 from the data named below.
Model comparison by AICc; parameters with asymptotic standard errors.</p>
{context_html}
<h2>Model comparison</h2>
<table>
<thead><tr><th>model</th><th>method</th><th class="num">n</th><th class="num">k</th>
<th class="num">R²</th><th class="num">AICc</th></tr></thead>
<tbody>{"".join(cmp_rows)}</tbody>
</table>
<h2>Parameters</h2>
<table>
<thead><tr><th>model</th><th>parameter</th><th class="num">estimate</th>
<th class="num">std. error</th></tr></thead>
<tbody>{"".join(param_rows)}</tbody>
</table>
{notes_html}
<p class="muted">A lower AICc is evidence for a model, not proof of a mechanism.
Linearised fits are shown for comparison but their AICc is not comparable with a
nonlinear fit's (different residual scale). Asymptotic standard errors are optimistic
for nonlinear fits with small n — use the bootstrap intervals where they matter.</p>
</body>
</html>
"""
