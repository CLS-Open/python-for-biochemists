"""One set of numbers, three renderings — and the renderings must agree.

The point of the report module is that a report cannot drift from the fit that
produced it. These tests read the three formats and check that they quote the same
estimates, because "the PDF said 9.6 and the JSON said 9.65" is the failure mode
this design exists to prevent.
"""

from __future__ import annotations

import json
import re

import pytest

from kinetics.fitting import FitResult, fit_michaelis_menten
from kinetics.report import results_to_html, results_to_json, results_to_markdown


@pytest.fixture
def fits(mm_data):
    S, v = mm_data(noise=0.2, seed=7)
    return [fit_michaelis_menten(S, v)]


def test_markdown_contains_every_parameter_estimate(fits):
    text = results_to_markdown(fits, title="T")
    for name, value in fits[0].params.items():
        assert name in text
        assert f"{value:.6g}" in text


def test_json_round_trips_to_the_same_numbers(fits):
    payload = json.loads(results_to_json(fits))
    assert payload["results"][0]["params"]["vmax"] == pytest.approx(fits[0].params["vmax"])
    assert payload["results"][0]["aicc"] == pytest.approx(round(fits[0].aicc, 4))


def test_html_quotes_the_same_estimate_as_json(fits):
    html = results_to_html(fits, title="T")
    payload = json.loads(results_to_json(fits))
    vmax = payload["results"][0]["params"]["vmax"]
    assert f"{vmax:.6g}" in html


def test_all_three_carry_the_notes(fits):
    note = fits[0].notes[0]
    assert note in results_to_markdown(fits)
    assert note in results_to_html(fits)
    assert note in json.loads(results_to_json(fits))["results"][0]["notes"]


def test_html_is_self_contained():
    """No remote assets: the report has to open on a machine with no network."""
    html = results_to_html([], title="Empty")
    assert "http://" not in html and "https://" not in html
    assert "<script" not in html.lower()
    assert "<link" not in html.lower()
    assert "<style>" in html


def test_html_escapes_untrusted_text():
    """A dataset name comes from a file path, so it is not trusted."""
    html = results_to_html([], context={"dataset": '<script>alert("x")</script>'})
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_markdown_escapes_table_breakers_from_context():
    text = results_to_markdown([], context={"note": "a|b"})
    assert "a|b" not in text


def test_undefined_aicc_renders_as_a_dash_not_a_crash():
    tiny = FitResult(
        model="m", params={"a": 1.0}, stderr={"a": 1.0}, n=3, k=2, rss=1.0, extra={"tss": 2.0}
    )
    for text in (results_to_markdown([tiny]), results_to_html([tiny]), results_to_json([tiny])):
        assert text.strip()
    assert "—" in results_to_markdown([tiny])
    assert json.loads(results_to_json([tiny]))["results"][0]["aicc"] is None


def test_nan_parameters_become_null_in_json():
    broken = FitResult(model="m", params={"a": float("nan")}, n=5, k=1, rss=1.0, extra={"tss": 2.0})
    payload = json.loads(results_to_json([broken]))  # must not raise
    assert payload["results"][0]["params"]["a"] is None


def test_markdown_table_has_one_row_per_parameter(fits):
    text = results_to_markdown(fits)
    rows = [
        line
        for line in text.splitlines()
        if line.startswith("| michaelis_menten |") and line.count("|") == 5
    ]
    assert len(rows) == len(fits[0].params)


def test_titles_and_context_are_rendered_both_ways(fits):
    md = results_to_markdown(fits, title="Assay 12", context={"protein": "3.2 mg"})
    html = results_to_html(fits, title="Assay 12", context={"protein": "3.2 mg"})
    assert "# Assay 12" in md
    assert "<title>Assay 12</title>" in html
    assert "3.2 mg" in md and "3.2 mg" in html


def test_no_placeholder_text_survives(fits):
    blob = results_to_markdown(fits) + results_to_html(fits) + results_to_json(fits)
    for marker in ("TODO", "FIXME", "XXX", "lorem"):
        assert marker not in blob
    assert not re.search(r"\bTBD\b", blob)


def test_linearised_fits_do_not_get_a_comparable_aicc(fits):
    """A transformed fit's AICc is computable but meaningless next to a nonlinear one."""
    from kinetics.fitting import compare_linearisations

    linearised = results_to_markdown(
        compare_linearisations([1.0, 2.0, 5.0, 10.0, 20.0], [1.1, 1.9, 5.2, 9.8, 20.5])
    )
    assert "not comparable" in linearised
    assert "not comparable with these" in linearised
    # ... and the nonlinear fit in the same table still has a number
    assert "| michaelis_menten | nonlinear |" in linearised


def test_non_finite_parameters_render_as_a_dash_not_nan(fits):
    from kinetics.fitting import FitResult

    broken = FitResult(
        model="m",
        params={"a": float("nan")},
        stderr={"a": float("nan")},
        n=5,
        k=1,
        rss=1.0,
        extra={"tss": 2.0},
    )
    text = results_to_markdown([broken])
    assert "nan" not in text.lower()
    assert "—" in text
