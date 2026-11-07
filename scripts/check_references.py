#!/usr/bin/env python3
"""Verify every cited DOI: does it resolve, and does its title match what we claim?

A fabricated reference is the most damaging kind of error an educational repository can
contain, because a learner will repeat it. This script makes that failure mode
mechanical: it reads every DOI out of ``docs/references.md`` (and the module reference
tables), asks Crossref, and compares the returned title with the title we printed.

Matching is deliberately fuzzy rather than exact: Crossref titles contain typos,
diacritics, and historic variants ("rcommendations" for "recommendations" is a real
example from the NC-IUB record). A normalised token overlap above the threshold counts as
a match; anything below it is reported for a human to look at.

Requires network. CI runs it on push to ``main`` and on a schedule, not on every PR from
a fork, so that an offline contributor is not blocked by Crossref's uptime.

    python scripts/check_references.py
    python scripts/check_references.py --offline     # only check internal consistency
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata
import urllib.error
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SOURCES = [REPO / "docs" / "references.md", *sorted((REPO / "docs" / "modules").glob("*.md"))]
USER_AGENT = "cls-open-python-for-biologists/0.1 (reference checker; mailto:noreply@cls-open.org)"
MIN_OVERLAP = 0.55

#: `[10.1234/abc.def](https://doi.org/10.1234/abc.def)` inside a table row, or a bare DOI.
#: Parentheses are legitimate inside a DOI -- `10.1016/0006-3002(63)91800-6` is a real one,
#: and an earlier version of this pattern stopped at the first `)`, reported a 404, and made
#: two correct references look fabricated. Balanced parentheses are kept; a trailing `)` that
#: closes a markdown link is stripped in :func:`clean_doi`.
ROW_RE = re.compile(r"^.*?\*\*(?P<title>.+?)\*\*.*?(?P<doi>10\.\d{4,9}/[^\s\]<>]+)", re.MULTILINE)


def clean_doi(raw: str) -> str:
    """Strip markdown/punctuation from the end of a DOI without eating balanced brackets."""
    doi = raw.rstrip(".,;:")
    while doi.endswith(")") and doi.count("(") < doi.count(")"):
        doi = doi[:-1]
    return doi


def normalise(text: str) -> set[str]:
    text = unicodedata.normalize("NFKD", text)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    return {w for w in re.findall(r"[a-z0-9]+", text.lower()) if len(w) > 2}


def overlap(claimed: str, actual: str) -> float:
    a, b = normalise(claimed), normalise(actual)
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def collect() -> list[tuple[Path, str, str]]:
    found: list[tuple[Path, str, str]] = []
    for path in SOURCES:
        if not path.is_file():
            continue
        for match in ROW_RE.finditer(path.read_text(encoding="utf-8")):
            title = match.group("title").strip()
            doi = clean_doi(match.group("doi"))
            if doi.lower().startswith("10."):
                found.append((path, title, doi))
    # de-duplicate by DOI, keeping first occurrence
    seen: dict[str, tuple[Path, str, str]] = {}
    for path, title, doi in found:
        seen.setdefault(doi.lower(), (path, title, doi))
    return list(seen.values())


def fetch(doi: str) -> dict:
    url = f"https://api.crossref.org/works/{doi}"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=30) as fh:
        return json.load(fh)["message"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="check structure only; no network")
    args = ap.parse_args()

    refs = collect()
    if not refs:
        print(
            "error: no references found — has docs/references.md changed format?", file=sys.stderr
        )
        return 2
    print(f"found {len(refs)} cited DOI(s) in {len(SOURCES)} file(s)")

    if args.offline:
        for _path, title, doi in refs:
            print(f"  {doi:45} {title[:58]}")
        print("\nOffline: resolution not checked.")
        return 0

    failures = 0
    for path, title, doi in refs:
        rel = path.relative_to(REPO)
        try:
            message = fetch(doi)
        except urllib.error.HTTPError as exc:
            print(f"FAIL {doi}: HTTP {exc.code}")
            failures += 1
            continue
        except (urllib.error.URLError, TimeoutError) as exc:
            print(f"SKIP {doi}: network error ({exc})")
            continue
        actual = (message.get("title") or ["(no title)"])[0]
        score = overlap(title, actual)
        if score >= MIN_OVERLAP:
            print(f"OK   {doi}  overlap={score:.2f}  {actual[:60]}")
        else:
            print(
                f"FAIL {doi}  claimed {title[:50]!r} but Crossref says {actual[:50]!r} (overlap {score:.2f})"
            )
            print(f"     in {rel}")
            failures += 1

    if failures:
        print(f"\n{failures} reference(s) did not verify — fix or remove before publishing")
        return 1
    print(f"\nOK: all {len(refs)} cited DOIs resolve and match their claimed titles.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
