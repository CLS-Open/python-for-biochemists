#!/usr/bin/env python3
"""Check internal documentation links and anchors.

Renaming a heading silently breaks every link to it. MkDocs `--strict` catches most of
that for the built site, but it does not run in a plain checkout, and it does not check
links inside ``README.md``, ``CONTRIBUTING.md`` or the project briefs, which are not part
of the site navigation. This script does.

Rules
-----
* Relative Markdown links must resolve to a file that exists in the repository.
* A ``#anchor`` must match a heading in the target file (GitHub's slugger rules).
* Absolute ``http(s)`` links are counted but not fetched: a network check belongs in a
  scheduled job, not in every push (and it would make the result depend on someone else's
  uptime).

    python scripts/check_docs_links.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from urllib.parse import unquote

REPO = Path(__file__).resolve().parents[1]
SCAN_SUFFIXES = {".md"}
SKIP_DIRS = {".git", ".venv", "venv", "__pycache__", "site", "node_modules", ".pytest_cache"}

LINK_RE = re.compile(r"(?<!\!)\[([^\]]*)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
HEADING_RE = re.compile(r"^(#{1,6})\s+(.*?)\s*$", re.MULTILINE)


def slug(title: str) -> str:
    """Approximate GitHub's heading-anchor slugger."""
    text = title.strip().lower()
    text = re.sub(r"`", "", text)
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"[^\w\s-]", "", text, flags=re.UNICODE)
    text = re.sub(r"\s+", "-", text.strip())
    return text


def anchors_of(path: Path) -> set[str]:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return set()
    return {slug(m.group(2)) for m in HEADING_RE.finditer(text)}


def main() -> int:
    files = [p for p in REPO.rglob("*.md") if not (set(p.parts) & SKIP_DIRS)]
    broken: list[str] = []
    internal = external = 0

    for path in sorted(files):
        text = path.read_text(encoding="utf-8", errors="ignore")
        for match in LINK_RE.finditer(text):
            target = match.group(2).strip()
            if target.startswith(("http://", "https://", "mailto:")):
                external += 1
                continue
            if target.startswith("<!--"):
                continue
            internal += 1
            file_part, _, anchor = target.partition("#")
            resolved = (path.parent / unquote(file_part)).resolve() if file_part else path.resolve()
            if file_part and not resolved.exists():
                broken.append(f"{path.relative_to(REPO)}: missing target {target!r}")
                continue
            if anchor and resolved.suffix == ".md" and slug(anchor) not in anchors_of(resolved):
                broken.append(f"{path.relative_to(REPO)}: missing anchor {target!r}")

    print(
        f"scanned {len(files)} markdown files: {internal} internal links, {external} external (not fetched)"
    )
    if broken:
        print(f"\n{len(broken)} broken link(s):")
        for item in broken:
            print(f"  - {item}")
        return 1
    print("OK: all internal links and anchors resolve.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
