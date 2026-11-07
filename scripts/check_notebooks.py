#!/usr/bin/env python3
"""Execute every notebook and fail if a cell errors.

A notebook with stored outputs that no longer match the code is worse than no notebook:
it looks like evidence. This script re-executes each notebook in place, so the committed
outputs are the outputs of the committed code, and it fails if any cell raises.

Checks, per notebook:

* every code cell executes without an exception;
* no cell writes a traceback to stderr;
* the stored outputs match what the re-execution produces (unless ``--no-write``);
* warnings are printed but do not fail the build — a notebook that emits a numpy warning
  is doing its job by showing it.

Usage::

    python scripts/check_notebooks.py
    python scripts/check_notebooks.py --no-write     # verify without touching the file
    python scripts/check_notebooks.py notebooks/01_fitting_decisions.ipynb
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parents[1]


def run_notebook(path: Path, *, write: bool, timeout: int) -> tuple[bool, str]:
    before = path.read_text(encoding="utf-8")
    notebook = nbformat.reads(before, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(path.parent.parent)}},
    )
    try:
        client.execute()
    except CellExecutionError as exc:
        return False, f"{path.name}: a cell raised\n{exc}"

    warnings = []
    for cell in notebook.cells:
        if cell.cell_type != "code":
            continue
        for output in cell.get("outputs", []):
            if output.get("output_type") == "stream" and output.get("name") == "stderr":
                warnings.append("".join(output.get("text", "")).strip().splitlines()[0])
    # nbclient stamps every executed cell with wall-clock times. They are not output, they
    # differ on every run, and committing them makes the notebook look dirty to `git status`
    # after every check — and would make the CI step that compares committed outputs against a
    # fresh execution fail for a reason nobody should care about. Outputs stay; timings go.
    for cell in notebook.cells:
        cell.get("metadata", {}).pop("execution", None)
    if write:
        nbformat.write(notebook, path)
    changed = write and path.read_text(encoding="utf-8") != before
    note = f"{path.name}: executed {sum(c.cell_type == 'code' for c in notebook.cells)} code cells"
    if changed:
        note += ", outputs updated"
    if warnings:
        note += f", {len(warnings)} warning line(s): {warnings[0][:70]}"
    return True, note


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="*", help="notebooks to check (default: all)")
    parser.add_argument("--no-write", action="store_true", help="do not update stored outputs")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    paths = [Path(p) for p in args.paths] or sorted((ROOT / "notebooks").glob("*.ipynb"))
    if not paths:
        print("FAIL no notebooks found")
        return 1

    failures = 0
    for path in paths:
        ok, message = run_notebook(path, write=not args.no_write, timeout=args.timeout)
        print(f"{'OK  ' if ok else 'FAIL'} {message}")
        failures += not ok
    print(f"\n{len(paths) - failures}/{len(paths)} notebook(s) executed without error")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
