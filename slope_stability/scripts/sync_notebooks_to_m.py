#!/usr/bin/env python3
"""Sync MATLAB/Octave scripts from same-named Jupyter notebooks."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent


def _cell_source(cell: dict) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(source)
    return str(source)


def _render_markdown(markdown: str) -> str:
    lines = markdown.splitlines()
    rendered = []
    for line in lines:
        stripped = line.rstrip()
        if not stripped:
            rendered.append("%")
            continue

        heading = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if heading:
            title = heading.group(2).strip()
            rendered.append(f"%% {title}" if title else "%%")
            continue

        rendered.append(f"% {stripped}")

    return "\n".join(rendered).rstrip()


def render_notebook(notebook_path: Path) -> str:
    data = json.loads(notebook_path.read_text(encoding="utf-8"))
    parts = [
        f"% Auto-generated from {notebook_path.name}.",
        "% Update the notebook first, then re-run scripts/sync_notebooks_to_m.py.",
    ]

    for cell in data.get("cells", []):
        cell_type = cell.get("cell_type")
        source = _cell_source(cell)

        if cell_type == "markdown":
            rendered = _render_markdown(source)
        elif cell_type == "code":
            rendered = source.rstrip()
        else:
            continue

        if rendered:
            parts.append(rendered)

    return "\n\n".join(parts).rstrip() + "\n"


def sync_notebook(notebook_path: Path) -> bool:
    matlab_path = notebook_path.with_suffix(".m")
    rendered = render_notebook(notebook_path)
    current = matlab_path.read_text(encoding="utf-8") if matlab_path.exists() else None
    if current == rendered:
        return False
    matlab_path.write_text(rendered, encoding="utf-8")
    return True


def find_notebooks(paths: list[str]) -> list[Path]:
    if paths:
        notebooks = []
        for raw_path in paths:
            path = Path(raw_path)
            if not path.is_absolute():
                path = (Path.cwd() / path).resolve()
            if path.suffix != ".ipynb":
                raise ValueError(f"Expected a notebook path, got: {raw_path}")
            notebooks.append(path)
        return sorted(notebooks)

    return sorted(ROOT_DIR.glob("*.ipynb"))


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "notebooks",
        nargs="*",
        help="Optional notebook paths to sync. Defaults to all notebooks in slope_stability/.",
    )
    args = parser.parse_args(argv)

    changed = 0
    notebooks = find_notebooks(args.notebooks)
    for notebook_path in notebooks:
        if sync_notebook(notebook_path):
            changed += 1
            print(f"updated {notebook_path.with_suffix('.m').name}")
        else:
            print(f"unchanged {notebook_path.with_suffix('.m').name}")

    print(f"sync complete: {changed} file(s) updated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
