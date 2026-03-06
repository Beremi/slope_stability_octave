#!/usr/bin/env python3
"""Execute a minimal Octave notebook cell and verify PNG plot output."""

from __future__ import annotations

import argparse
import base64
from typing import Iterable

import nbformat
from nbclient import NotebookClient


BASE_TEST_CELL = """%plot -f png -r 120
pkg load sparsersb;
try
  graphics_toolkit('qt');
catch err
  disp(['qt_error=', err.message]);
  graphics_toolkit('gnuplot');
end
set(0, 'defaultfigurevisible', 'off');
disp(['toolkit=', graphics_toolkit()]);
figure;
plot(1:10);
"""

LATEX_LABEL_CELL = """xlabel('Control variable - $\\omega$', 'interpreter', 'latex');
ylabel('strength reduction factor - $\\lambda$', 'interpreter', 'latex');
title('devcontainer plot test', 'interpreter', 'latex');
"""

PLAIN_LABEL_CELL = """title('devcontainer plot test');
"""


def flatten_text(outputs: Iterable[dict]) -> str:
    chunks: list[str] = []
    for output in outputs:
        output_type = output.get("output_type")
        if output_type == "stream":
            chunks.append(output.get("text", ""))
            continue
        if output_type not in {"display_data", "execute_result"}:
            continue
        data = output.get("data", {})
        text = data.get("text/plain")
        if text is not None:
            chunks.append(str(text))
    return "".join(chunks).strip()


def find_png(outputs: Iterable[dict]) -> bytes | None:
    for output in outputs:
        if output.get("output_type") not in {"display_data", "execute_result"}:
            continue
        data = output.get("data", {})
        encoded = data.get("image/png")
        if encoded:
            return base64.b64decode(encoded)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--kernel",
        default="octave-local-rsb",
        help="Jupyter kernel name to execute (default: %(default)s)",
    )
    parser.add_argument(
        "--workdir",
        default=".",
        help="Notebook execution directory (default: current directory)",
    )
    parser.add_argument(
        "--latex-labels",
        action="store_true",
        help="Exercise Octave's LaTeX text renderer in the notebook plot.",
    )
    args = parser.parse_args()

    test_cell = BASE_TEST_CELL
    if args.latex_labels:
        test_cell += LATEX_LABEL_CELL
    else:
        test_cell += PLAIN_LABEL_CELL
    test_cell += "drawnow;\n"

    nb = nbformat.v4.new_notebook()
    nb.cells = [nbformat.v4.new_code_cell(test_cell)]

    client = NotebookClient(
        nb,
        kernel_name=args.kernel,
        timeout=300,
        resources={"metadata": {"path": args.workdir}},
    )
    executed = client.execute()
    outputs = executed.cells[0].get("outputs", [])
    text = flatten_text(outputs)
    png = find_png(outputs)

    print("OUTPUT_TEXT_BEGIN")
    print(text)
    print("OUTPUT_TEXT_END")
    print(f"PNG_COUNT={1 if png else 0}")
    if not png:
        raise SystemExit("No PNG plot output was produced by the Jupyter kernel.")

    print(f"PNG_BYTES={len(png)}")
    print(f"PNG_SIGNATURE={png[:8].hex()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
