"""A row the code asks to wrap does wrap.

Input: every ``Layout(...)`` in delfin/. ipywidgets' Layout (8.1) has no
``flex_wrap`` trait: ``Layout(flex_wrap='wrap')`` only warns
(DeprecationWarning, unrecognized argument) and the row stays
``nowrap``. Measured on the Submit Job tab at 1366 px: the seven
structure-builder buttons (MANTA ... METALLOGEN, 150 px each) shrank to
86 px and cut their names ("ARCHIT...", "MOLAS..."), and the GOAT and
Fukui SUBMIT buttons read "SUBMIT...", in rows written with
``flex_wrap='wrap'``.

Output pinned: no Layout in delfin/ passes ``flex_wrap``; a wrapping row
says ``flex_flow='row wrap'``, the trait the widgets send to the page.
"""
from __future__ import annotations

import re
from pathlib import Path

import ipywidgets as widgets

ROOT = Path(__file__).resolve().parents[1] / "delfin"


def test_the_widgets_have_no_flex_wrap_trait():
    keys = widgets.Layout().keys
    assert "flex_flow" in keys
    assert "flex_wrap" not in keys


def test_no_layout_passes_flex_wrap():
    found = []
    for path in sorted(ROOT.rglob("*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for m in re.finditer(r"\bflex_wrap\s*=", text):
            found.append(f"{path.relative_to(ROOT)}:{text.count(chr(10), 0, m.start()) + 1}")
    assert not found, found
