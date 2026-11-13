"""The Submit Job tab's small job forms line up in columns.

Input: the GOAT, CO2 Coordinator, GUPPY and Fukui rows of tab_submit
(Charge / Solvent / Converter, PAL / MaxCore / Time Limit, ...). Measured
in the real page at 1366 px: each field sized to its own label, so the
inputs of two rows started at different places (Charge's at x=751, PAL's
at x=729; MaxCore's at 905, Solvent's at 916) and the same three fields
of GOAT, GUPPY and Fukui never lined up.

Output pinned: those rows carry the class ``submit-form-grid``; the tab's
CSS lays each one out as a grid of equal 190 px columns (three fit the
right-hand pane at 1366 px) with every field filling its cell and every
label at least 84 px wide, so labels, inputs and buttons start at the same
x in every row.
"""
from __future__ import annotations

from pathlib import Path

SOURCE = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
          / "tab_submit.py").read_text(encoding="utf-8")


def _css() -> str:
    return " ".join(SOURCE.split())


def test_every_job_form_row_is_a_grid_row():
    rows = SOURCE.count("_form_row(")
    # GOAT (two rows), CO2, GUPPY, Fukui (two rows); plus the definition.
    assert rows == 7, rows


def test_the_grid_has_equal_columns_and_filling_fields():
    css = _css()
    assert ".submit-split-pane .submit-form-grid { display: grid;" in css
    assert "grid-template-columns: repeat(auto-fill, 190px);" in css
    assert ".submit-split-pane .submit-form-grid > .jupyter-widgets { width: auto !important;" in css
    assert ".submit-split-pane .submit-form-grid .widget-label { min-width: 84px;" in css
