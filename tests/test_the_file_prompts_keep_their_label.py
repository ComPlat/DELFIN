"""The file browser's name prompts keep their label and fit their pane.

Input: the New folder / Rename / Duplicate prompt rows of the file
browsers (tab_calculations_browser: Calculations, Archive, Office).
Measured in the real page at 1366 px (page_check driver): the row's
"New folder name:" label shrank to 0-6 px -- the input took the room and
the HTML label has no minimum -- so the field stood alone without its
name, and Cancel ran 6 px past the pane's edge.

Output pinned: each label has a line of its own above the field
(flex 0 0 100% in a wrapping row, as forms label a field), and the field
(basis 140 px) and its two buttons share the line below (140 + 70 + 78 px
fit the pane at its default width).
"""
from __future__ import annotations

import re
from pathlib import Path

SOURCE = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
          / "tab_calculations_browser.py").read_text(encoding="utf-8")


def _call(name: str) -> str:
    start = SOURCE.index(f"    {name} = widgets.")
    end = SOURCE.index("\n    )\n", start)
    return SOURCE[start:end]


def test_each_prompt_label_keeps_its_width():
    for name in ("calc_rename_prompt_label", "calc_duplicate_prompt_label",
                 "calc_new_folder_prompt_label"):
        call = SOURCE[SOURCE.index(f"    {name} = widgets.HTML("):]
        call = call[:call.index("\n")]
        assert "flex='0 0 100%'" in call, name


def test_each_prompt_field_shares_its_line_with_the_two_buttons():
    # The left pane's own rule (.calc-left .widget-text, !important) sets
    # width:auto, so the input's intrinsic ~210 px decided the wrap and
    # Cancel took a line of its own (measured); the prompt rows' rule is
    # more specific and gives the field a 140 px basis.
    for name in ("calc_rename_prompt_row", "calc_duplicate_prompt_row",
                 "calc_new_folder_prompt_row"):
        assert name in SOURCE[SOURCE.index("for _prompt_row in ("):][:200], name
    assert "_prompt_row.add_class('calc-prompt-row')" in SOURCE
    assert ("'.calc-left .calc-prompt-row > .widget-text { flex:1 1 140px !important; "
            "min-width:140px !important; width:auto !important; }'") in SOURCE


def test_each_prompt_row_wraps_instead_of_overflowing():
    for name in ("calc_rename_prompt_row", "calc_duplicate_prompt_row",
                 "calc_new_folder_prompt_row"):
        call = SOURCE[SOURCE.index(f"    {name} = widgets.HBox("):]
        call = call[:call.index("\n    )\n")]
        assert re.search(r"flex_flow='row wrap'", call), name
