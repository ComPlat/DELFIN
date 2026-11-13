"""The file browser's name prompts keep their label and fit their pane.

Input: the New folder / Rename / Duplicate prompt rows of the file
browsers (tab_calculations_browser: Calculations, Archive, Office).
Measured in the real page at 1366 px (page_check driver): the row's
"New folder name:" label shrank to 0-6 px -- the input took the room and
the HTML label has no minimum -- so the field stood alone without its
name, and Cancel ran 6 px past the pane's edge.

Output pinned: each label has a line of its own above the field
(flex 0 0 100% in a wrapping row, as forms label a field), and the field
(at least 160 px) and its two buttons share the line below, which fits
the narrowest pane (160 + 70 + 78 px).
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


def test_each_prompt_field_keeps_a_usable_width():
    for name in ("calc_rename_prompt_input", "calc_duplicate_prompt_input",
                 "calc_new_folder_prompt_input"):
        assert "min_width='160px'" in _call(name), name


def test_each_prompt_row_wraps_instead_of_overflowing():
    for name in ("calc_rename_prompt_row", "calc_duplicate_prompt_row",
                 "calc_new_folder_prompt_row"):
        call = SOURCE[SOURCE.index(f"    {name} = widgets.HBox("):]
        call = call[:call.index("\n    )\n")]
        assert re.search(r"flex_flow='row wrap'", call), name
