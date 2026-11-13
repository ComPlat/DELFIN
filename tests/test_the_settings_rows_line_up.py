"""The Settings tab's rows line up: a label column, the field beside it.

Input: the rows of the Settings sections (tab_settings: HBoxes of a bold
label, a field, sometimes a button and a hint). Measured in the real page
at 1366 px with Runtime Backend open: every "<Tool> binary" field took a
line of its own below its label with its Scan button on a third line, and
the "Detected <Tool>" dropdowns started wherever their label ended (x = 161
to 186), so no two rows lined up.

Output pinned: the tab carries the class ``delfin-settings``; its CSS gives
the first label of a section row a minimum width (a label column), and a
full-width field in a row the room left beside it (flex basis 280 px,
growing), so label, field and button share one line where the row is
wide enough and wrap only where it is not.
"""
from __future__ import annotations

from pathlib import Path

SOURCE = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
          / "tab_settings.py").read_text(encoding="utf-8")


def test_the_tab_is_marked():
    assert "tab.add_class('delfin-settings')" in SOURCE


def test_a_label_column_and_fields_beside_it():
    css = " ".join(SOURCE.split())
    assert (".delfin-settings .widget-accordion .widget-hbox > .widget-html:first-child "
            "{ min-width: 190px; }") in css
    assert ('.delfin-settings .widget-accordion .widget-hbox > .widget-text[style*="width: 100%"] '
            "{ flex: 1 1 280px; }") in css
