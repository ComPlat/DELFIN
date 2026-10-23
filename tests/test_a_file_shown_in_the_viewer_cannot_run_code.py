"""What a structure file says is shown, never run.

The 3D viewers put file content -- an .xyz comment line, a molblock --
inside a ``<script>`` the dashboard renders in the user's authenticated
session. A comment line holding a backtick, ``${...}`` or ``</script>``
used to become code there. It is a string now, whatever it holds.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from delfin.dashboard.helpers import js_string_literal

HOSTILE = (
    "3\n"
    "`+fetch('/api/kernels')+`${document.cookie}</script>"
    "<img src=x onerror=alert(1)>&amp; \n"
    "C 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0\n"
)


def test_the_literal_holds_no_markup_and_reads_back_unchanged():
    out = js_string_literal(HOSTILE)
    assert out[0] == '"' and out[-1] == '"'
    for forbidden in ("<", ">", "&", " ", "\n"):
        assert forbidden not in out
    # JSON string syntax is a subset of JavaScript string syntax, so what
    # json.loads reads back is what the browser's parser reads back.
    assert json.loads(out) == HOSTILE


def test_none_and_non_strings_become_a_string():
    assert json.loads(js_string_literal(None)) == ""
    assert json.loads(js_string_literal(42)) == "42"


def _dashboard_sources():
    root = Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
    return [root / "tab_calculations_browser.py", root / "tab_remote_archive.py"]


def test_no_viewer_puts_file_content_in_a_template_literal():
    # A template literal around interpolated file content is exactly the
    # hole this closes; none may come back.
    pattern = re.compile(r"`\{(full_xyz|data|mol_block|reference_xyz|target_xyz)\}`")
    for src in _dashboard_sources():
        assert not pattern.search(src.read_text(encoding="utf-8")), src.name


def test_every_structure_payload_goes_through_the_literal_helper():
    names = ("mol_block", "data", "reference_xyz", "target_xyz", "full_xyz")
    plain = re.compile(r"json\.dumps\((?:str\()?(" + "|".join(names) + r")\b")
    for src in _dashboard_sources():
        text = src.read_text(encoding="utf-8")
        assert not plain.search(text), (src.name, plain.search(text).group(0))
        assert "js_string_literal(" in text, src.name
