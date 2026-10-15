"""The window's dolphin (tui/src/theme/logo.ts) is the logo's, reduced.

The PIXELS rows in logo.ts are pasted output of tui/scripts/logo_pixels.py
run on delfin/logo/DELFIN_logo.png with its defaults. These tests run the
script and compare, so a changed logo or a hand edit of the rows shows as
a difference rather than as a drawing nobody can reproduce.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOGO_TS = ROOT / "tui" / "src" / "theme" / "logo.ts"
SCRIPT = ROOT / "tui" / "scripts" / "logo_pixels.py"


def _pasted_rows() -> list[str]:
    text = LOGO_TS.read_text(encoding="utf-8")
    block = re.search(r"const PIXELS = \[(.*?)\]", text, re.S)
    assert block, "logo.ts has no PIXELS array"
    return re.findall(r'"([#.]+)"', block.group(1))


def _script_rows() -> list[str]:
    out = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True, check=True).stdout
    return re.findall(r'"([#.]+)"', out)


def test_the_pasted_rows_are_the_scripts_output():
    assert _pasted_rows() == _script_rows()


def test_the_rows_fill_eighteen_to_twenty_two_braille_cells_by_six_to_eight():
    rows = _pasted_rows()
    assert len({len(r) for r in rows}) == 1
    cols, lines = len(rows[0]) // 2, len(rows) // 4
    assert len(rows) % 4 == 0
    assert 18 <= cols <= 22 and 6 <= lines <= 8, (cols, lines)


def test_logo_ts_names_the_script_that_made_it():
    assert "tui/scripts/logo_pixels.py" in LOGO_TS.read_text(encoding="utf-8")
