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


def _array(text: str, name: str) -> list[str]:
    block = re.search(rf"const {name} = \[(.*?)\]", text, re.S)
    assert block, f"logo.ts has no {name}"
    return re.findall(r'"([.0-9a-f]+)"', block.group(1))


def _pasted_color() -> tuple[list[str], list[str]]:
    text = LOGO_TS.read_text(encoding="utf-8")
    return _array(text, "COLOR_PIXELS"), _array(text, "COLOR_256")


def _script_color() -> tuple[list[str], list[str]]:
    out = subprocess.run([sys.executable, str(SCRIPT), "--cell", "halfblock"],
                         capture_output=True, text=True, check=True).stdout
    rgb, idx = out.split("COLOR_PIXELS\n", 1)[1].split("COLOR_256\n", 1)
    return re.findall(r'"([.0-9a-f]+)"', rgb), re.findall(r'"([.0-9]+)"', idx)


def test_the_pasted_colour_rows_are_the_scripts_output():
    assert _pasted_color() == _script_color()


def _cells(row: str, width: int) -> list[str | None]:
    return [None if set(row[i:i + width]) == {"."} else row[i:i + width] for i in range(0, len(row), width)]


def test_the_colour_logo_is_thirty_cells_by_fifteen_of_averaged_pixels():
    rgb, idx = _pasted_color()
    assert len(rgb) == len(idx) == 30
    px = [_cells(r, 6) for r in rgb]
    ix = [_cells(r, 3) for r in idx]
    assert {len(r) for r in px} == {len(r) for r in ix} == {30}
    # The same pixels are outside the disc in both, and the corners are.
    assert [[c is None for c in r] for r in px] == [[c is None for c in r] for r in ix]
    for y, x in ((0, 0), (0, 29), (29, 0), (29, 29)):
        assert px[y][x] is None
    assert px[15][15] is not None
    colours = {c for r in px for c in r if c is not None}
    # Averaged, not thresholded: many tones between the disc's blue and white.
    assert len(colours) > 100
    rgbs = [tuple(int(c[i:i + 2], 16) for i in (0, 2, 4)) for c in colours]
    assert any(min(v) > 235 for v in rgbs)           # the dolphin's white
    assert sum(1 for r, _, b in rgbs if b > r + 40) > len(rgbs) // 2   # the disc's blue
    assert all(16 <= int(c) <= 255 for r in ix for c in r if c is not None)
