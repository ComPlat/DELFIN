"""The window's dolphin (tui/src/theme/logo.ts) is the logo's, reduced.

The PIXELS, OCTANT_PIXELS, SEXTANT_PIXELS and HALF_PIXELS rows and the
colour in logo.ts are pasted output of tui/scripts/logo_pixels.py run on
delfin/logo/DELFIN_logo.png (defaults; --cell octant, sextant, halfblock). These tests run the
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


MARK = {"octant": ("OCTANT_PIXELS", 2, 4), "sextant": ("SEXTANT_PIXELS", 2, 3), "halfblock": ("HALF_PIXELS", 1, 2)}


def _pasted_mark(cell: str) -> list[str]:
    block = re.search(rf"const {MARK[cell][0]} = \[(.*?)\]", LOGO_TS.read_text(encoding="utf-8"), re.S)
    assert block, f"logo.ts has no {MARK[cell][0]}"
    return re.findall(r'"([#.]+)"', block.group(1))


def _script_mark(cell: str) -> tuple[str, list[str]]:
    out = subprocess.run([sys.executable, str(SCRIPT), "--cell", cell], capture_output=True, text=True, check=True).stdout
    return out.split("\n", 1)[0], re.findall(r'"([#.]+)"', out)


def _script_module():
    sys.path.insert(0, str(SCRIPT.parent))
    try:
        import logo_pixels
    finally:
        sys.path.remove(str(SCRIPT.parent))
    return logo_pixels


def test_the_pasted_mark_rows_are_the_scripts_output_for_every_kind_of_cell():
    for cell in MARK:
        assert _pasted_mark(cell) == _script_mark(cell)[1], cell


def test_the_pasted_colour_is_the_scripts():
    text = LOGO_TS.read_text(encoding="utf-8")
    rgb = re.search(r'export const LOGO_RGB = "#([0-9a-f]{6})"', text)
    x256 = re.search(r"export const LOGO_256 = (\d+)", text)
    assert rgb and x256
    assert _script_mark("octant")[0] == f"COLOUR {rgb.group(1)} {x256.group(1)}"
    r, g, b = (int(rgb.group(1)[i:i + 2], 16) for i in (0, 2, 4))
    assert b > r + 40 and b > g        # a blue
    assert 16 <= int(x256.group(1)) <= 231


def test_the_mark_is_forty_cells_wide_and_keeps_its_aspect():
    rows = {}
    for cell, (_, sx, sy) in MARK.items():
        px = _pasted_mark(cell)
        assert len({len(r) for r in px}) == 1 and len(px[0]) == 40 * sx, cell
        assert len(px) % sy == 0, cell
        rows[cell] = len(px) // sy
        # Something in the first two and the last two columns of cells and
        # in the first row: cropped to the mark.
        assert "#" in "".join(px[:sy]), cell
        assert any("#" in r[:2 * sx] for r in px) and any("#" in r[-2 * sx:] for r in px), cell
    assert rows["octant"] == rows["sextant"] == 18
    assert rows["halfblock"] in (17, 18)


def test_every_octant_pattern_is_one_glyph_that_decodes_back():
    lp = _script_module()
    glyphs = [lp.octant_glyph(b) for b in range(256)]
    assert len(set(glyphs)) == 256
    # The 230 in U+1CD00-1CDE5 are named BLOCK OCTANT-<digits> by Unicode
    # 16; the rank rule reproduces these entries of UnicodeData.txt.
    named = {0x1CD00: "3", 0x1CD01: "23", 0x1CD02: "123", 0x1CD03: "4", 0x1CDE4: "1345678", 0x1CDE5: "2345678"}
    for cp, digits in named.items():
        assert glyphs[sum(1 << (int(d) - 1) for d in digits)] == chr(cp)
    in_block = sorted(ord(g) for g in glyphs if 0x1CD00 <= ord(g) <= 0x1CDE5)
    assert in_block == list(range(0x1CD00, 0x1CDE6))
    # The other 26, decoded from their names' geometry: (columns, rows of
    # four) each covers.
    elsewhere = {
        " ": "", "\u2588": "12345678", "\u2580": "1234", "\u2584": "5678", "\u258c": "1357", "\u2590": "2468",
        "\u2598": "13", "\u259d": "24", "\u2596": "57", "\u2597": "68", "\u259a": "1368", "\u259e": "2457",
        "\u259b": "123457", "\u259c": "123468", "\u2599": "135678", "\u259f": "245678",
        "\U0001cea8": "1", "\U0001ceab": "2", "\U0001cea3": "7", "\U0001cea0": "8",
        "\U0001fb82": "12", "\U0001fb85": "123456", "\u2582": "78", "\u2586": "345678",
        "\U0001fbe6": "35", "\U0001fbe7": "46",
    }
    decode = {g: b for b, g in enumerate(glyphs)}
    for g, digits in elsewhere.items():
        assert decode[g] == sum(1 << (int(d) - 1) for d in digits), hex(ord(g[0])) if g else "space"
    rest = [g for g in glyphs if g not in elsewhere]
    assert len(rest) == 230 and all(0x1CD00 <= ord(g) <= 0x1CDE5 for g in rest)


def test_every_sextant_pattern_is_one_glyph():
    lp = _script_module()
    glyphs = [lp.sextant_glyph(b) for b in range(64)]
    assert len(set(glyphs)) == 64
    assert [glyphs[0], glyphs[21], glyphs[42], glyphs[63]] == [" ", "\u258c", "\u2590", "\u2588"]
    others = [g for b, g in enumerate(glyphs) if b not in (0, 21, 42, 63)]
    assert [ord(g) for g in others] == list(range(0x1FB00, 0x1FB3C))
