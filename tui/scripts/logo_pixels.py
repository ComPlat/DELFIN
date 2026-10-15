#!/usr/bin/env python3
"""Pixel rows for the window's dolphin, derived from delfin/logo/DELFIN_logo.png.

Input: the logo (a white dolphin over two white waves on a blue disc).
Output, on stdout: the PIXELS rows of tui/src/theme/logo.ts ('#' set,
'.' clear); with --preview also the cells the window draws from them.

Method: the disc is the bounding box of the blue pixels (blue channel 40
above red). A pixel inside 0.97 of the disc's radius whose three channels
all exceed --white is white; the page around the disc is not part of the
figure. The white figure's bounding box is cropped, padded to the target
aspect, reduced by box averaging to --cols x --rows pixels, and a pixel is
set where at least --cover of its area was white. Box averaging rather
than a point sample: the fin, the flukes and the waves are a few source
pixels wide, and a point sample keeps or drops them by chance. Cell
columns and rows left empty at the edges are trimmed.

Cells: braille (U+2800-28FF, 2x4 dots per cell, dots square on screen) by
default; --cell quadrant gives 2x2 blocks (U+2596-259F, pixels twice as
tall as wide). Compared at about 20x8 cells: the quadrant grid (40x16
pixels) merges the body's upper and lower halves and the fin into one
blob; the braille grid (40x32) keeps the fin, the beak, the flukes and
both waves apart.

The mark (--cell octant, sextant, halfblock): the dolphin and its two
waves without the disc, for a terminal that draws colours. White is a
pixel with alpha over 128 and all three channels over 170; of its
8-connected parts the four largest (body, dorsal fin, the two waves) are
kept and cropped to their bounding box. That 0/255 mask is reduced by box
averaging to cols*sx x sh sub-pixels, sx x sy per cell (octant 2 x 4,
sextant 2 x 3, halfblock 1 x 2), with sh the height that keeps the
aspect on a cell twice as tall as wide, rounded up to whole cells; a
sub-pixel is set where at least 0.5 of its area is white (0.55 for
sextants: compared with 0.5, the strokes are thinner and the gaps
between the body, the fin and the tail stay open).
--cols 40 gives 40 x 18 cells (octant, sextant), 40 x 17 (halfblock).
Compared at 40 columns with the disc in half blocks (30 x 15, averaged
colours) and the disc in two-colour sextants: without the disc the
strokes are the foreground on the terminal's background, and the fin,
the tail and both waves stay apart. Octants are Unicode 16, sextants 13:
the window draws octants unless DELFIN_TUI_LOGO asks for an older kind.

Glyphs: an octant pattern (bit n-1 for octant n, octants 1-2 the top row
and 7-8 the bottom) is BLOCK OCTANT-<digits> (U+1CD00-1CDE5, 230 patterns
in ascending order) or, for the other 26, the space, the full, half and
quadrant blocks, U+2582, U+2586, U+1FB82, U+1FB85, U+1FBE6, U+1FBE7 and
U+1CEA0/1CEA3/1CEA8/1CEAB (the Unicode names of each are beside
OCTANT_ELSEWHERE). Sextants likewise (U+1FB00-1FB3B, 60 in order; the
halves and the full block not).

Colour: the per-channel median of the disc's blue pixels (opaque, not
white, blue 40 above red), and for 256 colours its nearest xterm colour
by OKLab distance among the greys and the cube colours within 8 degrees
of its OKLab hue (RGB distance over the whole cube gives lavender for
pale blues).

Output (octant, sextant, halfblock): a line 'COLOUR <hex> <xterm index>',
then the sub-pixel rows ('#' set, '.' clear).

Usage (from the repository root):
    python tui/scripts/logo_pixels.py                 # the braille rows in logo.ts
    python tui/scripts/logo_pixels.py --preview
    python tui/scripts/logo_pixels.py --cell octant     # OCTANT_PIXELS and the colour in logo.ts
    python tui/scripts/logo_pixels.py --cell octant --preview  # also the cells
    python tui/scripts/logo_pixels.py --cell sextant    # SEXTANT_PIXELS
    python tui/scripts/logo_pixels.py --cell halfblock  # HALF_PIXELS
    python tui/scripts/logo_pixels.py --cell quadrant --cols 40 --rows 16
"""

from __future__ import annotations

import argparse
import math
import statistics
from pathlib import Path

from PIL import Image

QUADRANTS = " ▘▝▀▖▌▞▛▗▚▐▜▄▙▟█"
# (dy, dx, bit) of the eight braille dots, in Unicode's numbering.
BRAILLE_DOTS = ((0, 0, 0x01), (1, 0, 0x02), (2, 0, 0x04), (0, 1, 0x08),
                (1, 1, 0x10), (2, 1, 0x20), (3, 0, 0x40), (3, 1, 0x80))
LOGO = Path(__file__).resolve().parents[2] / "delfin" / "logo" / "DELFIN_logo.png"


def white_mask(img: Image.Image, white: int) -> Image.Image:
    """An 'L' image: 255 where a pixel inside the disc is white, else 0."""
    rgb = img.convert("RGB")
    w, h = rgb.size
    px = rgb.load()
    blue = Image.new("L", (w, h), 0)
    bp = blue.load()
    for y in range(h):
        for x in range(w):
            r, _, b = px[x, y]
            if b > r + 40:
                bp[x, y] = 255
    box = blue.getbbox()
    if box is None:
        raise SystemExit("no blue disc found")
    cx, cy = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
    rad = 0.97 * min(box[2] - box[0], box[3] - box[1]) / 2
    out = Image.new("L", (w, h), 0)
    op = out.load()
    for y in range(h):
        for x in range(w):
            if min(px[x, y]) > white and (x - cx) ** 2 + (y - cy) ** 2 < rad * rad:
                op[x, y] = 255
    return out


def pixel_rows(path: Path, cols: int, rows: int, cover: float, white: int, cell: str) -> list[str]:
    mask = white_mask(Image.open(path), white)
    box = mask.getbbox()
    if box is None:
        raise SystemExit(f"{path}: no white pixel inside the disc")
    fig = mask.crop(box)
    # On screen a quadrant pixel is twice as tall as wide, a braille dot square.
    tall = 2.0 if cell == "quadrant" else 1.0
    want = cols / (rows * tall)  # target width / height in screen units
    w, h = fig.size
    if w / h < want:  # pad the sides
        nw = round(h * want)
        canvas = Image.new("L", (nw, h), 0)
        canvas.paste(fig, ((nw - w) // 2, 0))
    else:  # pad top and bottom
        nh = round(w / want)
        canvas = Image.new("L", (w, nh), 0)
        canvas.paste(fig, (0, (nh - h) // 2))
    small = canvas.resize((cols, rows), Image.Resampling.BOX)
    level = round(255 * cover)
    px = ["".join("#" if small.getpixel((x, y)) >= level else "." for x in range(cols)) for y in range(rows)]
    return trim(px, 2, 2 if cell == "quadrant" else 4)


def trim(px: list[str], cw: int, ch: int) -> list[str]:
    """Drop whole cells (cw pixels wide, ch tall) left empty at the edges."""
    def empty_cols(x: int) -> bool:
        return all(set(r[x:x + cw]) <= {"."} for r in px)

    def empty_rows(y: int) -> bool:
        return all(set(r) <= {"."} for r in px[y:y + ch])

    while px and empty_rows(0):
        px = px[ch:]
    while px and empty_rows(len(px) - ch):
        px = px[:-ch]
    while px and px[0] and empty_cols(0):
        px = [r[cw:] for r in px]
    while px and px[0] and empty_cols(len(px[0]) - cw):
        px = [r[:-cw] for r in px]
    return px


def _at(px: list[str], y: int, x: int) -> int:
    return 1 if y < len(px) and x < len(px[y]) and px[y][x] == "#" else 0


def as_quadrants(px: list[str]) -> list[str]:
    width = max(len(r) for r in px)
    return ["".join(QUADRANTS[_at(px, y, x) | _at(px, y, x + 1) << 1 | _at(px, y + 1, x) << 2
                              | _at(px, y + 1, x + 1) << 3] for x in range(0, width, 2))
            for y in range(0, len(px), 2)]


def as_braille(px: list[str]) -> list[str]:
    """An empty cell is a space, not U+2800: some fonts draw the blank
    pattern as faint dots, and trimEnd does not remove it."""
    width = max(len(r) for r in px)
    out = []
    for y in range(0, len(px), 4):
        line = ""
        for x in range(0, width, 2):
            bits = sum(m for dy, dx, m in BRAILLE_DOTS if _at(px, y + dy, x + dx))
            line += chr(0x2800 + bits) if bits else " "
        out.append(line)
    return out


XTERM_LEVELS = (0, 95, 135, 175, 215, 255)


def xterm_rgb(i: int) -> tuple[int, int, int]:
    """The colour of xterm index 16-255 (the 6x6x6 cube, then 24 greys)."""
    if i >= 232:
        v = 8 + 10 * (i - 232)
        return (v, v, v)
    j = i - 16
    return (XTERM_LEVELS[j // 36], XTERM_LEVELS[j // 6 % 6], XTERM_LEVELS[j % 6])


def oklab(c: tuple[float, float, float]) -> tuple[float, float, float]:
    def lin(u: float) -> float:
        u /= 255
        return u / 12.92 if u <= 0.04045 else ((u + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(v) for v in c)
    l_ = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m_ = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s_ = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
            1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
            0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_)


def _hue(lab: tuple[float, float, float]) -> float:
    return math.degrees(math.atan2(lab[2], lab[1]))


def _nearest(c: tuple[float, float, float], among: list[int]) -> int:
    o = oklab(c)
    return min(among, key=lambda i: sum((p - q) ** 2 for p, q in zip(o, oklab(xterm_rgb(i)))))


def logo_blue(path: Path) -> tuple[int, int, int]:
    """The disc's blue: the per-channel median of its opaque pixels that are
    not white and are blue (blue channel 40 above red)."""
    img = Image.open(path).convert("RGBA")
    blue = [p[:3] for p in img.getdata() if p[3] > 128 and min(p[:3]) <= MARK_WHITE and p[2] > p[0] + 40]
    return tuple(int(statistics.median(p[i] for p in blue)) for i in range(3))  # type: ignore[return-value]


def blue256(blue: tuple[int, int, int]) -> int:
    """The nearest xterm colour (OKLab) among the greys and the cube colours
    within 8 degrees of OKLab hue of `blue` (no lavender, no green-grey)."""
    hue = _hue(oklab(blue))
    allowed = [i for i in range(16, 256)
               if len(set(xterm_rgb(i))) == 1 or abs(_hue(oklab(xterm_rgb(i))) - hue) < 8]
    return _nearest(blue, allowed)


# The mark: the dolphin and its two waves, without the disc.
MARK_WHITE = 170
MARK_PARTS = 4
# Sub-pixels per cell (columns, rows) and the share of a sub-pixel's area
# that sets it, per kind of cell.
MARK_CELLS = {"octant": (2, 4, 0.5), "sextant": (2, 3, 0.55), "halfblock": (1, 2, 0.5)}


def mark_mask(path: Path) -> Image.Image:
    """An 'L' image of the four largest white connected parts (8-connected)
    of the logo, cropped to their bounding box: 255 set, 0 clear."""
    img = Image.open(path).convert("RGBA")
    w, h = img.size
    px = img.load()
    white = bytearray(w * h)
    for y in range(h):
        for x in range(w):
            r, g, b, a = px[x, y]
            if a > 128 and min(r, g, b) > MARK_WHITE:
                white[y * w + x] = 1
    label = [0] * (w * h)
    sizes: dict[int, int] = {}
    n = 0
    for start in range(w * h):
        if not white[start] or label[start]:
            continue
        n += 1
        label[start] = n
        stack = [start]
        size = 0
        while stack:
            i = stack.pop()
            size += 1
            y, x = divmod(i, w)
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    yy, xx = y + dy, x + dx
                    if 0 <= yy < h and 0 <= xx < w:
                        j = yy * w + xx
                        if white[j] and not label[j]:
                            label[j] = n
                            stack.append(j)
        sizes[n] = size
    keep = set(sorted(sizes, key=sizes.get, reverse=True)[:MARK_PARTS])
    out = Image.new("L", (w, h), 0)
    out.putdata([255 if label[i] in keep else 0 for i in range(w * h)])
    box = out.getbbox()
    if box is None:
        raise SystemExit(f"{path}: no white part")
    return out.crop(box)


def mark_rows(mask: Image.Image, cols: int, cell: str) -> list[str]:
    """The mark as sub-pixel rows ('#' set, '.' clear): cols cells wide,
    as many rows as keep its aspect on a cell twice as tall as wide,
    rounded up to whole cells."""
    sx, sy, cover = MARK_CELLS[cell]
    mw, mh = mask.size
    sw = cols * sx
    # A sub-pixel is (1/sx) cell wide and (2/sy) cell widths tall.
    sh = round(mh / mw * sw * (1 / sx) / (2 / sy))
    sh += -sh % sy
    small = mask.resize((sw, sh), Image.Resampling.BOX)
    level = round(255 * cover)
    return ["".join("#" if small.getpixel((x, y)) >= level else "." for x in range(sw)) for y in range(sh)]


# The octant patterns not in U+1CD00-1CDE5, by their set octants (bit n-1
# for octant n; octants numbered 1-2 in the top row, 7-8 in the bottom):
# the space, the full, half and quadrant blocks, and the eighth-row
# blocks the Unicode octant proposal names for them.
OCTANT_ELSEWHERE = {
    0x00: " ", 0xFF: "\u2588", 0x0F: "\u2580", 0xF0: "\u2584", 0x55: "\u258c", 0xAA: "\u2590",
    0x05: "\u2598", 0x0A: "\u259d", 0x50: "\u2596", 0xA0: "\u2597", 0xA5: "\u259a", 0x5A: "\u259e",
    0x5F: "\u259b", 0xAF: "\u259c", 0xF5: "\u2599", 0xFA: "\u259f",
    0x01: "\U0001cea8",  # LEFT HALF UPPER ONE QUARTER BLOCK
    0x02: "\U0001ceab",  # RIGHT HALF UPPER ONE QUARTER BLOCK
    0x40: "\U0001cea3",  # LEFT HALF LOWER ONE QUARTER BLOCK
    0x80: "\U0001cea0",  # RIGHT HALF LOWER ONE QUARTER BLOCK
    0x03: "\U0001fb82",  # UPPER ONE QUARTER BLOCK
    0x3F: "\U0001fb85",  # UPPER THREE QUARTERS BLOCK
    0xC0: "\u2582",      # LOWER ONE QUARTER BLOCK
    0xFC: "\u2586",      # LOWER THREE QUARTERS BLOCK
    0x14: "\U0001fbe6",  # MIDDLE LEFT ONE QUARTER BLOCK
    0x28: "\U0001fbe7",  # MIDDLE RIGHT ONE QUARTER BLOCK
}
SEXTANT_ELSEWHERE = {0: " ", 21: "\u258c", 42: "\u2590", 63: "\u2588"}
HALF_BLOCKS = (" ", "\u2580", "\u2584", "\u2588")


def octant_glyph(bits: int) -> str:
    """BLOCK OCTANT-<digits> (U+1CD00-1CDE5) lists the 230 others in
    ascending order of `bits`."""
    if bits in OCTANT_ELSEWHERE:
        return OCTANT_ELSEWHERE[bits]
    return chr(0x1CD00 + sum(1 for b in range(bits) if b not in OCTANT_ELSEWHERE))


def sextant_glyph(bits: int) -> str:
    """BLOCK SEXTANT-<digits> (U+1FB00-1FB3B), the same scheme with three
    rows: 60 in order, the left and right halves and the full block not."""
    if bits in SEXTANT_ELSEWHERE:
        return SEXTANT_ELSEWHERE[bits]
    return chr(0x1FB00 + bits - 1 - (bits > 21) - (bits > 42))


def as_cells(px: list[str], cell: str) -> list[str]:
    """Sub-pixel rows as lines of glyphs."""
    sx, sy, _ = MARK_CELLS[cell]
    glyph = {"octant": octant_glyph, "sextant": sextant_glyph, "halfblock": HALF_BLOCKS.__getitem__}[cell]
    out = []
    for y in range(0, len(px), sy):
        line = ""
        for x in range(0, len(px[0]), sx):
            bits = sum(1 << (dy * sx + dx) for dy in range(sy) for dx in range(sx) if px[y + dy][x + dx] == "#")
            line += glyph(bits)
        out.append(line)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cell", choices=("braille", "quadrant", "octant", "sextant", "halfblock"), default="braille")
    ap.add_argument("--cols", type=int, default=None,
                    help="braille, quadrant: pixel columns (default 40); octant, sextant, halfblock: cells (default 40)")
    ap.add_argument("--rows", type=int, default=32, help="braille, quadrant: pixel rows")
    ap.add_argument("--cover", type=float, default=0.45, help="braille, quadrant: white fraction of a pixel's area that sets it")
    ap.add_argument("--white", type=int, default=200, help="braille, quadrant: level all three channels exceed in a white pixel")
    ap.add_argument("--preview", action="store_true", help="also print the cells the window draws")
    ap.add_argument("--logo", type=Path, default=LOGO)
    a = ap.parse_args()
    cols = a.cols if a.cols is not None else 40
    if a.cell in MARK_CELLS:
        px = mark_rows(mark_mask(a.logo), cols, a.cell)
        blue = logo_blue(a.logo)
        print("COLOUR %02x%02x%02x %d" % (*blue, blue256(blue)))
        for row in px:
            print(f'  "{row}",')
        if a.preview:
            print()
            for line in as_cells(px, a.cell):
                print(line.rstrip() + "|")
        return
    px = pixel_rows(a.logo, cols, a.rows, a.cover, a.white, a.cell)
    for row in px:
        print(f'  "{row}",')
    if a.preview:
        print()
        for line in (as_braille(px) if a.cell == "braille" else as_quadrants(px)):
            print(line.rstrip() + "|")


if __name__ == "__main__":
    main()
