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

Colour (--cell halfblock): the logo itself, for a terminal that draws
colours. The PNG's alpha channel is the disc (opaque inside, transparent
outside); its bounding box is cropped, padded to a square and reduced to
--size x --size pixels by Lanczos resampling of the premultiplied RGBA
(RGBa), so every pixel is the area average of the source under it and a
rim pixel's colour is not darkened by the transparent black beyond it.
No threshold: the dolphin's strokes, a pixel or less wide here, keep their
mixed blue-white tone, and that tone is what draws the fin, the beak and
the waves. A pixel with alpha under 50% after the reduction is outside the
disc and drawn as the terminal's own background. Each cell is an
upper-half block U+2580 with the upper pixel as its foreground and the
lower one as its background, so pixels are square on screen: --size 30 is
30 x 15 cells. Compared at 24, 30 and 36 pixels: 30 is the smallest at
which the dolphin reads as one; the earlier thresholded white/blue figure
(26 x 13 cells, five tone steps) did not read as a dolphin at all.

256 colours: each pixel's nearest xterm colour by OKLab distance, among
a ramp: the colours (cube 16-231 or grey 232-255) nearest to 33 blends
from the disc's median blue to white, where only greys and cube colours
within 8 degrees of OKLab hue of the disc's blue are candidates. Rejected:
RGB distance over the whole cube and grey ramp (pale blues become 146,
189, lavender, or 152, a green-grey); OKLab distance over the whole set
(lavender 103 and 146 remain, and the disc's blue, between two cube
levels, flips between 67, 68 and 74 from pixel to pixel, a speckled rim).

Output (halfblock): COLOR_PIXELS, one string per pixel row, six hex
digits per pixel ('......' outside the disc); COLOR_256, the same pixels
as three-digit xterm indices ('...' outside).

Usage (from the repository root):
    python tui/scripts/logo_pixels.py                 # the braille rows in logo.ts
    python tui/scripts/logo_pixels.py --preview
    python tui/scripts/logo_pixels.py --cell halfblock  # the colour rows in logo.ts
    python tui/scripts/logo_pixels.py --cell halfblock --preview  # also as ANSI colour cells
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


def color_pixels(path: Path, size: int) -> list[list[tuple[int, int, int] | None]]:
    """The disc as `size` x `size` area-averaged pixels; None outside it."""
    img = Image.open(path).convert("RGBA")
    box = img.getchannel("A").point(lambda v: 255 if v >= 128 else 0).getbbox()
    if box is None:
        raise SystemExit(f"{path}: no opaque pixel")
    disc = img.crop(box)
    w, h = disc.size
    side = max(w, h)
    square = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    square.paste(disc, ((side - w) // 2, (side - h) // 2))
    small = square.convert("RGBa").resize((size, size), Image.Resampling.LANCZOS).convert("RGBA")
    px = small.load()
    return [[px[x, y][:3] if px[x, y][3] >= 128 else None for x in range(size)] for y in range(size)]


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


def ramp256(pixels: list[list[tuple[int, int, int] | None]]) -> list[int]:
    """The xterm indices from the disc's median blue to white (see the
    module's docstring): greys and cube colours of the blue's hue."""
    disc = [p for row in pixels for p in row if p is not None]
    blue = tuple(statistics.median(p[i] for p in disc) for i in range(3))
    hue = _hue(oklab(blue))
    allowed = [i for i in range(16, 256)
               if len(set(xterm_rgb(i))) == 1 or abs(_hue(oklab(xterm_rgb(i))) - hue) < 8]
    ramp: list[int] = []
    for k in range(33):
        t = k / 32
        i = _nearest(tuple(v + (255 - v) * t for v in blue), allowed)
        if i not in ramp:
            ramp.append(i)
    return ramp


def color_rows(path: Path, size: int) -> tuple[list[str], list[str]]:
    """COLOR_PIXELS and COLOR_256 rows (see the module's docstring)."""
    pixels = color_pixels(path, size)
    ramp = ramp256(pixels)
    rgb = ["".join("......" if p is None else "%02x%02x%02x" % p for p in row) for row in pixels]
    idx = ["".join("..." if p is None else f"{_nearest(p, ramp):03d}" for p in row) for row in pixels]
    return rgb, idx


def as_halfblocks(rows: list[str], width: int) -> list[str]:
    """ANSI colour lines of COLOR_PIXELS (width 6) or COLOR_256 (width 3)
    rows: the preview."""
    def cells(row: str) -> list[str | None]:
        return [None if set(row[i:i + width]) == {"."} else row[i:i + width] for i in range(0, len(row), width)]

    def sgr(c: str, bg: bool) -> str:
        base = 48 if bg else 38
        if width == 3:
            return f"{base};5;{int(c)}"
        return f"{base};2;{int(c[0:2], 16)};{int(c[2:4], 16)};{int(c[4:6], 16)}"
    out = []
    for y in range(0, len(rows), 2):
        top = cells(rows[y])
        low = cells(rows[y + 1]) if y + 1 < len(rows) else [None] * len(top)
        line = ""
        for a, b in zip(top, low):
            if a is None and b is None:
                line += "\x1b[0m "
            elif a is None:
                line += f"\x1b[0m\x1b[{sgr(b, False)}m▄"
            elif b is None:
                line += f"\x1b[0m\x1b[{sgr(a, False)}m▀"
            else:
                line += f"\x1b[{sgr(a, False)};{sgr(b, True)}m▀"
        out.append(line + "\x1b[0m")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cell", choices=("braille", "quadrant", "halfblock"), default="braille")
    ap.add_argument("--cols", type=int, default=40, help="pixel columns (two per braille or quadrant cell)")
    ap.add_argument("--rows", type=int, default=32, help="pixel rows (four per braille cell, two per quadrant cell)")
    ap.add_argument("--size", type=int, default=30,
                    help="halfblock: pixels per side of the disc (one column, half a row each)")
    ap.add_argument("--cover", type=float, default=0.45, help="white fraction of a pixel's area that sets it")
    ap.add_argument("--white", type=int, default=200, help="level all three channels exceed in a white pixel")
    ap.add_argument("--preview", action="store_true", help="also print the cells the window draws")
    ap.add_argument("--logo", type=Path, default=LOGO)
    a = ap.parse_args()
    if a.cell == "halfblock":
        rgb, idx = color_rows(a.logo, a.size)
        print("COLOR_PIXELS")
        for row in rgb:
            print(f'  "{row}",')
        print("COLOR_256")
        for row in idx:
            print(f'  "{row}",')
        if a.preview:
            print()
            print("\n".join(as_halfblocks(rgb, 6)))
            print()
            print("\n".join(as_halfblocks(idx, 3)))
        return
    px = pixel_rows(a.logo, a.cols, a.rows, a.cover, a.white, a.cell)
    for row in px:
        print(f'  "{row}",')
    if a.preview:
        print()
        for line in (as_braille(px) if a.cell == "braille" else as_quadrants(px)):
            print(line.rstrip() + "|")


if __name__ == "__main__":
    main()
