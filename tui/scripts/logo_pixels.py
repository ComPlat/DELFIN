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

Usage (from the repository root):
    python tui/scripts/logo_pixels.py                 # the rows in logo.ts
    python tui/scripts/logo_pixels.py --preview
    python tui/scripts/logo_pixels.py --cell quadrant --cols 40 --rows 16
"""

from __future__ import annotations

import argparse
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cell", choices=("braille", "quadrant"), default="braille")
    ap.add_argument("--cols", type=int, default=40, help="pixel columns (two per cell)")
    ap.add_argument("--rows", type=int, default=32, help="pixel rows (four per braille cell, two per quadrant cell)")
    ap.add_argument("--cover", type=float, default=0.45, help="white fraction of a pixel's area that sets it")
    ap.add_argument("--white", type=int, default=200, help="level all three channels exceed in a white pixel")
    ap.add_argument("--preview", action="store_true", help="also print the cells the window draws")
    ap.add_argument("--logo", type=Path, default=LOGO)
    a = ap.parse_args()
    px = pixel_rows(a.logo, a.cols, a.rows, a.cover, a.white, a.cell)
    for row in px:
        print(f'  "{row}",')
    if a.preview:
        print()
        for line in (as_braille(px) if a.cell == "braille" else as_quadrants(px)):
            print(line.rstrip() + "|")


if __name__ == "__main__":
    main()
