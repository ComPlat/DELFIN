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

Colour (--cell halfblock): the disc itself, for a terminal that draws
colours. Each cell is an upper-half block U+2580 with the upper pixel as
its foreground and the lower one as its background, so pixels are square
on screen: --cols 26 --rows 26 is 26 x 13 cells. A pixel belongs to the
figure where at least half its area lies inside the disc (the page
around it stays the terminal's own background). Its tone is the white
fraction of its area, stretched from --lo..--hi to 0..1 (a box average
alone turns the strokes, a pixel or less wide here, into pale blue) and
cut to --levels steps between the disc's mean blue and white. Output: the
COLOR_PIXELS rows ('.' outside the disc, '0'.. the step) and the
COLOR_PALETTE hex colours of logo.ts. Compared at 22x11, 24x12 and 26x13
cells, rendered to an image: 26x13 is the smallest at which the fin, the
beak and both waves stay apart.

Usage (from the repository root):
    python tui/scripts/logo_pixels.py                 # the braille rows in logo.ts
    python tui/scripts/logo_pixels.py --preview
    python tui/scripts/logo_pixels.py --cell halfblock  # the colour rows in logo.ts
    python tui/scripts/logo_pixels.py --cell quadrant --cols 40 --rows 16
"""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image, ImageDraw

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


def disc_box(rgb: Image.Image) -> tuple[float, float, float]:
    """Centre and radius of the blue disc (blue channel 40 above red)."""
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
    return (box[0] + box[2]) / 2, (box[1] + box[3]) / 2, min(box[2] - box[0], box[3] - box[1]) / 2


def color_rows(path: Path, cols: int, rows: int, white: int, lo: float, hi: float,
               levels: int) -> tuple[list[str], list[str]]:
    """The disc as `rows` x `cols` pixels: '.' outside it, else the tone
    step '0'..; and the palette, one hex colour per step (the disc's mean
    blue to white)."""
    rgb = Image.open(path).convert("RGB")
    cx, cy, rad = disc_box(rgb)
    box = (round(cx - rad), round(cy - rad), round(cx + rad), round(cy + rad))
    crop = rgb.crop(box)
    side = crop.size[0]
    inside = Image.new("L", crop.size, 0)
    ImageDraw.Draw(inside).ellipse((0, 0, side - 1, crop.size[1] - 1), fill=255)
    # White only well inside the rim: the page's white beyond the disc's
    # antialiased edge would read as a pale ring.
    rim = Image.new("L", crop.size, 0)
    m = 0.03 * side / 2
    ImageDraw.Draw(rim).ellipse((m, m, side - 1 - m, crop.size[1] - 1 - m), fill=255)
    cp, ip, rp = crop.load(), inside.load(), rim.load()
    wmask = Image.new("L", crop.size, 0)
    wp = wmask.load()
    total = [0, 0, 0]
    count = 0
    for y in range(crop.size[1]):
        for x in range(side):
            r, g, b = cp[x, y]
            if rp[x, y] and min(r, g, b) > white:
                wp[x, y] = 255
            elif ip[x, y] and b > r + 40:
                total[0] += r
                total[1] += g
                total[2] += b
                count += 1
    blue = [round(c / max(1, count)) for c in total]
    tone = wmask.resize((cols, rows), Image.Resampling.BOX)
    area = inside.resize((cols, rows), Image.Resampling.BOX)
    out = []
    for y in range(rows):
        line = ""
        for x in range(cols):
            if area.getpixel((x, y)) < 128:
                line += "."
                continue
            v = min(1.0, max(0.0, (tone.getpixel((x, y)) / 255 - lo) / (hi - lo)))
            line += str(round(v * (levels - 1)))
        out.append(line)
    palette = []
    for i in range(levels):
        t = i / (levels - 1)
        palette.append("#" + "".join(f"{round(c + (255 - c) * t):02x}" for c in blue))
    return out, palette


def as_halfblocks(px: list[str], palette: list[str]) -> list[str]:
    """ANSI truecolour lines of the colour rows (the preview)."""
    def rgb(ch: str) -> str:
        h = palette[int(ch)]
        return f"{int(h[1:3], 16)};{int(h[3:5], 16)};{int(h[5:7], 16)}"
    out = []
    for y in range(0, len(px), 2):
        line = ""
        for x in range(len(px[y])):
            top, low = px[y][x], px[y + 1][x] if y + 1 < len(px) else "."
            if top == "." and low == ".":
                line += "\x1b[0m "
            elif low == ".":
                line += f"\x1b[0m\x1b[38;2;{rgb(top)}m\u2580"
            elif top == ".":
                line += f"\x1b[0m\x1b[38;2;{rgb(low)}m\u2584"
            else:
                line += f"\x1b[38;2;{rgb(top)};48;2;{rgb(low)}m\u2580"
        out.append(line + "\x1b[0m")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cell", choices=("braille", "quadrant", "halfblock"), default="braille")
    ap.add_argument("--cols", type=int, default=None,
                    help="pixel columns (two per braille or quadrant cell, one per halfblock cell); 40, halfblock 26")
    ap.add_argument("--rows", type=int, default=None,
                    help="pixel rows (four per braille cell, two per quadrant or halfblock cell); 32, halfblock 26")
    ap.add_argument("--lo", type=float, default=0.1, help="halfblock: white fraction drawn as the disc's blue")
    ap.add_argument("--hi", type=float, default=0.6, help="halfblock: white fraction drawn as white")
    ap.add_argument("--levels", type=int, default=5, help="halfblock: tone steps from blue to white")
    ap.add_argument("--cover", type=float, default=0.45, help="white fraction of a pixel's area that sets it")
    ap.add_argument("--white", type=int, default=200, help="level all three channels exceed in a white pixel")
    ap.add_argument("--preview", action="store_true", help="also print the cells the window draws")
    ap.add_argument("--logo", type=Path, default=LOGO)
    a = ap.parse_args()
    if a.cell == "halfblock":
        px, palette = color_rows(a.logo, a.cols or 26, a.rows or 26, a.white, a.lo, a.hi, a.levels)
        print("palette:", ", ".join(f'"{c}"' for c in palette))
        for row in px:
            print(f'  "{row}",')
        if a.preview:
            print()
            print("\n".join(as_halfblocks(px, palette)))
        return
    a.cols = a.cols or 40
    a.rows = a.rows or 32
    px = pixel_rows(a.logo, a.cols, a.rows, a.cover, a.white, a.cell)
    for row in px:
        print(f'  "{row}",')
    if a.preview:
        print()
        for line in (as_braille(px) if a.cell == "braille" else as_quadrants(px)):
            print(line.rstrip() + "|")


if __name__ == "__main__":
    main()
