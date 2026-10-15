/**
 * The dolphin at the top of a session: the DELFIN logo (a dolphin leaping
 * over two waves, delfin/logo/DELFIN_logo.png) reduced to braille cells,
 * drawn in the accent -- the logo's blue, or the terminal's own colour
 * without colour (NO_COLOR).
 *
 * The pixels are the output of
 *     python tui/scripts/logo_pixels.py
 * (its defaults: --cell braille --cols 40 --rows 32 --cover 0.45
 * --white 200), pasted; tests/test_the_tui_dolphin_is_drawn_from_the_logo.py
 * runs the script and compares. Four by two they become braille cells
 * (U+2801-28FF), eight rows of 19 columns. Braille rather than quadrant
 * blocks (2x2): at this size the quadrant grid merges the body, the fin
 * and the beak into one blob; braille has twice the rows of pixels. An
 * empty cell is a space, not the blank pattern U+2800, which some fonts
 * draw as faint dots.
 */

const PIXELS = [
  ".............#####.....########.......",
  "...........######################.....",
  ".............####.################....",
  ".............#####################....",
  ".............########################.",
  ".............#########..############..",
  "............#######..############.....",
  "............####...#########..........",
  "...........####..#########............",
  "..........###..#####..####............",
  ".........###..###.....###.............",
  ".........##..##......###..............",
  "........##..##......###...............",
  "........##.#........#.................",
  ".......##.#...........................",
  ".......#.#............................",
  "......###.............................",
  "......#.#.............................",
  "......##..............................",
  "....####..............................",
  "..#######.............................",
  ".#########............................",
  "##.....###.......############.........",
  ".........#.....#########...####.......",
  ".............#######...........#......",
  "...........#######......#####.........",
  ".#......########.....####.............",
  "..#############....####...............",
  "....#########....####.................",
  "................###...................",
  "..............####....................",
  "..........#####.......................",
]

/** The eight dots of a braille cell, as (row, column, bit). */
const DOTS: ReadonlyArray<readonly [number, number, number]> = [
  [0, 0, 0x01], [1, 0, 0x02], [2, 0, 0x04], [0, 1, 0x08],
  [1, 1, 0x10], [2, 1, 0x20], [3, 0, 0x40], [3, 1, 0x80],
]

/** Pixel rows ('#' set) as braille cells, four rows by two columns each. */
export function braille(pixels: readonly string[]): string[] {
  const at = (y: number, x: number) => pixels[y]?.[x] === "#"
  const width = Math.max(...pixels.map((r) => r.length))
  const out: string[] = []
  for (let y = 0; y < pixels.length; y += 4) {
    let line = ""
    for (let x = 0; x < width; x += 2) {
      let bits = 0
      for (const [dy, dx, bit] of DOTS) if (at(y + dy, x + dx)) bits |= bit
      line += bits ? String.fromCodePoint(0x2800 + bits) : " "
    }
    out.push(line)
  }
  return out
}

export const DOLPHIN: readonly string[] = braille(PIXELS)

/** Its width in columns (every row is as wide). */
export const DOLPHIN_WIDTH = DOLPHIN[0]!.length
