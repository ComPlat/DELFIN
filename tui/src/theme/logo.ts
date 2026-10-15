/**
 * The dolphin at the top of a session, as Claude Code shows its own mark:
 * a small one, diving to the right -- tail up behind, fin, its back
 * curving down to the eye and the beak -- drawn in the accent (the
 * terminal's own colour without colour, NO_COLOR).
 *
 * The pixels are the drawing; two by two they become quadrant blocks
 * (U+2596-259F, which every terminal font the window supports has), six
 * rows of 17 columns.
 */

const PIXELS = [
  "..........##......................",
  "...........###....................",
  ".##.......#######.................",
  "..##..###############.............",
  "...####################...........",
  "..#....###################........",
  "............################......",
  ".................############.....",
  ".....................#####.##.....",
  "........................######....",
  "...........................#####..",
  "..............................###.",
]

const QUADRANTS = " ▘▝▀▖▌▞▛▗▚▐▜▄▙▟█"

/** Pixel rows ('#' set) as block characters, two by two. */
export function quadrants(pixels: readonly string[]): string[] {
  const at = (y: number, x: number) => (pixels[y]?.[x] === "#" ? 1 : 0)
  const width = Math.max(...pixels.map((r) => r.length))
  const out: string[] = []
  for (let y = 0; y < pixels.length; y += 2) {
    let line = ""
    for (let x = 0; x < width; x += 2) {
      line += QUADRANTS[at(y, x) | (at(y, x + 1) << 1) | (at(y + 1, x) << 2) | (at(y + 1, x + 1) << 3)]
    }
    out.push(line)
  }
  return out
}

export const DOLPHIN: readonly string[] = quadrants(PIXELS)

/** Its width in columns (every row is as wide). */
export const DOLPHIN_WIDTH = DOLPHIN[0]!.length
