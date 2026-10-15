/**
 * The dolphin at the top of a session: the DELFIN logo's dolphin and its
 * two waves (delfin/logo/DELFIN_logo.png), without the disc.
 *
 * Where the terminal draws colours (logoDepthFor): the mark in the logo's
 * blue (LOGO_RGB in 24-bit colour, COLORTERM truecolor or 24bit;
 * LOGO_256 in 256 colours, TERM *256color*) on the terminal's own
 * background, in one of three kinds of cell (logoCellsFor):
 *   sextant    2 x 3 sub-pixels a cell (U+1FB00-1FB3B, the halves and
 *              the full block); the default
 *   octant     2 x 4 (U+1CD00-1CDE5 and the block elements that hold the
 *              other 26 patterns); DELFIN_TUI_LOGO=octant
 *   halfblock  1 x 2 (U+2580, U+2584, U+2588); DELFIN_TUI_LOGO=halfblock
 * Sextants are Unicode 13 (2020), octants Unicode 16 (2024): a terminal
 * or font without them draws boxes or question marks (octants: the
 * VS Code terminal); half blocks are in every font with block elements.
 * Sextants for every terminal rather than a guess from TERM_PROGRAM,
 * KITTY_WINDOW_ID or TERM: the environment names the terminal, not its
 * font, and inside tmux or over ssh it names neither.
 *
 * In three sizes (pickLogo), each in every kind of cell: L 30 x 13
 * cells, M 20 x 9, S 16 x 7. The banner takes the largest whose minimum
 * of terminal rows is met and whose row (mark, LOGO_GAP, BANNER_TEXT_MIN
 * for the words) fits the width; the words alone below S. The minimum
 * rows keep the mark at about a third of the terminal's height or less.
 * OCTANT_PIXELS_<cols>, SEXTANT_PIXELS_<cols> and HALF_PIXELS_<cols>
 * ('#' set) and the colour are the output of
 *     python tui/scripts/logo_pixels.py --cell octant|sextant|halfblock --cols 30|20|16
 * pasted (the method is in the script); the tests run it and compare.
 *
 * Otherwise -- NO_COLOR, the mono theme, a terminal that says neither --
 * braille cells drawn in the accent: the logo's blue, or the terminal's
 * own colour without colour.
 *
 * The braille pixels (PIXELS) are the output of
 *     python tui/scripts/logo_pixels.py
 * (its defaults: --cell braille --cols 40 --rows 32 --cover 0.45
 * --white 200), pasted; tests/test_the_tui_dolphin_is_drawn_from_the_logo.py
 * runs the script and compares. Four by two they become braille cells
 * (U+2801-28FF), eight rows of 19 columns, drawn where the terminal has
 * BRAILLE_MIN_ROWS rows (8 rows, a third of 24). Braille rather than quadrant
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

/** The mark's colour: the median of the disc's blue pixels in the PNG,
 *  and the nearest xterm colour of its hue (the script's COLOUR line). */
export const LOGO_RGB = "#628ebe"
export const LOGO_256 = 67

/** The mark in octant sub-pixels, 60 x 52 (30 x 13 cells). */
const OCTANT_PIXELS_30 = [
  ".........................##.............##########..........",
  "....................############....################........",
  "..................#############..#####################......",
  ".....................########..########################.....",
  "......................######.##########################.....",
  "......................####..########################.###....",
  "......................###.##########################.######.",
  "......................##.###################################",
  ".....................#..###############....#################",
  "......................############.....###################..",
  ".....................##########.....##################......",
  "....................#########....################...........",
  "...................########....##############...............",
  "...................######....##############.................",
  "..................######...################.................",
  ".................#####....######....######..................",
  "................#####....#####......######..................",
  "................####...#####........#####...................",
  "...............####...####.........#####....................",
  "...............###...###...........####.....................",
  "..............###...###...........####......................",
  "..............###..###...........###........................",
  ".............###..##............#...........................",
  ".............##..##.........................................",
  "............###.##..........................................",
  "............##..#...........................................",
  "............#..#............................................",
  "...........##.##............................................",
  "...........#..#.............................................",
  "...........#.#..............................................",
  "..........####..............................................",
  ".........##.#...............................................",
  "......########..............................................",
  "....###########.............................................",
  "...#############............................................",
  "..########.######................#########..................",
  ".####.......#####............#################..............",
  "...............###.........######################...........",
  ".................#.......############..........####.........",
  ".......................###########................###.......",
  ".....................###########....................##......",
  "...................###########.........##########.....#.....",
  ".#...............###########........#######.................",
  "...###.......##############.......######....................",
  "....#####################.......######......................",
  "......#################.......######........................",
  "........#############........#####..........................",
  "............#####..........######...........................",
  "..........................#####.............................",
  "........................######..............................",
  ".....................######.................................",
  "................########....................................",
]

/** The mark in octant sub-pixels, 40 x 36 (20 x 9 cells). */
const OCTANT_PIXELS_20 = [
  "...............####.......########......",
  "............#########..#############....",
  "..............######################....",
  "...............###.##################...",
  "...............########################.",
  "...............#########################",
  "...............########...#############.",
  "..............######...############.....",
  ".............#####...##########.........",
  "............####...##########...........",
  "............###...####..####............",
  "...........###..####....####............",
  "..........###..###......###.............",
  "..........##..##.......###..............",
  ".........##..##.......###...............",
  ".........##.##........#.................",
  "........##.##...........................",
  "........#..#............................",
  "........#.#.............................",
  ".......#.#..............................",
  ".......#................................",
  ".......##...............................",
  "....#####...............................",
  "..#########.............................",
  ".##########...........######............",
  ".#.......###.......#############........",
  "...........#.....#########.....###......",
  "...............#######............#.....",
  ".............#######.......#####........",
  "...........########.....####............",
  "..###############.....####..............",
  "....###########.....####................",
  "......#######......####.................",
  ".................####...................",
  "...............####.....................",
  "...........######.......................",
]

/** The mark in octant sub-pixels, 32 x 28 (16 x 7 cells). */
const OCTANT_PIXELS_16 = [
  "............####....#######.....",
  "..........###################...",
  "............##################..",
  "............###################.",
  "............####################",
  "...........######..##########...",
  "...........####..########.......",
  "..........###..########.........",
  ".........###..###..###..........",
  "........###.###....###..........",
  "........##.##.....###...........",
  ".......##.##......#.............",
  ".......#.#......................",
  "......##........................",
  "......#.#.......................",
  "......##........................",
  ".....##.........................",
  "...####.........................",
  ".########.......................",
  ".###..###.......########........",
  "..............#######...###.....",
  "............######.........#....",
  "..........######....#####.......",
  "..##..########....###...........",
  "...#########....###.............",
  "......####.....###..............",
  ".............###................",
  "..........####..................",
]

/** The mark in sextant sub-pixels, 60 x 39 (30 x 13 cells). */
const SEXTANT_PIXELS_30 = [
  "........................####............##########..........",
  "...................#############...##################.......",
  ".....................#########..######################......",
  "......................#####..##########################.....",
  "......................###..#########################.###....",
  "......................##.###################################",
  ".....................#.##############.....#################.",
  "......................##########.....###################....",
  ".....................########.....################..........",
  "...................#######.....#############................",
  "..................######....###############.................",
  ".................#####....#######...######..................",
  "................#####...#####.......######..................",
  "................####...####........######...................",
  "...............###...####..........####.....................",
  "..............###...###...........####......................",
  ".............###...##............##.........................",
  ".............##..##.........................................",
  "............##..##..........................................",
  "............##.##...........................................",
  "...........##.##............................................",
  "...........#..#.............................................",
  "...........#.#..............................................",
  "..........#.#...............................................",
  ".....#########..............................................",
  "...#############............................................",
  "..########.######..................######...................",
  ".##..........####............##################.............",
  "................##........##############.....#####..........",
  ".......................###########................##........",
  ".....................##########......................#......",
  "..................###########.........#######...............",
  "..###.........#############.......######....................",
  "....#####################.......#####.......................",
  ".......###############........#####.........................",
  "...........#######..........#####...........................",
  ".........................######.............................",
  "......................######................................",
  ".................########...................................",
]

/** The mark in sextant sub-pixels, 40 x 27 (20 x 9 cells). */
const SEXTANT_PIXELS_20 = [
  "...............#####.....#########......",
  ".............#######################....",
  "...............######################...",
  "...............########################.",
  "...............##########.##############",
  "..............######....###########.....",
  ".............#####...#########..........",
  "............####..##########............",
  "...........###...###....####............",
  "..........###..##.......###.............",
  "..........##..#........##...............",
  ".........##.##........#.................",
  "........##.#............................",
  "........#.#.............................",
  ".......#.#..............................",
  ".......#................................",
  ".....####...............................",
  "..########..............................",
  ".#####.####...........######............",
  "..........##......###############.......",
  "...............########...........#.....",
  ".............#######......#####.........",
  "..#.......########.....####.............",
  "....############.....###................",
  ".......#####.......###..................",
  "................####....................",
  ".............####.......................",
]

/** The mark in sextant sub-pixels, 32 x 21 (16 x 7 cells). */
const SEXTANT_PIXELS_16 = [
  "...........#####....#######.....",
  "...........##################...",
  "............##################..",
  "............####################",
  "...........#####..#########.....",
  "..........###...#######.........",
  ".........##..###...###..........",
  "........##..##.....##...........",
  ".......##.##......##............",
  ".......#.#......................",
  "......#.#.......................",
  "......##........................",
  ".....##.........................",
  "..######........................",
  ".###..###........######.........",
  ".............#######.....##.....",
  "...........######...............",
  "........######....###...........",
  "....########....###.............",
  "..............###...............",
  "...........###..................",
]

/** The mark in half-block sub-pixels, 30 x 26 (30 x 13 cells). */
const HALF_PIXELS_30 = [
  "...........####....#######....",
  "..........#################...",
  "...........#################..",
  "...........###################",
  "...........#######.##########.",
  "..........#####..#########....",
  "..........###..#######........",
  ".........##..####.###.........",
  "........##..##....###.........",
  ".......##..##....###..........",
  ".......##.#......##...........",
  "......##.#....................",
  "......#.#.....................",
  "......##......................",
  ".....#........................",
  ".....##.......................",
  "..#####.......................",
  ".#######......................",
  ".......##.....##########......",
  "............######......##....",
  "..........#####......##.......",
  "........######...####.........",
  "..##########....##............",
  ".....#####....###.............",
  "............###...............",
  ".........####.................",
]

/** The mark in half-block sub-pixels, 20 x 18 (20 x 9 cells). */
const HALF_PIXELS_20 = [
  ".......###..#####...",
  ".......###########..",
  ".......#############",
  ".......####.#######.",
  "......###.#####.....",
  "......#.##..##......",
  ".....#.##...#.......",
  ".....##....#........",
  "....#...............",
  "....#...............",
  "....................",
  ".####...............",
  ".##.##....#####.....",
  "........####........",
  "......####...##.....",
  "..######...#........",
  "....#....##.........",
  ".......##...........",
]

/** The mark in half-block sub-pixels, 16 x 14 (16 x 7 cells). */
const HALF_PIXELS_16 = [
  ".....###.#####..",
  "......#########.",
  "......#########.",
  ".....##.####....",
  "....####..#.....",
  "....##...#......",
  "...#............",
  "...#............",
  "..##............",
  ".###.....##.....",
  "......####......",
  "....###..##.....",
  "..####..#.......",
  "......#.........",
]

/** How the logo is drawn: the mark in 24-bit or 256 colours, or braille. */
export type LogoDepth = "rgb" | "256" | "mono"

/** Which, from the terminal's environment and the theme: braille in the
 *  mono theme (NO_COLOR) whatever the terminal says; DELFIN_TUI_LOGO
 *  (rgb, 256, mono) overrides a terminal that says too much or too little;
 *  its cell kinds (sextant, octant, halfblock) leave the depth to the
 *  terminal (logoCellsFor). */
export function logoDepthFor(env: Record<string, string | undefined>, theme: string): LogoDepth {
  if (theme === "mono" || env.NO_COLOR) return "mono"
  const asked = env.DELFIN_TUI_LOGO ?? ""
  if (asked === "rgb" || asked === "256" || asked === "mono") return asked
  const colorterm = (env.COLORTERM ?? "").toLowerCase()
  if (colorterm === "truecolor" || colorterm === "24bit") return "rgb"
  if (/256col/i.test(env.TERM ?? "")) return "256"
  return "mono"
}

/** The kind of cell the mark is drawn in (see the top of the file). */
export type LogoCells = "octant" | "sextant" | "halfblock"

/** Sextants unless DELFIN_TUI_LOGO asks for octants (sharper, where the
 *  font has them) or half blocks (a terminal that shows sextants as
 *  boxes or question marks). */
export function logoCellsFor(env: Record<string, string | undefined>): LogoCells {
  const asked = env.DELFIN_TUI_LOGO ?? ""
  return asked === "octant" || asked === "halfblock" ? asked : "sextant"
}

let depth: LogoDepth = "mono"
let cellKind: LogoCells = "sextant"

/** Set once by the window at its start (main.tsx); tests leave braille. */
export function setLogoDepth(d: LogoDepth): void {
  depth = d
}

export function logoDepth(): LogoDepth {
  return depth
}

/** Set once by the window at its start, beside the depth. */
export function setLogoCells(c: LogoCells): void {
  cellKind = c
}

export function logoCells(): LogoCells {
  return cellKind
}

/** Octant patterns (bit n-1 for octant n; 1-2 the top row, 7-8 the
 *  bottom) that are not in U+1CD00-1CDE5: the space, the full, half and
 *  quadrant blocks, and the eighth-row blocks the Unicode octant names
 *  leave them to. */
const OCTANT_ELSEWHERE: ReadonlyMap<number, string> = new Map([
  [0x00, " "], [0xff, "\u2588"], [0x0f, "\u2580"], [0xf0, "\u2584"], [0x55, "\u258c"], [0xaa, "\u2590"],
  [0x05, "\u2598"], [0x0a, "\u259d"], [0x50, "\u2596"], [0xa0, "\u2597"], [0xa5, "\u259a"], [0x5a, "\u259e"],
  [0x5f, "\u259b"], [0xaf, "\u259c"], [0xf5, "\u2599"], [0xfa, "\u259f"],
  [0x01, "\u{1cea8}"], [0x02, "\u{1ceab}"], [0x40, "\u{1cea3}"], [0x80, "\u{1cea0}"],
  [0x03, "\u{1fb82}"], [0x3f, "\u{1fb85}"], [0xc0, "\u2582"], [0xfc, "\u2586"],
  [0x14, "\u{1fbe6}"], [0x28, "\u{1fbe7}"],
])

/** The glyph of an octant pattern: BLOCK OCTANT-<digits> (U+1CD00-1CDE5)
 *  holds the 230 others in ascending order of the pattern. */
export function octant(bits: number): string {
  const other = OCTANT_ELSEWHERE.get(bits)
  if (other !== undefined) return other
  let rank = 0
  for (let b = 0; b < bits; b++) if (!OCTANT_ELSEWHERE.has(b)) rank++
  return String.fromCodePoint(0x1cd00 + rank)
}

/** The glyph of a sextant pattern (bit 1 upper left ... 32 lower right):
 *  U+1FB00-1FB3B hold 60 in order; empty, left half, right half and
 *  full are the space, U+258C, U+2590 and U+2588. */
export function sextant(bits: number): string {
  if (bits === 0) return " "
  if (bits === 21) return "\u258c"
  if (bits === 42) return "\u2590"
  if (bits === 63) return "\u2588"
  return String.fromCodePoint(0x1fb00 + bits - 1 - (bits > 21 ? 1 : 0) - (bits > 42 ? 1 : 0))
}

/** The glyph of a half-block pattern (1 upper, 2 lower). */
export function halfBlock(bits: number): string {
  return [" ", "\u2580", "\u2584", "\u2588"][bits]!
}

/** Sub-pixel rows ('#' set) as lines of cells `sx` x `sy` sub-pixels each. */
export function cellLines(pixels: readonly string[], sx: number, sy: number, glyph: (bits: number) => string): string[] {
  const out: string[] = []
  for (let y = 0; y < pixels.length; y += sy) {
    let line = ""
    for (let x = 0; x < pixels[0]!.length; x += sx) {
      let bits = 0
      for (let dy = 0; dy < sy; dy++) for (let dx = 0; dx < sx; dx++) if (pixels[y + dy]?.[x + dx] === "#") bits |= 1 << (dy * sx + dx)
      line += glyph(bits)
    }
    out.push(line)
  }
  return out
}

/** The mark's sizes, largest first: its width in columns, and the fewest
 *  terminal rows it is drawn in (its rows, 13, 9 and 7 in every kind of
 *  cell, are then at most about a third: 13/40, 9/28, 7/20). */
export type LogoSize = "L" | "M" | "S"
export const LOGO_SIZES: ReadonlyArray<{ size: LogoSize; cols: number; minRows: number }> = [
  { size: "L", cols: 30, minRows: 40 },
  { size: "M", cols: 20, minRows: 28 },
  { size: "S", cols: 16, minRows: 20 },
]

/** The mark's lines per kind of cell and size. */
export const MARK: Readonly<Record<LogoCells, Readonly<Record<LogoSize, readonly string[]>>>> = {
  octant: {
    L: cellLines(OCTANT_PIXELS_30, 2, 4, octant),
    M: cellLines(OCTANT_PIXELS_20, 2, 4, octant),
    S: cellLines(OCTANT_PIXELS_16, 2, 4, octant),
  },
  sextant: {
    L: cellLines(SEXTANT_PIXELS_30, 2, 3, sextant),
    M: cellLines(SEXTANT_PIXELS_20, 2, 3, sextant),
    S: cellLines(SEXTANT_PIXELS_16, 2, 3, sextant),
  },
  halfblock: {
    L: cellLines(HALF_PIXELS_30, 1, 2, halfBlock),
    M: cellLines(HALF_PIXELS_20, 1, 2, halfBlock),
    S: cellLines(HALF_PIXELS_16, 1, 2, halfBlock),
  },
}

/** The banner: the picture, this much space, then the words, which need
 *  about this many columns beside it. */
export const LOGO_GAP = 3
export const BANNER_TEXT_MIN = 32
/** The fewest terminal rows the braille outline (8 rows) is drawn in. */
export const BRAILLE_MIN_ROWS = 24

/** What the banner draws beside its words: a size of the mark (colour) or
 *  the braille outline (mono) -- the largest whose minimum of terminal
 *  rows `rows` is met and whose row fits `width` (the banner's width: the
 *  window less the gutter) -- or null, the words alone. */
export function pickLogo(d: LogoDepth, width: number, rows: number): LogoSize | "braille" | null {
  const fits = (cols: number) => width >= cols + LOGO_GAP + BANNER_TEXT_MIN
  if (d === "mono") return rows >= BRAILLE_MIN_ROWS && fits(DOLPHIN_WIDTH) ? "braille" : null
  for (const s of LOGO_SIZES) if (rows >= s.minRows && fits(s.cols)) return s.size
  return null
}

/** The width in columns of what pickLogo picked (0 for none). */
export function pickedWidth(p: LogoSize | "braille" | null): number {
  if (p === null) return 0
  if (p === "braille") return DOLPHIN_WIDTH
  return LOGO_SIZES.find((s) => s.size === p)!.cols
}

/** One run of cells with the same colours: `fg`/`bg` a hex colour or an
 *  xterm index, or null for the terminal's own. */
export interface LogoRun {
  text: string
  fg: string | number | null
  bg: string | number | null
}

/** The runs for a depth, size and kind of cell (mono has none: it draws
 *  DOLPHIN): each line one run in the logo's blue on the terminal's
 *  background. */
export function logoRuns(d: LogoDepth, size: LogoSize, c: LogoCells = cellKind): readonly LogoRun[][] {
  if (d === "mono") return []
  const fg = d === "256" ? LOGO_256 : LOGO_RGB
  return MARK[c][size].map((text) => [{ text, fg, bg: null }])
}
