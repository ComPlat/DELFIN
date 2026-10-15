/**
 * The dolphin at the top of a session: the DELFIN logo (a dolphin leaping
 * over two waves, delfin/logo/DELFIN_logo.png).
 *
 * Where the terminal draws colours (logoDepthFor): the logo as a picture,
 * the white dolphin and waves on the blue disc, in upper-half blocks --
 * each cell two pixels, the upper one its foreground, the lower one its
 * background -- 26 x 13 cells. Its pixels and palette are the output of
 *     python tui/scripts/logo_pixels.py --cell halfblock
 * pasted (COLOR_PIXELS, COLOR_PALETTE). In 24-bit colour (COLORTERM
 * truecolor or 24bit) the palette as it is; in 256 colours (TERM
 * *256color*) each colour as the nearest of the xterm palette's cube and
 * grey ramp (indices 16-255; 0-15 are the user's to redefine).
 *
 * Otherwise -- NO_COLOR, the mono theme, a terminal that says neither --
 * braille cells drawn in the accent: the logo's blue, or the terminal's
 * own colour without colour.
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

/** The disc's pixels: '.' outside it, else the step into COLOR_PALETTE. */
const COLOR_PIXELS = [
  ".........00000000.........",
  ".......000000000000.......",
  ".....0000000000000000.....",
  "....000000000000000000....",
  "...00000000000000000000...",
  "..0000000000000000000000..",
  "..0000000244412444430000..",
  ".000000001444444444440000.",
  ".000000000444444444444100.",
  "00000000014443344444320000",
  "00000000044234444100000000",
  "00000000441440243000000000",
  "00000002424000440000000000",
  "00000004230001100000000000",
  "00000013200000000000000000",
  "00000033000000000000000000",
  "00001442000000000000000000",
  ".002444420000233310000000.",
  ".000000020144443233200000.",
  "..0000001444400242101000..",
  "..0034444441044200000000..",
  "...00034310244000000000...",
  "....000013440000000000....",
  ".....0000000000000000.....",
  ".......000000000000.......",
  ".........00000000.........",
]

/** The disc's blue to white, one colour per step. */
const COLOR_PALETTE = ["#628ebe", "#89aace", "#b0c6de", "#d8e3ef", "#ffffff"]

/** How the logo is drawn: a picture in 24-bit or 256 colours, or braille. */
export type LogoDepth = "rgb" | "256" | "mono"

/** Which, from the terminal's environment and the theme: braille in the
 *  mono theme (NO_COLOR) whatever the terminal says; DELFIN_TUI_LOGO
 *  (rgb, 256, mono) overrides a terminal that says too much or too little. */
export function logoDepthFor(env: Record<string, string | undefined>, theme: string): LogoDepth {
  if (theme === "mono" || env.NO_COLOR) return "mono"
  const asked = env.DELFIN_TUI_LOGO ?? ""
  if (asked === "rgb" || asked === "256" || asked === "mono") return asked
  const colorterm = (env.COLORTERM ?? "").toLowerCase()
  if (colorterm === "truecolor" || colorterm === "24bit") return "rgb"
  if (/256col/i.test(env.TERM ?? "")) return "256"
  return "mono"
}

let depth: LogoDepth = "mono"

/** Set once by the window at its start (main.tsx); tests leave braille. */
export function setLogoDepth(d: LogoDepth): void {
  depth = d
}

export function logoDepth(): LogoDepth {
  return depth
}

const LEVELS = [0, 95, 135, 175, 215, 255]

/** The xterm 256-colour index nearest to `hex` among the 6x6x6 cube
 *  (16-231) and the grey ramp (232-255), by squared distance -- among those
 *  that keep its hue's order: where two channels differ by more than 16,
 *  the candidate's differ the same way. By distance alone the logo's pale
 *  blue #b0c6de becomes 152 (#afd7d7), a green-grey. */
export function nearest256(hex: string): number {
  const src = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16)) as [number, number, number]
  const keepsHue = (c: readonly number[]) => [[0, 1], [1, 2], [0, 2]].every(([i, j]) => {
    const d = src[j!]! - src[i!]!
    return Math.abs(d) <= 16 || Math.sign(c[j!]! - c[i!]!) === Math.sign(d)
  })
  let index = -1
  let best = Infinity
  for (let i = 16; i < 256; i++) {
    const c = i >= 232 ? [8 + 10 * (i - 232), 8 + 10 * (i - 232), 8 + 10 * (i - 232)]
      : [LEVELS[Math.floor((i - 16) / 36)]!, LEVELS[Math.floor((i - 16) / 6) % 6]!, LEVELS[(i - 16) % 6]!]
    const d = (c[0]! - src[0]) ** 2 + (c[1]! - src[1]) ** 2 + (c[2]! - src[2]) ** 2
    if (d < best && keepsHue(c)) {
      best = d
      index = i
    }
  }
  return index < 0 ? 16 : index
}

/** One run of cells with the same colours: `fg`/`bg` a palette step, or
 *  null for the terminal's own (outside the disc). */
export interface LogoRun {
  text: string
  fg: number | null
  bg: number | null
}

/** The colour rows as runs of cells: ▀ with the upper pixel as foreground
 *  and the lower as background; ▄ where only the lower is on the disc;
 *  a space where neither is. */
export function halfBlocks(pixels: readonly string[]): LogoRun[][] {
  const out: LogoRun[][] = []
  for (let y = 0; y < pixels.length; y += 2) {
    const runs: LogoRun[] = []
    const top = pixels[y]!
    const low = pixels[y + 1] ?? ""
    for (let x = 0; x < top.length; x++) {
      const a = top[x] === "." || top[x] === undefined ? null : Number(top[x])
      const b = low[x] === "." || low[x] === undefined ? null : Number(low[x])
      const cell = a === null && b === null ? { text: " ", fg: null, bg: null }
        : a === null ? { text: "\u2584", fg: b, bg: null }
        : { text: "\u2580", fg: a, bg: b }
      const last = runs[runs.length - 1]
      if (last && last.fg === cell.fg && last.bg === cell.bg && last.text[0] === cell.text) last.text += cell.text
      else runs.push(cell)
    }
    out.push(runs)
  }
  return out
}

export const LOGO_RUNS: readonly LogoRun[][] = halfBlocks(COLOR_PIXELS)

/** Its width in columns. */
export const LOGO_WIDTH = COLOR_PIXELS[0]!.length

/** The palette for a depth: hex colours, or xterm indices. */
export function logoPalette(d: LogoDepth): (string | number)[] {
  return d === "256" ? COLOR_PALETTE.map(nearest256) : [...COLOR_PALETTE]
}
