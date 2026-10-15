/**
 * The dolphin at the top of a session: the DELFIN logo (a dolphin leaping
 * over two waves, delfin/logo/DELFIN_logo.png).
 *
 * Where the terminal draws colours (logoDepthFor): the logo as a picture,
 * the disc with its white dolphin and waves, in upper-half blocks -- each
 * cell two pixels, the upper one its foreground, the lower one its
 * background -- 30 x 15 cells. The pixels are area averages of the PNG
 * (no threshold: a stroke narrower than a pixel keeps its blue-white
 * mix, which is what draws the fin and the waves); outside the disc the
 * terminal's own background. They are the output of
 *     python tui/scripts/logo_pixels.py --cell halfblock
 * pasted: COLOR_PIXELS in 24-bit colour (COLORTERM truecolor or 24bit),
 * COLOR_256 in 256 colours (TERM *256color*), the script's nearest xterm
 * colour (indices 16-255; 0-15 are the user's to redefine) among the
 * greys and the cube colours of the disc's hue.
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

/** The disc's pixels, six hex digits each: '......' outside the disc. */
const COLOR_PIXELS = [
  "..................................................................6894c46792c36692c26390bf6390bf6591c06692c36793c3..................................................................",
  "......................................................6693c36593c46796c86593c56391c16390c06290c06391c26593c56695c76593c46692c2......................................................",
  "..........................................6593c4699acd6492c3618ebd628ebf618fbf628fbe628fbe628fbe628ebe618fbe618dbe608dbd6492c36999cc6492c3..........................................",
  "..............................6591c06898cb6391c2618dbc618dbd628fbf628fbf628fbe628fbe628fbf628ebe628fbe628fbe628fbe628fbe618dbd618ebd6391c26998cb6792c2..............................",
  "........................6491c16999ce608dbd618ebd628ebe618ebd628ebe628fbe628fbe628fbf628fbf628fbf628fbf628fbe628fbf628fbf628fbe628fbf618fbf618fbe6a9ace6491c0........................",
  "..................638ebe6999ce608dbd618ebe628fbe628ebe628ebe628ebe6490bf6490bf6490bf6490bf6390bf6390bf6490bf6490bf628fbf628ebe6390bf6491c0638fbe608dbd6a9ace648fbe..................",
  "..................6999cd618ebd628ebe628fbe628fbe628ebe628ebe5b89bb5585b85888ba5887ba5787ba5e8cbd5b89bb5686b95a89bb5f8cbc618dbd5c8bbc5787ba608dbe638fbe618ebd6998cc..................",
  "............6390c06290c1618ebe628fbe628fbe628ebe628ebd618dbc82a5cab3c8dec6d5e6d0ddeab8cce1789ec67ea3c9aec4ddcedbeadfe7f0e3eaf2d7e2edafc5dd6993c05f8cbd6491bf6391c26490c1............",
  "......6590be699ace618ebc628fbe628fbe628fbe628ebe628ebd638fbe81a4cad1ddebffffffeef3f7d4e0ede2e9f2fcfcfdffffffffffffffffffffffffffffffffffffd7e2ed5e8bbc5686b9628fbe699ace6591c0......",
  "......628ebd6491c2628ebe628ebe618ebd628ebe628fbf628ebe638fbf5182b782a4cadbe4efcedbeaf6f8fafffffffffffffdfcfcf4f6f8f2f4f8f7f8faf8f9fbd6e1ede5ebf3d9e3ee99b6d45e8dbd6593c4638fbe......",
  "......6693c5618ebd628fbe628fbe628ebe618ebe618ebe628ebe648fbf5a88ba9fb9d6e9eff5fffffffffffed3dfecb5cae0bbcee3d1ddebe9eef4ffffffffffffffffffeaeff5d8e3efa3bcd85c8abb618ebd6694c6......",
  "658fbe6695c7628dbe618ebe618ebe628ebe618ebe618ebe638fbf5b89bb769cc5f1f5f8ffffffc0d2e58dadd0aec5dde5ecf3f8f9fafffffefafafcccdae9a6bfda84a7cc6e97c35d8abb5787b9628ebe618ebe6695c86690c0",
  "628ebd6392c3618ebd628ebe628fbe628fbe628ebe628fbf608cbd628ebde9eef5f5f7fa89aacd9bb7d5ecf1f7f0f4f8c2d3e5eff2f7ffffff9db8d64f80b65888ba5a89bb5e8cbd628ebe6490bf628fbe618ebe6492c36490be",
  "628ebd628fc0618ebe628ebe628fbf628fbe628fbf6490bf5484b8c5d5e6f0f3f87a9fc7c5d5e6e0e8f194b2d35485b8779dc6ffffffd8e2ee5e8bbc638fbe6591c06490c06390bf628fbe618ebe628fbe628fbe6390c0648fbd",
  "628ebe628ebe628fbe628ebe628fbe628fbe6390bf5a89bb83a6cbf7f8fa89a9cccedbe9a9c1dc5888ba5686ba6692c0e5ecf3d4e0ed6792c15e8cbd638fbf628fbf628fbf628fbe628fbe628fbe628fbe628fbe628fbe648fbe",
  "628ebe628fbf628fbf628ebe618ebe618ebe638fbf5786b9cedbe9aac1dabbcde184a6cc5183b86490bf5f8dbd88aace95b3d35c8abc5f8cbd638fbf628ebe628fbf628fbe628ebe618ebe628fbe628ebe618ebd628fbf638fbd",
  "6390be6390c1628ebe618ebe618ebe638fbf5c8abb7da1c9c6d6e6acc3db7ba0c85686b96591c0618ebe628fbe5d8bbc5787ba638fbe6390bf628fbe628fbf6290bf628fbe628ebe628ebe628ebe618ebe618dbd6390c16490be",
  "638fbe6593c4618fbe618ebe638fbf6490bf5182b799b5d3c1d1e37fa2c95787ba6591c0628ebe618ebe618ebe628fbe6490bf628fbe628fbe628ebe628fbf628fbe628ebe628ebe628ebe628fbe628fbe618fbe6593c56490be",
  "6690bf6796c8618ebe638fbe5988ba5d8abb779dc6bfd0e396b2d25282b76591c0618ebe628ebe628fbe6490bf6490bf628fbf628ebe628ebe6490bf6491c06390bf628fbe628ebe628fbe628ebe628fbe628ebe6796c86693bf",
  "......6490c1638fbe5b89bb8aabcee7edf4fffffffffefee6ecf47ea1c85b89bb6490bf6490bf5f8dbe5887ba5888ba5c8abc608cbc5e8bbc5988ba5686b95b8abc638fbf638fbf628fbe628ebe628fbe618ebe6391c1......",
  "......648fbe6191c4789ec7d4e0edc7d6e7a5bed998b5d4e3eaf3ebf0f56691bf5e8bbd5787ba6d97c39cb8d6cad8e8e9eef4f1f4f7e4ebf2ccd9e8abc2db7fa3c95d8bbc5e8cbd6390bf628fbe628ebf6695c76490bf......",
  "............6797cb6993c15a89bb5484b85988bb5a89bb5a88ba8babce6792c06691c0afc5ddf1f4f8ffffffffffffcbd9e994b2d37fa3ca88a9cd99b6d4abc2dba2bcd86f97c35e8bbc638fbf618ebd6999cc............",
  "............618cbc6e9aca789dc66490c05888bb5887ba5b8abc6993c1a7c0daeef2f7ffffffffffffb8cbe1648fbe5d8bbc96b3d3b8cce1adc4dd92b0d1799fc77ba0c88aaacd6d96c35f8cbc6897ca6390be............",
  "..................6692c38eafd2c1d2e4cbd9e9cddae9e1e8f1fffffeffffffffffffcfddeb749bc55283b7a6bfdaf1f5f9cedbea89abcf6a95c26490bf6692c0618ebe5e8cbd6793c06291c26593c3..................",
  "........................5f8fc46c97c3a2bcd7d2deece3ebf2d7e2eeb3c8df7a9fc85082b781a5cbe8eef4eaf0f68baccf5585b95989bb618ebe6290bf628fbf628fbf628fbe618fc16897c8........................",
  "..............................6291c35687bc5485b75f8cbc6a94c16e97c28daed0ccdae9e7eef5a9c1db5f8dbd5888bb6491c06490bf628fbf628fbf6290bf618ebe6391c26896c8..............................",
  "....................................6894c46998cb6490bf799fc796b3d2a3bdd98daed06390bf5686ba628fbf6491c0628fbe628fbf628fbf618fbe618ebe6797c96794c5....................................",
  "..........................................6792c26796c96394c85a89bc5686ba5888ba618ebe6491c0628fbf628fbe618fbe628ebf628fbf6695c76999cb6795c4..........................................",
  "......................................................6994c46693c26998c96795c66492c36390c16390c16391c26594c66796c86492c26894c4......................................................",
  "........................................................................6792c26591bf6491be648fbe6592c06794c36a96c5..................................................................",
]

/** The same pixels as xterm 256-colour indices, three digits each: '...'
 *  outside the disc. */
const COLOR_256 = [
  ".................................067067067067067067067067.................................",
  "...........................067067067067067067067067067067067067...........................",
  ".....................067067067067067067067067067067067067067067067067.....................",
  "...............067067067067067067067067067067067067067067067067067067067067...............",
  "............067067067067067067067067067067067067067067067067067067067067067067............",
  ".........067067067067067067067067067067067067067067067067067067067067067067067067.........",
  ".........067067067067067067067067067067067067067067067067067067067067067067067067.........",
  "......067067067067067067067067110251252253251110110251253254255254251067067067067067......",
  "...067067067067067067067067067110253231255253254231231231231231231231254067067067067067...",
  "...067067067067067067067067067067110254253231231231231255255231231254255254110067067067...",
  "...067067067067067067067067067067110255231231253251252253255231231231255254250067067067...",
  "067067067067067067067067067067110255231252110251255231231231188250110067067067067067067067",
  "067067067067067067067067067067255231110110255255252255231110067067067067067067067067067067",
  "067067067067067067067067067252255110252254110067110231254067067067067067067067067067067067",
  "067067067067067067067067110231110253250067067067255253067067067067067067067067067067067067",
  "067067067067067067067067253250251110067067067110110067067067067067067067067067067067067067",
  "067067067067067067067110188250110067067067067067067067067067067067067067067067067067067067",
  "067067067067067067067110252110067067067067067067067067067067067067067067067067067067067067",
  "067067067067067067110252110067067067067067067067067067067067067067067067067067067067067067",
  "...067067067110255231231255110067067067067067067067067067067067067067067067067067067067...",
  "...067067110253188250110255255067067067067110188255255255188250110067067067067067067067...",
  "......067067067067067067067110067067251255231231188110110110110250250067067067067067......",
  "......067067110067067067067067250255231231251067067110251251110110110110067067067067......",
  ".........067110252188188254231231231253110067250255253110067067067067067067067067.........",
  "............067067250253255254251110067110255255110067067067067067067067067067............",
  "...............067067067067067067110188255250067067067067067067067067067067...............",
  "..................067067067110110250110067067067067067067067067067067067..................",
  ".....................067067067067067067067067067067067067067067067067.....................",
  "...........................067067067067067067067067067067067067...........................",
  "....................................067067067067067067067.................................",
]

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

/** One run of cells with the same colours: `fg`/`bg` a hex colour or an
 *  xterm index, or null for the terminal's own (outside the disc). */
export interface LogoRun {
  text: string
  fg: string | number | null
  bg: string | number | null
}

/** Rows of pixel colours as runs of cells: ▀ with the upper pixel as
 *  foreground and the lower as background; ▄ where only the lower is on
 *  the disc; a space where neither is. */
export function halfBlocks(pixels: readonly (readonly (string | number | null)[])[]): LogoRun[][] {
  const out: LogoRun[][] = []
  for (let y = 0; y < pixels.length; y += 2) {
    const runs: LogoRun[] = []
    const top = pixels[y]!
    const low = pixels[y + 1] ?? []
    for (let x = 0; x < top.length; x++) {
      const a = top[x] ?? null
      const b = low[x] ?? null
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

/** A pasted row cut into pixels `width` characters wide: all dots is
 *  outside the disc (null). */
function cells<T>(row: string, width: number, read: (s: string) => T): (T | null)[] {
  const out: (T | null)[] = []
  for (let i = 0; i < row.length; i += width) {
    const s = row.slice(i, i + width)
    out.push(/^\.+$/.test(s) ? null : read(s))
  }
  return out
}

/** The disc as runs: hex colours for 24-bit, xterm indices for 256. */
export const LOGO_RUNS_RGB: readonly LogoRun[][] = halfBlocks(COLOR_PIXELS.map((r) => cells(r, 6, (s) => "#" + s)))
export const LOGO_RUNS_256: readonly LogoRun[][] = halfBlocks(COLOR_256.map((r) => cells(r, 3, Number)))

/** Its width in columns, and its rows. */
export const LOGO_WIDTH = COLOR_PIXELS[0]!.length / 6
export const LOGO_ROWS = LOGO_RUNS_RGB.length

/** The runs for a depth (mono has none: it draws DOLPHIN). */
export function logoRuns(d: LogoDepth): readonly LogoRun[][] {
  return d === "256" ? LOGO_RUNS_256 : d === "rgb" ? LOGO_RUNS_RGB : []
}
