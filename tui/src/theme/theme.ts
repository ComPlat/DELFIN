/**
 * The look: "DELFIN blue, restrained" (docs/TUI_REWRITE.md).
 *
 * Neutral text; the logo's steel blue as the one accent; green, amber and
 * red only for state (done, waiting, failed). Every colour the UI draws
 * comes from a theme's tokens, so a light terminal gets a light theme and
 * a user's `~/.delfin/tui.json` can override any token.
 *
 * "mono" draws in the terminal's own colours (SGR 39/49: no colour of the
 * window's own), for NO_COLOR (no-color.org): bold, italics and underlines
 * still carry the structure.
 */
import { parseColor, RGBA, SyntaxStyle, type ColorInput } from "@opentui/core"

export interface Theme {
  name: "dark" | "light" | "mono"
  accent: ColorInput        // the one accent: prompt marker, rules, selection
  accentSoft: ColorInput    // headings, links, the softer half of the accent
  text: ColorInput          // body text
  bright: ColorInput        // what the user typed, emphasis
  muted: ColorInput         // secondary text: arguments, status details
  faint: ColorInput         // rules, hints, line numbers
  ok: ColorInput            // done
  warn: ColorInput          // waiting for the user, notices
  err: ColorInput           // failed, refused
  code: ColorInput          // inline code, types
  string: ColorInput        // string literals in code
  added: ColorInput         // diff: added line background
  removed: ColorInput       // diff: removed line background
}

export const DARK: Theme = {
  name: "dark",
  accent: "#6290C0", accentSoft: "#8FB2DA", text: "#D5DAE1", bright: "#F1F4F8",
  muted: "#8B94A3", faint: "#555D6B", ok: "#7FB685", warn: "#D9A65B", err: "#D46F6A",
  code: "#C9B98F", string: "#A9C28E", added: "#1D2B22", removed: "#2E1F22",
}

export const LIGHT: Theme = {
  name: "light",
  accent: "#3A6EA5", accentSoft: "#4F7FB3", text: "#2B3038", bright: "#0F1216",
  muted: "#5D6573", faint: "#9AA1AD", ok: "#3E7D47", warn: "#9A6A1C", err: "#B23B36",
  code: "#7A5E1E", string: "#4E6B2E", added: "#E3F1E5", removed: "#F6E3E3",
}

/** The terminal's own foreground; no background of the window's own. */
const own = RGBA.defaultForeground()
const none = RGBA.fromValues(0, 0, 0, 0)

export const MONO: Theme = {
  name: "mono",
  accent: own, accentSoft: own, text: own, bright: own, muted: own, faint: own,
  ok: own, warn: own, err: own, code: own, string: own, added: none, removed: none,
}

/** NO_COLOR (no-color.org): set and not empty. */
export function noColor(env: Record<string, string | undefined> = process.env): boolean {
  return Boolean(env.NO_COLOR)
}

/** A theme by name, with the user's token overrides on top (none in mono:
 *  it has no colours of its own). */
export function themeFor(mode: "dark" | "light" | "mono" | null | undefined,
                         overrides: Partial<Record<keyof Theme, string>> = {}): Theme {
  if (mode === "mono") return MONO
  const base = mode === "light" ? LIGHT : DARK
  const clean: Partial<Theme> = {}
  for (const [key, value] of Object.entries(overrides)) {
    if (key in base && key !== "name" && typeof value === "string" && /^#[0-9a-fA-F]{6}$/.test(value)) {
      ;(clean as Record<string, string>)[key] = value
    }
  }
  return { ...base, ...clean, name: base.name }
}

const heading = (color: ColorInput) => ({ fg: parseColor(color), bold: true })

/** The syntax styles markdown and code are drawn with. */
export function markdownStyles(t: Theme): Parameters<typeof SyntaxStyle.fromStyles>[0] {
  return {
    default: { fg: parseColor(t.text) },
    "markup.heading": heading(t.accentSoft),
    "markup.heading.1": heading(t.accent),
    "markup.heading.2": heading(t.accentSoft),
    "markup.heading.3": heading(t.accentSoft),
    "markup.heading.4": heading(t.text),
    "markup.heading.5": heading(t.text),
    "markup.heading.6": heading(t.text),
    "markup.strong": { fg: parseColor(t.bright), bold: true },
    "markup.italic": { fg: parseColor(t.text), italic: true },
    "markup.list": { fg: parseColor(t.accent) },
    "markup.raw": { fg: parseColor(t.code) },
    "markup.link": { fg: parseColor(t.accentSoft), underline: true },
    "markup.link.label": { fg: parseColor(t.accentSoft), underline: true },
    "markup.link.url": { fg: parseColor(t.muted), underline: true },
    conceal: { fg: parseColor(t.faint) },
    keyword: { fg: parseColor(t.accentSoft), bold: true },
    string: { fg: parseColor(t.string) },
    comment: { fg: parseColor(t.faint), italic: true },
    number: { fg: parseColor(t.warn) },
    function: { fg: parseColor(t.bright) },
    type: { fg: parseColor(t.code) },
    variable: { fg: parseColor(t.text) },
    property: { fg: parseColor(t.text) },
  }
}

export function diffColors(t: Theme) {
  return { text: t.text, lineNumber: t.faint, added: t.added, removed: t.removed,
           addedSign: t.ok, removedSign: t.err }
}

export const SPINNER = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
