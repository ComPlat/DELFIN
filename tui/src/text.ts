/**
 * Text measured the way the terminal draws it: columns, not characters
 * (wide characters take two, combining marks none -- Bun.stringWidth).
 */

export function width(text: string): number {
  return Bun.stringWidth(text)
}

/** `text` cut to `max` columns, with an ellipsis when cut. */
export function fit(text: string, max: number): string {
  if (max <= 0) return ""
  if (width(text) <= max) return text
  let out = ""
  for (const ch of text) {
    if (width(out + ch) > max - 1) break
    out += ch
  }
  return out + "…"
}

/** Word-wrapped lines of `text`, each at most `max` columns. */
export function wrap(text: string, max: number): string[] {
  const out: string[] = []
  for (const para of text.split("\n")) {
    let line = ""
    for (const word of para.split(/(\s+)/)) {
      if (!word) continue
      if (width(line + word) <= max) {
        line += word
        continue
      }
      if (line.trim()) out.push(line.trimEnd())
      line = word.trimStart()
      while (width(line) > max) {
        const cut = fit(line, max + 1).slice(0, -1)       // without the ellipsis
        out.push(cut)
        line = line.slice(cut.length)
      }
    }
    out.push(line.trimEnd())
  }
  return out
}

