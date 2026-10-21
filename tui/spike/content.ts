/**
 * Fake session content for the spike, optionally tagged for exact checks.
 *
 * In test mode every line that reaches the scrollback carries a unique tag,
 * "{t00001}", "{t00002}", ...; the spike records which tags it wrote, in
 * order, and the harness compares that list with the tags it finds in the
 * terminal afterwards. A row lost, duplicated or out of place shows up as
 * a difference -- the check phase 0's first harness could not make.
 *
 * Braces, not brackets: markdown reads "[t00001]" as a shortcut link and
 * the renderer conceals the brackets.
 */

export class Tagger {
  private n = 0

  constructor(readonly enabled: boolean) {}

  /** A fresh tag (or "" when tagging is off). */
  next(): string {
    if (!this.enabled) return ""
    this.n += 1
    return `{t${String(this.n).padStart(5, "0")}}`
  }

  /** `text` with a fresh tag appended. */
  line(text: string): string {
    const tag = this.next()
    return tag ? `${text} ${tag}` : text
  }
}

export const TAG_RE = /\{t(\d{5})\}/g

export function tagsIn(text: string): string[] {
  return [...text.matchAll(TAG_RE)].map((m) => m[0])
}

/** A tiny deterministic random source, so a fuzz run can be repeated. */
export function rng(seed: number): () => number {
  let s = seed >>> 0 || 1
  return () => {
    s ^= s << 13
    s ^= s >>> 17
    s ^= s << 5
    return ((s >>> 0) % 1_000_000) / 1_000_000
  }
}

/** A markdown answer: heading, paragraph, list, table, code, closing line. */
export function answer(t: Tagger, turn: number): string {
  return [
    `## What was wrong ${t.next()}`.trimEnd(),
    "",
    t.line("ORCA 6.1 prints the final energy with an extra column,"),
    t.line("so the parser picked up the **wrong field** every time."),
    "",
    `- ${t.line("`delfin/xtb.py` now reads the value with a regex")}`,
    `- ${t.line("the ORCA 5 layout still parses")}`,
    `- ${t.line(`three new tests cover both layouts (turn ${turn})`)}`,
    "",
    `| Output ${t.next()} | Before | After |`,
    "| --- | --- | --- |",
    `| ORCA 5.0.4 ${t.next()} | ok | ok |`,
    `| ORCA 6.1.1 ${t.next()} | wrong field | ok |`,
    "",
    "```python",
    t.line("energy = float(match.group(1))  #"),
    t.line("assert energy < 0  #"),
    "```",
    "",
    t.line("Run the tests to check it yourself."),
    "",
  ].join("\n")
}

/** A unified diff whose lines all carry tags (in a trailing comment). */
export function diff(t: Tagger): string {
  const c = (code: string) => (t.enabled ? `${code}  # ${t.next()}` : code)
  return [
    "--- a/delfin/xtb.py",
    "+++ b/delfin/xtb.py",
    "@@ -41,5 +41,8 @@ def read_final_energy(path):",
    ` ${c("    for line in lines:")}`,
    ` ${c('        if "FINAL SINGLE POINT ENERGY" in line:')}`,
    `-${c("            energy = float(line.split()[3])")}`,
    `+${c("            match = _ENERGY_RE.search(line)")}`,
    `+${c("            if match is None:")}`,
    `+${c('                raise ValueError(f"no energy on line: {line!r}")')}`,
    `+${c("            energy = float(match.group(1))")}`,
    ` ${c("            break")}`,
    ` ${c("    return energy")}`,
    "",
  ].join("\n")
}
