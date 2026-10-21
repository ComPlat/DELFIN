/**
 * A pipe table in an answer, as the window draws it.
 *
 * Input: a markdown block as blocks.ts commits it, and the columns there
 * are. Output: the rows on screen. Reproduction (OpenTUI 0.6.4's markdown
 * table): the alignment row was ignored; at 40-60 columns a number was cut
 * in two ("-613." / "4171289", "-0." / "2515"); a table indented under a
 * list item was drawn as its raw pipes.
 *
 * Semantics pinned here: recognition is the dashboard's and the plain
 * terminal's (delfin/agent/markdown_table.py: the alignment row must have
 * the header's cell count, escaped \| and a pipe in a code span stay in
 * their cell, a column of numbers defaults to the right); the table is
 * drawn in a box with a rule under the header; a column's alignment row
 * is honoured; a number is never split -- text columns are wrapped first,
 * at words, then inside words; a table that cannot fit even so is shown
 * as the model wrote it.
 */
import { describe, expect, test } from "bun:test"
import { createTestRenderer } from "@opentui/core/testing"
import { blockRows } from "../../src/footer/live"
import { drawTable, tableSegments, isNumeric, splitRow, alignment } from "../../src/scrollback/table"
import { commitSurface, markdownBlock } from "../../src/scrollback/writers"
import { DARK, GUTTER, markdownStyles } from "../../src/theme/theme"

const plain = (drawn: ReturnType<typeof drawTable>) =>
  drawn ? drawn.map((line) => line.map((s) => s.text).join("")) : null

const ALIGNED = [
  "| Rechnung | Methode | Energie / Eh | ΔE / kJ mol⁻¹ | Status |",
  "|:---|:---:|---:|---:|:---:|",
  "| run_a | PBE0/def2-SVP | -613.4171289 | 0.00 | ok |",
  "| run_b | B3LYP/def2-TZVP | -613.4012515 | 41.73 | fehlgeschlagen |",
  "| run_c | r2SCAN-3c | -612.9 | 1234.5 | ok |",
]

const NUMBERS = [
  "| Konformer | E / Eh | Besetzung / % | Kommentar zum Konformer und seiner Geometrie |",
  "|---|---|---|---|",
  "| c1 | -0.2515 | 46.36 | der stabilste, mit intramolekularer Wasserstoffbrücke |",
  "| c2 | -10.100800 | 30.4 | zweiter |",
]

async function rendered(md: string, width: number): Promise<string[]> {
  const setup = await createTestRenderer({ width, height: 40, screenMode: "split-footer", footerHeight: 4,
    externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false })
  try {
    await commitSurface(setup.renderer, markdownBlock(md, { styles: markdownStyles(DARK), marginTop: 1, gutter: GUTTER }), 300)
    await setup.flush()
    await setup.renderOnce()
    return setup.externalOutput.takeText().split("\n").map((l) => l.replace(/\s+$/, ""))
  } finally {
    await setup.dispose()
  }
}

describe("recognition, as markdown_table.py has it", () => {
  test("cells split at unescaped pipes outside code spans", () => {
    expect(splitRow("| `a | b` | x \\| y | z |")).toEqual(["`a | b`", "x \\| y", "z"])
  })
  test("the alignment row must have the header's cell count", () => {
    expect(alignment("| a | b |", "|:--|--:|")).toEqual(["l", "r"])
    expect(alignment("| a | b |", "|:-:|")).toBeNull()
    expect(alignment("a | b", "---")).toBeNull()
  })
  test("a number is a number with its sign, separators, exponent and unit-less percent", () => {
    for (const n of ["-613.4171289", "−2515", "1.2e-3", "46.36 %", "1,234.5", "`0.00`"]) expect(isNumeric(n)).toBe(true)
    for (const n of ["def2-SVP", "ok", "41.73 kJ"]) expect(isNumeric(n)).toBe(false)
  })
})

describe("drawing", () => {
  test("the alignment row is honoured: left, centre, right", () => {
    const rows = plain(drawTable(ALIGNED, 100))!
    expect(rows[0]!.startsWith("┌")).toBe(true)
    expect(rows.at(-1)!.startsWith("└")).toBe(true)
    expect(rows[2]!.startsWith("├")).toBe(true)                 // the rule under the header
    const runA = rows.find((r) => r.includes("run_a"))!
    const runB = rows.find((r) => r.includes("run_b"))!
    // Right: the last digits stand in one column.
    expect(runA.indexOf("0.00") + 4).toBe(runB.indexOf("41.73") + 5)
    // Centre: "ok" in the middle of the cell "fehlgeschlagen" is as wide as.
    const cellB = runB.split("│")[5]!
    const cellA = runA.split("│")[5]!
    expect(cellA.length).toBe(cellB.length)
    expect(cellA.trim()).toBe("ok")
    expect(cellA.indexOf("ok")).toBe(Math.floor((cellA.length - 2) / 2))
    // Left.
    expect(runA.split("│")[1]).toBe(" run_a    ")
  })

  test("every row is as wide as the box and no wider than the columns", () => {
    for (const cols of [100, 60, 46, 40]) {
      const rows = plain(drawTable(NUMBERS, cols))
      if (!rows) continue
      const w = Bun.stringWidth(rows[0]!)
      for (const r of rows) expect(Bun.stringWidth(r)).toBe(w)
      expect(w).toBeLessThanOrEqual(cols)
    }
  })

  test("a number is never split: text columns wrap first", () => {
    for (const cols of [60, 50, 44]) {
      for (const table of [ALIGNED, NUMBERS]) {
        const rows = plain(drawTable(table, cols))
        if (!rows) continue
        const text = rows.join("\n")
        for (const n of ["-613.4171289", "-613.4012515", "-0.2515", "-10.100800", "1234.5", "46.36", "41.73"]) {
          if (table.join("").includes(n)) expect(text).toContain(n)
        }
      }
    }
    // At 60 columns the comment column wraps, the numbers stay whole.
    const rows = plain(drawTable(NUMBERS, 60))!
    expect(rows.some((r) => r.includes("-10.100800"))).toBe(true)
    expect(rows.filter((r) => r.includes("│")).length).toBeGreaterThan(4)
  })

  test("a table that cannot fit is not drawn", () => {
    expect(drawTable(ALIGNED, 30)).toBeNull()
  })

  test("inline markup is concealed in a cell, and its width is the text's", () => {
    const rows = plain(drawTable(["| a | b |", "|---|---|", "| **fett** | `x \\| y` |"], 80))!
    expect(rows.find((r) => r.includes("fett"))).toBe("│ fett │ x | y │")
  })
})

describe("segments of a block", () => {
  test("a table indented under a list item is found, with its indent", () => {
    const segs = tableSegments("- Lauf A, die Energien:\n\n  | E / Eh | ΔE |\n  |---|---:|\n  | -1.5 | 0.0 |\n- Lauf B")
    expect(segs.map((s) => s.kind)).toEqual(["markdown", "table", "markdown"])
    expect(segs[1]!.indent).toBe(2)
    expect(segs[1]!.gapBefore).toBe(true)
    expect(segs[2]!.text).toBe("- Lauf B")
  })
  test("a block without a table is one markdown segment, unchanged", () => {
    const md = "Text with a | pipe\n\n```\n| a | b |\n|---|---|\n```"
    expect(tableSegments(md)).toEqual([{ kind: "markdown", text: md, indent: 0, gapBefore: false }])
  })
})

describe("in the window", () => {
  test("a table under a list item is drawn as a table", async () => {
    const out = await rendered("1. Lauf A:\n   | E / Eh | ΔE |\n   |---|---:|\n   | -1.5 | 0.0 |\n2. Lauf B", 60)
    const text = out.join("\n")
    expect(text).toContain("┌")
    expect(text).not.toContain("|---")
    expect(text).toContain("2. Lauf B")
    // Under the item's text, not at the left edge.
    const top = out.find((l) => l.includes("┌"))!
    expect(top.indexOf("┌")).toBeGreaterThan(GUTTER)
  })

  test("at 40 columns no number is cut in two", async () => {
    const text = (await rendered(NUMBERS.join("\n"), 40)).join("\n")
    for (const n of ["-0.2515", "-10.100800", "46.36"]) expect(text).toContain(n)
  })

  test("a table too wide for the window is shown as written", async () => {
    const out = await rendered(ALIGNED.join("\n"), 30)
    expect(out.join("\n")).toContain("|:---|:---:|")
  })

  test("the footer counts the rows the table is drawn in", async () => {
    for (const width of [100, 60, 40]) {
      const out = await rendered(NUMBERS.join("\n"), width)
      const drawn = out.filter((l) => l.trim()).length
      expect(blockRows(NUMBERS.join("\n"), width - GUTTER)).toBe(drawn)
    }
  })
})

describe("streamed into the window", () => {
  const ANSWER = "Ergebnis:\n\n" + ALIGNED.join("\n") + "\n\n1. Lauf A:\n   | E / Eh | ΔE |\n   |---|---:|\n   | -1.5 | 0.0 |\n2. Lauf B\n"
  for (const full of [false, true]) {
    test(`${full ? "full screen" : "inline"}: the rows and the footer draw the tables, not their pipes`, async () => {
      const { mountWindow, FIXTURE } = await import("../harness")
      const w = await mountWindow({ full, width: 100, height: 40 })
      try {
        let seq = 20_000
        const ev = (type: string, fields: Record<string, unknown> = {}) => ({ type, session: "s1", ts: 0, ...fields }) as never
        w.store.push(ev("session.info", { seq: 0, ...FIXTURE.session }))
        w.store.push(ev("user.message", { seq: ++seq, text: "tabelle bitte", source: "user" }))
        w.store.push(ev("turn.started", { seq: ++seq, turnId: "t1", source: "user" }))
        w.store.push(ev("text.started", { partId: "p1" }))
        let shown = ""
        let committed = ""
        const external = (w.t as unknown as { externalOutput?: { takeText(): string } }).externalOutput
        const take = () => {
          const rows = external?.takeText() ?? ""
          committed += rows
          return rows
        }
        // Every frame from the one where the table's alignment row is
        // complete: before it, nothing tells the header line is a table.
        const from = ANSWER.indexOf(ALIGNED[1]!) + ALIGNED[1]!.length + 1
        for (let i = 0; i < ANSWER.length; i += 40) {
          w.store.push(ev("text.delta", { partId: "p1", delta: ANSWER.slice(i, i + 40) }))
          await w.settle()
          const now = take() + "\n" + w.t.captureCharFrame()
          if (i + 40 >= from) shown += now
        }
        w.store.push(ev("text.ended", { partId: "p1", text: ANSWER }))
        w.store.push(ev("turn.ended", { seq: ++seq, turnId: "t1", status: "ok" }))
        await w.settle()
        const end = take() + "\n" + w.t.captureCharFrame()
        shown += end
        expect(shown).not.toContain("|:---|")
        expect(shown).not.toContain("|---|")
        expect(shown).not.toContain("| -1.5 |")
        // What stays: the full screen's last frame, the rows written inline.
        const after = full ? end : committed
        expect(after).toContain("┌")
        // The right-aligned column: its last digits in one column.
        const rows = after.split("\n")
        const a = rows.find((r) => r.includes("run_a") && r.includes("0.00"))!
        const b = rows.find((r) => r.includes("run_b") && r.includes("41.73"))!
        expect(a.indexOf("0.00") + 4).toBe(b.indexOf("41.73") + 5)
        expect(after).toContain("2. Lauf B")
      } finally {
        await w.close()
      }
    })
  }
})
