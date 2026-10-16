/**
 * The approval posture at the bottom: plan
 * and accept edits each in a colour of their own, with a mark, in bold;
 * the unattended bypass in red; default plain. The prompt marker takes the
 * posture's colour. Without colour (mono) the marks still tell them apart.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { TextAttributes, parseColor, rgbToHex, type CapturedSpan } from "@opentui/core"
import { DARK, MONO, postureLook } from "../../src/theme/theme"
import { mountWindow } from "../harness"

let w: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => { await w?.close(); w = null })

const hex = (c: Parameters<typeof rgbToHex>[0]) => rgbToHex(c).toLowerCase()
const same = (a: Parameters<typeof rgbToHex>[0], b: string) => hex(a) === hex(parseColor(b))

async function spansIn(posture: string): Promise<CapturedSpan[]> {
  w = await mountWindow({ width: 100 })
  w.store.update((d) => { d.sessions.s1!.approval = posture })
  await w.settle()
  return w.t.captureSpans().lines.flatMap((l) => l.spans)
}

describe("the posture at the bottom", () => {
  for (const [posture, text, colour] of [["plan", "⏸ plan", DARK.plan], ["acceptEdits", "⏵⏵ accept edits", DARK.edits]] as const) {
    test(`${posture}: its mark and its colour, in bold, and the prompt marker in it`, async () => {
      const spans = await spansIn(posture)
      const label = spans.find((s) => s.text.includes(text))!
      expect(label).toBeDefined()
      expect(same(label.fg, colour as string)).toBe(true)
      expect(label.attributes & TextAttributes.BOLD).toBeTruthy()
      const marker = spans.find((s) => s.text.startsWith("❯"))!
      expect(same(marker.fg, colour as string)).toBe(true)
    })
  }

  test("default: plain, and the marker in the accent", async () => {
    const spans = await spansIn("default")
    const hints = spans.find((s) => s.text.includes("default"))!
    expect(same(hints.fg, DARK.faint as string)).toBe(true)
    expect(spans.some((s) => /⏸|⏵/.test(s.text))).toBe(false)
    const marker = spans.find((s) => s.text.startsWith("❯"))!
    expect(same(marker.fg, DARK.accent as string)).toBe(true)
  })

  test("the bypass is red; without colour the marks stay", () => {
    expect(postureLook("bypassPermissions", DARK).color).toBe(DARK.err)
    for (const posture of ["plan", "acceptEdits", "bypassPermissions"]) {
      const look = postureLook(posture, MONO)
      expect(look.marked).toBe(true)
      expect(look.text).toMatch(/^(⏸|⏵⏵) /)
    }
    expect(postureLook("default", MONO).marked).toBe(false)
  })
})
