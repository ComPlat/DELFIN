/**
 * B04 probe: an approval in a short terminal. footerRows caps the footer at
 * the terminal's height less one; what is drawn when the panel needs more?
 */
import { afterEach, expect, test } from "bun:test"
import { FIXTURE, mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const edit = FIXTURE.events.find((e: any) => e.type === "permission.asked" && e.kind === "confirm")
const upToEdit = FIXTURE.events.slice(0, FIXTURE.events.indexOf(edit) + 1)

for (const height of [24, 18, 14]) {
  test(`an approval at 80x${height}`, async () => {
    w = await mount({ width: 80, height })
    for (const e of upToEdit) w.store.push(e)
    await w.settle()
    const frame = w.frame().split("\n")
    console.log(`--- 80x${height}, footerHeight=${w.t.renderer.footerHeight}`)
    console.log(frame.join("\n"))
    const shown = frame.join("\n")
    expect(shown).toContain("[y] yes")
    expect(shown).toContain("❯")
  })
}
