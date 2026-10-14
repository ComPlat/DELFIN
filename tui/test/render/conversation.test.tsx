/**
 * A recorded conversation (fixtures/fake-conversation.json, made by
 * scripts/record_fixture.py from a real agent service) rendered headless
 * through the same app and session wiring as the window. The scrollback it
 * commits must be the snapshot (fixtures/fake-conversation.expected.txt);
 * the tmux test (test/e2e/check_tui_tmux.py) runs the same prompts live
 * and must show the same text. `UPDATE_SNAPSHOTS=1 bun test` rewrites it.
 */
import { describe, expect, test } from "bun:test"
import { createTestRenderer } from "@opentui/core/testing"
import { render } from "@opentui/solid"
import { App } from "../../src/app"
import { Keymap } from "../../src/keymap/keymap"
import { startSessionView } from "../../src/session"
import { ensureSession } from "../../src/store/state"
import { createAppStore } from "../../src/store/store"
import { DARK } from "../../src/theme/theme"

const FIXTURES = `${import.meta.dir}/../fixtures`

/** What changes from run to run, made fixed (shared with the tmux test). */
export function normalise(text: string): string {
  return text
    .split("\n")
    .map((line) => line.replace(/\s+$/, ""))
    .map((line) => line.replace(/✓ \d+(\.\d+)?s/g, "✓ N.Ns").replace(/✗ \d+(\.\d+)?s/g, "✗ N.Ns")
      .replace(/↑[\d.]+[kM]? ↓[\d.]+[kM]?/g, "↑N ↓N"))
    .join("\n")
    .replace(/\n{3,}/g, "\n\n")
    .trim()
}

async function renderFixture(width: number): Promise<{ text: string; footer: string }> {
  const fixture = await Bun.file(`${FIXTURES}/fake-conversation.json`).json()
  const setup = await createTestRenderer({
    width, height: 30, screenMode: "split-footer", footerHeight: 4,
    externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false,
  })
  const store = createAppStore()
  const theme = () => DARK
  const view = startSessionView(setup.renderer, store, theme, { showThinking: false, replayTurns: 20 })
  store.update((d) => { ensureSession(d, "s1", fixture.session) })
  view.show("s1")
  await render(() => (
    <App renderer={setup.renderer} store={store} keymap={new Keymap()} theme={theme}
         client={{ request: async () => ({}) }} onQuit={() => {}} onRedraw={() => {}} />
  ), setup.renderer)
  store.push({ type: "session.info", session: "s1", ts: 0, seq: 0, ...fixture.session } as never)
  for (const event of fixture.events) store.push(event)
  store.flush()
  await view.queue.drain()
  await setup.flush()
  await setup.renderOnce()
  const text = setup.externalOutput.takeText()
  const footer = setup.captureCharFrame()
  view.stop()
  await setup.dispose()
  return { text, footer }
}

describe("a recorded conversation", () => {
  test("commits the scrollback the snapshot holds", async () => {
    const { text } = await renderFixture(100)
    const expected = Bun.file(`${FIXTURES}/fake-conversation.expected.txt`)
    if (process.env.UPDATE_SNAPSHOTS === "1" || !(await expected.exists())) {
      await Bun.write(expected, normalise(text) + "\n")
    }
    expect(normalise(text)).toBe((await expected.text()).trim())
  })

  test("shows every part of the conversation, in order", async () => {
    const { text } = await renderFixture(100)
    const order = ["DELFIN  fake/fake-model", "❯ hello there", "⏺ read_file", "The workspace holds", "✓",
                   "❯ please edit the parser", "approved: edit_file delfin/xtb.py", "⏺ edit_file",
                   "match = _ENERGY_RE.search(line)", "regular expression",
                   "❯ ask me which solvent", "answered: Which solvent model?", "Using SMD."]
    let at = 0
    for (const part of order) {
      const found = text.indexOf(part, at)
      expect(found).toBeGreaterThanOrEqual(at)
      at = found
    }
  })

  test("leaves the footer idle, with the composer and the posture", async () => {
    const { footer } = await renderFixture(100)
    expect(footer).toContain("❯ ask DELFIN")
    expect(footer).toContain("default")
    expect(footer).not.toContain("to interrupt")
  })

  test("wraps to a narrow window without losing a part", async () => {
    const { text } = await renderFixture(60)
    for (const part of ["❯ hello there", "⏺ edit_file", "Using SMD."]) expect(text).toContain(part)
    const widest = Math.max(...text.split("\n").map((l) => [...l.replace(/\s+$/, "")].length))
    expect(widest).toBeLessThanOrEqual(60)
  })
})
