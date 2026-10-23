/**
 * B11: keypad Enter under the kitty keyboard protocol arrives as
 * "kpenter". The keymap's submit is ["return"] only, so the composer's own
 * Enter logic is skipped and the TEXTAREA's alias (kpenter -> enter ->
 * return -> submit) sends the prompt: with a ctrl+r search open, the match
 * shown goes out at once ("enter takes it into the prompt, not out").
 */
import { afterEach, expect, test } from "bun:test"
import { createTestRenderer } from "@opentui/core/testing"
import { render } from "@opentui/solid"
import { App } from "../../../src/app"
import { History } from "../../../src/composer/history"
import { Keymap } from "../../../src/keymap/keymap"
import { ensureSession } from "../../../src/store/state"
import { createAppStore } from "../../../src/store/store"
import { DARK } from "../../../src/theme/theme"

let setup: Awaited<ReturnType<typeof createTestRenderer>> | null = null
afterEach(async () => { await setup?.renderer.destroy(); setup = null })

test("keypad enter after a trailing backslash", async () => {
  setup = await createTestRenderer({ width: 80, height: 24, screenMode: "split-footer", footerHeight: 4,
    externalOutputMode: "capture-stdout", useMouse: false, exitOnCtrlC: false, kittyKeyboard: true, useKittyKeyboard: {} } as any)
  const store = createAppStore()
  store.update((d) => { const s = ensureSession(d, "s1", { model: "m", approval: "default" }); s.state = "idle"; d.current = "s1" })
  const calls: any[] = []
  const client = { request: async (method: string, params: any = {}) => { calls.push({ method, params }); return { accepted: true } } }
  await render(() => (
    <App renderer={setup!.renderer} store={store} keymap={new Keymap()} theme={() => DARK} client={client}
         history={new History(null, [{ text: "rm -rf build and rerun everything" }])} onQuit={() => {}} onRedraw={() => {}} />
  ), setup.renderer)
  await setup.renderOnce()
  const seen: string[] = []
  setup.renderer.keyInput.on("keypress", (k: any) => seen.push(`${k.name}${k.ctrl ? "+ctrl" : ""}`))
  await setup.mockInput.typeText("first line\\")          // a trailing backslash: enter makes a new line
  ;(setup.renderer as any).stdin.emit("data", Buffer.from("\x1b[57414u"))   // keypad Enter, kitty CSI u
  for (let i = 0; i < 3; i++) { await Bun.sleep(10); await setup.renderOnce() }
  console.log(setup.captureCharFrame().split("\n").filter((l) => l.trim()).join("\n"))
  console.log("keys seen:", seen.join(" "), "\nprompts:", JSON.stringify(calls.filter((c) => c.method === "session.prompt")))
  // Expected (as with the main Enter key): a new line in place of the backslash.
  expect(setup.captureCharFrame()).not.toContain("first line\\")
})
