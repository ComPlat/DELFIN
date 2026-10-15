/**
 * B16: ctrl+<letter> and alt+<letter> answer an approval as <letter>.
 * approval.tsx keyChar() falls back to key.name when the sequence is a
 * control character (ctrl+y = "\x19" -> "y"), and handleKey never looks at
 * the modifiers. The REPL reads raw bytes: "\x19" is not an offered key and
 * is ignored. So ctrl+y approves, ctrl+a aborts the turn, ctrl+e turns off
 * write approvals for the session, ctrl+n refuses -- the readline/emacs
 * keys a person presses out of habit in a prompt.
 */
import { afterEach, expect, test } from "bun:test"
import { FIXTURE, mount } from "./h"

let w: Awaited<ReturnType<typeof mount>> | null = null
afterEach(async () => { await w?.close(); w = null })

const edit = FIXTURE.events.find((e: any) => e.type === "permission.asked" && e.kind === "confirm")
const bash = {
  type: "permission.asked", session: "s1", ts: 0, seq: 100, requestId: "c1", kind: "confirm", tool: "bash",
  preview: "$ rm -rf build/ dist/", flags: { shell: true, write: false, protected: false, outsideRead: false },
  options: [{ key: "y", label: "yes" }, { key: "n", label: "no" }, { key: "A", label: "always this exact command" },
            { key: "a", label: "abort" }, { key: "?", label: "explain" }],
}

for (const [label, key, mods, question] of [
  ["ctrl+y", "y", { ctrl: true }, bash],
  ["ctrl+a", "a", { ctrl: true }, bash],
  ["ctrl+n", "n", { ctrl: true }, bash],
  ["alt+y", "y", { meta: true }, bash],
  ["ctrl+e (edit question)", "e", { ctrl: true }, edit],
] as const) {
  test(`${label} with a question open`, async () => {
    w = await mount({ answerPauseMs: 350 })
    w.store.push(question as never)
    await w.settle()
    await Bun.sleep(1000)                                  // a deliberate key, long after it opened
    w.keys.pressKey(key, mods as any)
    await w.settle()
    console.log(label, "->", JSON.stringify(w.replies()))
    expect(w.replies()).toEqual([])
  })
}
