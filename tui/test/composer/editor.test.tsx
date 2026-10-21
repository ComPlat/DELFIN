/**
 * /editor, /keys and /exit in the window: the REPL's rules for writing a
 * prompt in $VISUAL or $EDITOR (a private scratch draft; nothing sent on a
 * failed, empty or unchanged draft), the window's own key table, and
 * leaving with the session that /exit ended.
 */
import { afterEach, describe, expect, test } from "bun:test"
import { chmodSync, mkdtempSync, writeFileSync } from "node:fs"
import { tmpdir } from "node:os"
import { join } from "node:path"
import { editPrompt, editorArgv, splitCommand } from "../../src/composer/editor"
import { mountWindow } from "../harness"

const terminal = () => {
  const calls: string[] = []
  return { calls, suspend: () => { calls.push("suspend") }, resume: () => { calls.push("resume") } }
}

/** An "editor": a script that does `body` to the file it is given ($1). */
function fakeEditor(body: string): string {
  const dir = mkdtempSync(join(tmpdir(), "fake-editor-"))
  const path = join(dir, "edit.sh")
  writeFileSync(path, `#!/bin/sh\n${body}\n`)
  chmodSync(path, 0o755)
  return path
}

describe("writing a prompt in the editor", () => {
  test("the variable is split as a shell would, never handed to one", () => {
    expect(splitCommand("code -w")).toEqual(["code", "-w"])
    expect(splitCommand(`emacs -nw "my file"`)).toEqual(["emacs", "-nw", "my file"])
    expect(splitCommand("vim 'x")).toBe(null)
    expect(editorArgv({ VISUAL: " ", EDITOR: "nano -R" })).toEqual({ argv: ["nano", "-R"], name: "EDITOR" })
    expect(editorArgv({})).toBe(null)
  })

  test("what the editor wrote is sent; the terminal is handed over and back", async () => {
    const t = terminal()
    const editor = fakeEditor(`printf 'line one\\nline two\\n' > "$1"`)
    const out = await editPrompt(t, "", { EDITOR: editor })
    expect(out.text).toBe("line one\nline two")
    expect(out.note).toContain("sending 2 line(s)")
    expect(t.calls).toEqual(["suspend", "resume"])
  })

  test("a failed, an empty and an unchanged draft send nothing", async () => {
    expect((await editPrompt(terminal(), "", { EDITOR: fakeEditor("exit 3") })).note).toContain("exited 3")
    expect((await editPrompt(terminal(), "", { EDITOR: fakeEditor(`: > "$1"`) })).note).toContain("empty")
    const same = await editPrompt(terminal(), "half a thought", { EDITOR: fakeEditor("true") })
    expect(same.text).toBeUndefined()
    expect(same.note).toContain("unchanged")
    expect((await editPrompt(terminal(), "", {})).note).toContain("no editor is configured")
  })

  test("the draft is private and gone afterwards", async () => {
    const editor = fakeEditor(`stat -c %a "$1" > "$1.mode"; cp "$1.mode" /dev/stdout >/dev/null; echo done >> "$1"; echo "$1" > ${tmpdir()}/fake-editor-last`)
    const out = await editPrompt(terminal(), "", { EDITOR: editor })
    expect(out.text).toBe("done")
    const last = (await Bun.file(`${tmpdir()}/fake-editor-last`).text()).trim()
    expect(await Bun.file(last).exists()).toBe(false)
  })
})

let open: Awaited<ReturnType<typeof mountWindow>> | null = null
afterEach(async () => {
  await open?.close()
  open = null
})

describe("the window's own commands", () => {
  test("/keys shows the window's keys, as bound", async () => {
    const w = open = await mountWindow({ width: 120 })
    await w.t.mockInput.typeText("/keys")
    w.t.mockInput.pressEscape()
    await Bun.sleep(40)
    w.t.mockInput.pressEnter()
    await w.settle()
    const text = w.t.externalOutput.takeText()
    expect(text).toContain("keys in this window")
    expect(text).toContain("ctrl+p")
    expect(text).toContain("shift+enter / alt+enter / ctrl+j")
    expect(w.calls.some((c) => c.method === "session.prompt")).toBe(false)
  })

  test("/exit ends the session, and the window leaves with it", async () => {
    const w = open = await mountWindow()
    w.store.push({ type: "session.ended", session: "s1", ts: 0, seq: 1, reason: "exit", exitCode: 0 } as never)
    await w.settle()
    expect(w.quits()).toBe(1)
  })
})
