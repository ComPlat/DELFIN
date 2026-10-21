/**
 * /editor: write the next prompt in $VISUAL or $EDITOR -- the REPL's
 * command (repl_commands._editor), run by the window itself, because the
 * session has no terminal to open an editor in.
 *
 * The REPL's rules: the variable is split the way a shell splits it but
 * never handed to a shell; the draft is a private (0600) scratch file,
 * pre-filled with what followed /editor; a non-zero exit, an empty file,
 * or a file that came back unchanged sends nothing; the file is removed
 * however it ends. The window hands the terminal to the editor
 * (renderer.suspend) and takes it back after.
 */
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs"
import { tmpdir } from "node:os"
import { join } from "node:path"

export interface EditorOutcome {
  /** The prompt to send, when there is one. */
  text?: string
  /** What to say: why nothing was sent, or what is being sent. */
  note: string
}

/** A command line split the way a shell would split it ("code -w",
 *  'emacs -nw'); null when the quotes do not balance. */
export function splitCommand(value: string): string[] | null {
  const out: string[] = []
  let word = ""
  let quote: string | null = null
  let any = false
  for (let i = 0; i < value.length; i++) {
    const ch = value[i]!
    if (quote) {
      if (ch === quote) quote = null
      else if (ch === "\\" && quote === '"' && i + 1 < value.length) word += value[++i]
      else word += ch
      continue
    }
    if (ch === "'" || ch === '"') {
      quote = ch
      any = true
    } else if (ch === "\\" && i + 1 < value.length) {
      word += value[++i]
      any = true
    } else if (/\s/.test(ch)) {
      if (word || any) out.push(word)
      word = ""
      any = false
    } else word += ch
  }
  if (quote) return null
  if (word || any) out.push(word)
  return out
}

export function editorArgv(env: Record<string, string | undefined> = process.env): { argv: string[]; name: string } | null {
  for (const name of ["VISUAL", "EDITOR"]) {
    const value = String(env[name] ?? "").trim()
    if (!value) continue
    const argv = splitCommand(value)
    if (argv && argv.length) return { argv, name }
  }
  return null
}

export interface Terminal {
  suspend(): Promise<void> | void
  resume(): Promise<void> | void
}

/** Run the editor over a scratch draft and decide what came back. */
export async function editPrompt(terminal: Terminal, template: string,
                                 env: Record<string, string | undefined> = process.env): Promise<EditorOutcome> {
  const editor = editorArgv(env)
  if (!editor) {
    return { note: "no editor is configured: set $VISUAL or $EDITOR (e.g. export EDITOR=nano) — nothing was sent" }
  }
  const dir = mkdtempSync(join(env.TMPDIR || tmpdir(), "delfin-prompt-"))
  const path = join(dir, "prompt.md")
  const start = template.trim() ? `${template.trim()}\n` : ""
  try {
    writeFileSync(path, start, { mode: 0o600 })
    await terminal.suspend()
    let code: number
    try {
      const child = Bun.spawn([...editor.argv, path], { stdio: ["inherit", "inherit", "inherit"] })
      code = await child.exited
    } catch (error) {
      return { note: `could not start ${editor.argv[0]} ($${editor.name}): ${(error as Error).message} — nothing was sent` }
    } finally {
      await terminal.resume()
    }
    if (code !== 0) return { note: `${editor.argv[0]} exited ${code} — nothing was sent` }
    let text: string
    try {
      text = readFileSync(path, "utf8")
    } catch (error) {
      return { note: `could not read the draft back: ${(error as Error).message} — nothing was sent` }
    }
    if (!text.trim()) return { note: "the draft is empty — nothing was sent" }
    if (text.trim() === start.trim()) return { note: "the draft is unchanged — nothing was sent" }
    const body = text.trim()
    return { text: body, note: `sending ${body.split("\n").length} line(s) from ${editor.argv[0]}` }
  } finally {
    rmSync(dir, { recursive: true, force: true })
  }
}
