/**
 * What a tool call looks like in the scrollback: its headline and its
 * result -- the REPL's rules (delfin/agent/repl_render.py), so a session
 * reads the same in both:
 *
 *   ⏺ edit_file  delfin/xtb.py              the call and the argument that matters
 *     ⎿ Edited delfin/xtb.py (+3 -1)        an edit: what it did, then its diff
 *     ⎿ blocked: outside the workspace      a call that did not succeed
 *     ⎿ 34 lines, 2.1 kB                    anything else: its size,
 *         … the last three lines            and where the answer usually is
 *
 * Pure functions of the call's events, so the presenter decides once and
 * the row renderer only draws.
 */

/** The argument a tool's headline shows (repl_render._HEADLINE_ARG). */
export const HEADLINE_ARG: Record<string, readonly string[]> = {
  bash: ["command"], bash_background: ["command"], bash_output: ["job_id"], bash_kill: ["job_id"],
  read_file: ["path", "file_path"], write_file: ["path", "file_path"], edit_file: ["path", "file_path"],
  multi_edit: ["path", "file_path"], notebook_read: ["path", "file_path"], notebook_edit: ["path", "file_path"],
  apply_patch: ["path", "file_path"], grep_file: ["pattern", "query"], list_files: ["path", "pattern"],
  find_definition: ["name", "symbol"], find_references: ["name", "symbol"], read_document: ["path", "file_path"],
  edit_sheet: ["path", "file_path"], create_docx: ["path", "file_path"], create_pdf: ["path", "file_path"],
  draft_email: ["subject", "to"], publish_report: ["title", "path"], search_docs: ["query"],
  read_section: ["section", "doc"], search_calcs: ["query"], get_calc_info: ["calc_id", "name"],
  web_search: ["query"], web_fetch: ["url"], run_tests: ["path", "target"],
  subagent: ["description", "task", "prompt"], orchestrate: ["description", "goal"], skill: ["name"],
  task_create: ["title", "subject"], task_update: ["task_id", "title"], remember: ["text", "name"],
  history_search: ["query"], schedule_wakeup: ["reason"], cron_create: ["name", "schedule"],
  enter_worktree: ["branch", "name"], ask_user_question: ["question"],
}

export const FALLBACK_KEYS = ["command", "path", "file_path", "pattern", "query", "url", "name",
                              "title", "description", "task", "text"] as const

/** Lines of a result shown under its summary. */
export const RESULT_TAIL_LINES = 3
/** Rows of an edit's diff shown; the rest is counted. */
export const EDIT_DIFF_ROWS = 24

const EDIT_TOOLS = new Set(["write_file", "edit_file", "multi_edit", "apply_patch"])
const CARRYING_FIELDS = ["stdout", "output", "text", "result"] as const
const HUNK = /^@@\s+-(\d+)(?:,(\d+))?\s+\+(\d+)(?:,(\d+))?\s+@@/
const PATCH_FILE = /^\*\*\*\s+(?:Update|Add|Delete) File:\s*(.+?)\s*$/

type Fields = Record<string, any>

/** `mcp__delfin-docs__read_file` -> `delfin-docs:read_file`: the server
 *  stays in view (a call that ran outside DELFIN's gate). */
export function shortToolName(name: string): string {
  const n = String(name ?? "")
  if (n.startsWith("mcp__")) {
    const parts = n.split("__")
    if (parts.length >= 3 && parts[1] && parts.slice(2).join("__")) return `${parts[1]}:${parts.slice(2).join("__")}`
  }
  return n
}

export function bareName(name: string): string {
  const short = shortToolName(name)
  return short.includes(":") ? short.split(":").slice(1).join(":") : short
}

import { redact } from "./redact"

/** No control characters (a tool's output is not trusted to steer the terminal). */
export function stripControl(text: string): string {
  return String(text ?? "").replace(/\x1b\[[0-9;?]*[ -/]*[@-~]/g, "").replace(/[\x00-\x08\x0b-\x1f\x7f]/g, "")
}

function oneLine(text: string): string {
  return stripControl(text).replace(/\s+/g, " ").trim()
}

/** The call's name and the argument that matters, on one line. */
export function toolHeadline(name: string, input: Fields | string | null | undefined): { name: string; arg: string } {
  const raw = oneLine(String(name || "tool")) || "tool"
  let args: Fields = {}
  if (typeof input === "string") {
    try {
      const parsed = JSON.parse(input)
      if (parsed && typeof parsed === "object" && !Array.isArray(parsed)) args = parsed
    } catch {
      args = {}
    }
  } else if (input && typeof input === "object") args = input
  for (const key of [...(HEADLINE_ARG[bareName(raw)] ?? []), ...FALLBACK_KEYS]) {
    const value = args[key]
    if ((typeof value === "string" || typeof value === "number") && String(value).trim()) {
      // A command is masked as the REPL masks it (only commands: a path
      // that names "auth" must still say which file is touched).
      const shown = key === "command" ? redact(String(value)) : String(value)
      return { name: shortToolName(raw), arg: oneLine(shown) }
    }
  }
  const n = Object.keys(args).length
  return { name: shortToolName(raw), arg: n ? `(${n} args)` : "" }
}

/** The output a one-line JSON envelope carries ({"exit_code":0,"stdout":…}),
 *  or null: the answer rides in its one carrying field. */
export function envelopeCarrying(output: string): string | null {
  if (output.includes("\n") || !output.trim()) return null
  let payload: unknown
  try {
    payload = JSON.parse(output)
  } catch {
    return null
  }
  if (!payload || typeof payload !== "object" || Array.isArray(payload)) return null
  const values = CARRYING_FIELDS.filter((f) => f in (payload as Fields)).map((f) => (payload as Fields)[f])
  if (values.length !== 1 || typeof values[0] !== "string") return null
  return values[0].trim() ? values[0] : null
}

/** The carrying field of a one-line envelope that was CUT short (the
 *  engine forwards the head of a result): [value, complete], or null. */
export function cutEnvelopeCarrying(output: string): [string, boolean] | null {
  if (output.includes("\n") || !output.trimStart().startsWith("{")) return null
  const keys = [...output.matchAll(new RegExp(`"(${CARRYING_FIELDS.join("|")})"\\s*:\\s*"`, "g"))]
  if (keys.length !== 1) return null
  const start = keys[0]!.index! + keys[0]![0].length
  let i = start
  let complete = false
  while (i < output.length) {
    const ch = output[i]
    if (ch === "\\") {
      i += 2
      continue
    }
    if (ch === '"') {
      complete = true
      break
    }
    i++
  }
  const literal = output.slice(start, Math.min(i, output.length))
  for (let trim = 0; trim <= 6; trim++) {
    try {
      const value = JSON.parse(`"${literal.slice(0, literal.length - trim)}"`) as string
      return value.trim() ? [value, complete] : null
    } catch {
      continue
    }
  }
  return null
}

/** A result as its output, not its packaging (what "expand the last
 *  result" shows). */
export function resultBody(output: string): string {
  const carrying = envelopeCarrying(output ?? "")
  if (carrying !== null) return carrying
  const cut = cutEnvelopeCarrying(output ?? "")
  if (cut) return cut[1] ? cut[0] : cut[0].replace(/\n+$/, "") + "\n… (cut for display)"
  return output ?? ""
}

/** The unified (or apply_patch) part of an editor's result. */
export function extractDiff(text: string): string[] {
  const lines = stripControl(text ?? "").split("\n")
  let start = -1
  let patchForm = false
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i]!
    if (line === "*** Begin Patch") {
      start = i
      patchForm = true
      break
    }
    if (line.startsWith("diff --git ")) {
      start = i
      break
    }
    if (line.startsWith("--- ") && lines[i + 1]?.startsWith("+++ ")) {
      start = i
      break
    }
  }
  if (start < 0) return []
  let end = lines.length
  if (patchForm) {
    for (let i = start + 1; i < lines.length; i++) {
      if (lines[i] === "*** End Patch") {
        end = i + 1
        break
      }
    }
  } else {
    // An editor's advice follows the diff after a truly empty row (an empty
    // context row of a diff is " ", never "").
    for (let i = start + 1; i < lines.length - 1; i++) {
      if (lines[i] === "" && /^\s*(Tip:|Note:|NOTE:)/.test(lines[i + 1]!)) {
        end = i
        break
      }
    }
  }
  return lines.slice(start, end)
}

function diffPath(value: string): string {
  let v = (value ?? "").trim().split("\t")[0]!
  if (!v || v === "/dev/null") return ""
  if (v.startsWith("a/") || v.startsWith("b/")) v = v.slice(2)
  return v
}

export function diffPaths(lines: string[], output: string, input: Fields): string[] {
  const paths: string[] = []
  const keep = (value: string) => {
    const v = diffPath(value)
    if (v && !paths.includes(v)) paths.push(v)
  }
  for (const line of lines) {
    const m = PATCH_FILE.exec(line)
    if (m) keep(m[1]!)
    else if (line.startsWith("+++ ")) keep(line.slice(4))
    else if (line.startsWith("--- ") && paths.length === 0) keep(line.slice(4))
  }
  try {
    const payload = JSON.parse(output ?? "")
    if (payload && typeof payload === "object") for (const p of payload.files_touched ?? []) keep(String(p))
  } catch {
    // not JSON
  }
  for (const key of ["path", "file_path"]) {
    const v = input?.[key]
    if (typeof v === "string" && v.trim()) {
      keep(v)
      break
    }
  }
  return paths
}

export function applyPatchSucceeded(output: string, meta: Fields, input: Fields): boolean {
  if (input?.check_only) return false
  let payload: any = null
  try {
    payload = JSON.parse(output ?? "")
  } catch {
    payload = null
  }
  if (payload && typeof payload === "object" && !Array.isArray(payload)) {
    if (payload.error || payload.check_only) return false
    const status = String(payload.status ?? "").toLowerCase()
    if (["staged", "denied", "error", "failed", "check_failed"].includes(status)) return false
    if (status) return ["ok", "success", "applied"].includes(status)
    return Boolean(payload.files_touched?.length)
  }
  if (/error|failed|refused|denied/.test((output ?? "").toLowerCase())) return false
  return meta?.ok !== false
}

export function humanSize(n: number): string {
  if (n < 1024) return `${n} B`
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(1)} kB`
  return `${(n / 1024 / 1024).toFixed(1)} MB`
}

export type DiffRowKind = "add" | "delete" | "context" | "hunk" | "file" | "meta"
export interface DiffRow {
  kind: DiffRowKind
  text: string
}

/** Diff lines as rows with the line numbers that help (repl_render's). */
export function numberedDiffRows(lines: string[], paths: string[]): DiffRow[] {
  const rows: DiffRow[] = []
  let oldLine: number | null = null
  let newLine: number | null = null
  const multi = paths.length > 1
  for (const line of lines) {
    const pf = PATCH_FILE.exec(line)
    if (pf) {
      oldLine = newLine = null
      if (multi) rows.push({ kind: "file", text: diffPath(pf[1]!) })
      continue
    }
    if (line === "*** Begin Patch" || line === "*** End Patch") continue
    if (line.startsWith("+++ ")) {
      const path = diffPath(line.slice(4))
      if (multi && path) rows.push({ kind: "file", text: path })
      continue
    }
    if (/^(--- |diff --git |index |new file mode |deleted file mode |similarity index |rename from |rename to )/.test(line)) continue
    if (line.startsWith("@@")) {
      const m = HUNK.exec(line)
      oldLine = m ? Number(m[1]) : null
      newLine = m ? Number(m[3]) : null
      rows.push({ kind: "hunk", text: line })
      continue
    }
    if (!line) continue
    if (line.startsWith("\\ No newline")) {
      rows.push({ kind: "meta", text: line })
      continue
    }
    let kind: DiffRowKind = "context"
    let number: number | null = null
    if (line.startsWith("+")) {
      kind = "add"
      number = newLine
      if (newLine !== null) newLine++
    } else if (line.startsWith("-")) {
      kind = "delete"
      number = oldLine
      if (oldLine !== null) oldLine++
    } else if (line.startsWith(" ")) {
      number = newLine ?? oldLine
      if (oldLine !== null) oldLine++
      if (newLine !== null) newLine++
    } else kind = "meta"
    rows.push({ kind, text: number !== null ? `${String(number).padStart(5)} ${line}` : line })
  }
  return rows
}

export type ToolResult =
  | { kind: "blocked"; reason: string }
  | { kind: "edit"; action: "Edited" | "Created" | "Deleted"; label: string; adds: number; dels: number;
      rows: DiffRow[]; omitted: number }
  | { kind: "text"; summary: string; lines: string[]; ok: boolean; exitCode?: number }

/** What came back, and whether it worked (repl_render.tool_result_line). */
export function toolResult(name: string, output: string, meta: Fields = {}, input: Fields = {}): ToolResult {
  const out = stripControl(output ?? "")
  if (meta?.ok === false) {
    return { kind: "blocked", reason: oneLine(String(meta.error ?? "")) || "the call did not succeed" }
  }
  const edit = editResult(name, out, meta, input)
  if (edit) return edit
  if (!out && !Object.keys(meta ?? {}).length) return { kind: "text", summary: "(no output)", lines: [], ok: true }
  // A shell's envelope: its exit code, and the end of what it printed.
  const trimmed = out.trim()
  if (trimmed.startsWith("{") && /"exit_code"\s*:/.test(trimmed) && !trimmed.includes("\n")) {
    try {
      const run = JSON.parse(trimmed)
      const code = Number(run.exit_code ?? 0)
      const printed = (String(run.stdout ?? "") + (run.stderr ? `\n${run.stderr}` : "")).split("\n").filter((l) => l.trim())
      return { kind: "text", summary: `exit ${code}${printed.length ? ` · ${printed.length} line${printed.length === 1 ? "" : "s"}` : ""}`,
               lines: printed.slice(-RESULT_TAIL_LINES).map(oneLine), ok: code === 0, exitCode: code }
    } catch {
      // cut short: handled below
    }
  }
  const carrying = envelopeCarrying(out)
  if (carrying !== null) {
    const lines = carrying.split("\n").filter((l) => l.trim())
    return { kind: "text", summary: lines.length ? `${lines.length} line${lines.length === 1 ? "" : "s"}` : "(no output)",
             lines: lines.slice(-RESULT_TAIL_LINES).map(oneLine), ok: true }
  }
  const cut = cutEnvelopeCarrying(out)
  if (cut) {
    const [value, complete] = cut
    const lines = value.split("\n").filter((l) => l.trim())
    if (complete) {
      return { kind: "text", summary: `${lines.length} line${lines.length === 1 ? "" : "s"}`,
               lines: lines.slice(-RESULT_TAIL_LINES).map(oneLine), ok: true }
    }
    // The output itself was cut: its last lines here are not its last lines.
    const head = lines.slice(0, RESULT_TAIL_LINES).map(oneLine)
    return { kind: "text", summary: `${lines.length}+ lines, cut for display`,
             lines: head, ok: true }
  }
  const lines = out.replace(/\n+$/, "").split("\n")
  const count = out ? lines.length : 0
  const chars = Number(meta?.chars ?? out.length)
  const exact = "chars" in (meta ?? {}) && chars >= out.length
  const truncated = Boolean(meta?.truncated)
  const word = count === 1 ? "line" : "lines"
  let summary = exact && !truncated ? `${count} ${word}, ${humanSize(chars)}` : `${count}+ ${word}, ${humanSize(chars)}`
  if (truncated) summary += " (truncated)"
  const notes = oneLine(String(meta?.notes ?? ""))
  if (notes) summary += ` · ${notes.slice(0, 60)}`
  return { kind: "text", summary, lines: lines.filter((l) => l.trim()).slice(-RESULT_TAIL_LINES).map(oneLine), ok: true }
}

function editResult(name: string, output: string, meta: Fields, input: Fields): ToolResult | null {
  const bare = bareName(name)
  if (!EDIT_TOOLS.has(bare)) return null
  let lines = extractDiff(output)
  if (bare === "apply_patch") {
    if (!applyPatchSucceeded(output, meta, input)) return null
    lines = extractDiff(String(input?.diff ?? input?.patch ?? ""))
  }
  if (!lines.length) return null
  const adds = lines.filter((l) => l.startsWith("+") && !l.startsWith("+++")).length
  const dels = lines.filter((l) => l.startsWith("-") && !l.startsWith("---")).length
  const paths = diffPaths(lines, output, input)
  let action: "Edited" | "Created" | "Deleted" = "Edited"
  if (lines.some((l) => l.startsWith("--- /dev/null") || l.startsWith("*** Add File:"))) action = "Created"
  else if (lines.some((l) => l.startsWith("+++ /dev/null") || l.startsWith("*** Delete File:"))) action = "Deleted"
  const label = paths.length === 1 ? paths[0]! : paths.length ? `${paths.length} files` : "file"
  const all = numberedDiffRows(lines, paths)
  return { kind: "edit", action, label, adds, dels, rows: all.slice(0, EDIT_DIFF_ROWS),
           omitted: Math.max(0, all.length - EDIT_DIFF_ROWS) }
}

export function filetypeOf(path: string): string {
  const ext = path.split(".").pop()?.toLowerCase() ?? ""
  return ({ py: "python", ts: "typescript", tsx: "tsx", js: "javascript", md: "markdown",
            sh: "bash", json: "json", zig: "zig" } as Record<string, string>)[ext] ?? "text"
}
