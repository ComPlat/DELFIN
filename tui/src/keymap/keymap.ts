/**
 * Every key binding in one registry, like OpenCode's keymap.
 *
 * The defaults are today's keys (the REPL's, which follow Claude Code's);
 * `~/.delfin/tui.json` may rebind any action ("keys": {"steer": ["ctrl+g"]}).
 * A binding is a key name with modifiers: "ctrl+g", "shift+tab", "escape",
 * "left", "alt+return". The footer asks `matches("steer", key)`; nothing
 * else in the TUI compares key names.
 */
import type { KeyEvent } from "@opentui/core"

export type Action =
  | "submit" | "newline" | "interrupt" | "steer" | "cycleMode" | "expand"
  | "tasks" | "redraw" | "agents" | "palette" | "quit" | "historyPrev"
  | "historyNext" | "historySearch" | "agentsUp" | "agentsDown" | "agentsOpen"
  | "agentsPeek" | "agentsBack" | "agentsStop" | "agentsRename" | "agentsPin"

export const DEFAULT_KEYS: Readonly<Record<Action, readonly string[]>> = {
  submit: ["return"],
  // ctrl+j arrives as a bare line feed where the terminal does not report
  // modifiers (no kitty keyboard protocol, e.g. tmux without extended keys).
  newline: ["shift+return", "alt+return", "ctrl+j", "linefeed"],
  interrupt: ["escape"],
  steer: ["ctrl+g"],
  cycleMode: ["shift+tab"],
  expand: ["ctrl+o"],
  tasks: ["ctrl+t"],
  redraw: ["ctrl+l"],
  agents: ["left"],                 // with nothing typed
  palette: ["ctrl+p"],
  quit: ["ctrl+d"],                 // with nothing typed
  historyPrev: ["up"],
  historyNext: ["down"],
  historySearch: ["ctrl+r"],
  // In the agent view letters are typing (a new session's prompt): its
  // keys are arrows, space and ctrl keys only.
  agentsUp: ["up"],
  agentsDown: ["down"],
  agentsOpen: ["right", "return"],  // with nothing typed
  agentsPeek: ["space"],            // with nothing typed
  agentsBack: ["escape", "left"],   // left with nothing typed
  agentsStop: ["ctrl+x"],
  agentsRename: ["ctrl+r"],
  agentsPin: ["ctrl+t"],
}

export interface Binding {
  name: string
  ctrl: boolean
  alt: boolean
  shift: boolean
}

const ALIASES: Record<string, string> = {
  esc: "escape", enter: "return", ret: "return", kpenter: "return", del: "delete", bs: "backspace",
  " ": "space", spacebar: "space", pgup: "pageup", pgdn: "pagedown",
}

export function parseBinding(text: string): Binding | null {
  const parts = text.trim().toLowerCase().split("+").map((p) => p.trim()).filter(Boolean)
  if (parts.length === 0) return null
  const binding: Binding = { name: "", ctrl: false, alt: false, shift: false }
  for (const part of parts.slice(0, -1)) {
    if (part === "ctrl" || part === "control") binding.ctrl = true
    else if (part === "alt" || part === "meta" || part === "option") binding.alt = true
    else if (part === "shift") binding.shift = true
    else return null
  }
  const last = parts[parts.length - 1]!
  binding.name = ALIASES[last] ?? last
  return binding
}

/** How a key reads in a hint: short words and arrows. */
const SHOWN: Record<string, string> = {
  escape: "esc", return: "enter", left: "←", right: "→", up: "↑", down: "↓",
}

export function describe(binding: Binding): string {
  return [binding.ctrl && "ctrl", binding.alt && "alt", binding.shift && "shift",
          SHOWN[binding.name] ?? binding.name].filter(Boolean).join("+")
}

export class Keymap {
  private readonly bindings = new Map<Action, Binding[]>()
  readonly problems: string[] = []

  constructor(overrides: Record<string, string[]> = {}) {
    for (const [action, keys] of Object.entries(DEFAULT_KEYS) as [Action, readonly string[]][]) {
      this.bindings.set(action, keys.map((k) => parseBinding(k)!))
    }
    for (const [action, keys] of Object.entries(overrides)) {
      if (!this.bindings.has(action as Action)) {
        this.problems.push(`keys.${action}: no such action`)
        continue
      }
      const parsed = keys.map((k) => [k, parseBinding(k)] as const)
      const bad = parsed.filter(([, b]) => b === null).map(([k]) => k)
      if (bad.length) {
        this.problems.push(`keys.${action}: cannot read ${bad.map((k) => JSON.stringify(k)).join(", ")}`)
        continue
      }
      this.bindings.set(action as Action, parsed.map(([, b]) => b!))
    }
  }

  keys(action: Action): Binding[] {
    return this.bindings.get(action) ?? []
  }

  /** The first binding of `action`, as the hint line shows it. */
  label(action: Action): string {
    const first = this.keys(action)[0]
    return first ? describe(first) : ""
  }

  matches(action: Action, key: Pick<KeyEvent, "name" | "ctrl" | "meta" | "shift"> & { option?: boolean }): boolean {
    const name = ALIASES[key.name] ?? key.name
    const alt = Boolean(key.meta || key.option)
    return this.keys(action).some((b) => b.name === name && b.ctrl === Boolean(key.ctrl)
      && b.alt === alt && b.shift === Boolean(key.shift))
  }
}

/** What each key does, for /keys: the REPL's table (repl_commands._keys),
 *  in its order, then the keys only the TUI has. A row names an action
 *  (its keys come from the keymap, so a rebinding shows) or a fixed key. */
export const KEY_HELP: readonly { action?: Action; keys?: string; does: string }[] = [
  { action: "interrupt", does: "End this turn. A running tool call finishes first. With a question open: refuse it, as n does." },
  { action: "submit", does: "Send what you typed; during a turn it is queued and goes out when the turn ends." },
  { action: "steer", does: "Send what you typed INTO the running turn; the model reads it on its next round." },
  { action: "cycleMode", does: "Next approval posture (never onto bypass)." },
  { action: "expand", does: "Show the last tool result in full." },
  { action: "tasks", does: "Show or hide the open task list at the end of each turn." },
  { action: "redraw", does: "Write the session again at the terminal's width." },
  { keys: "backspace / ctrl+u", does: "Edit or clear the line you are typing (ctrl+w a word, ctrl+k to the end)." },
  { action: "historyPrev", does: "On the prompt's first line: the previous prompt from the history." },
  { action: "historyNext", does: "On the prompt's last line: the next prompt, back to what you were typing." },
  { keys: "tab", does: "Complete a /command or an @path; with a list open, take the selected one." },
  { action: "historySearch", does: "Search the history backwards; again for the next older match; enter takes it." },
  { action: "quit", does: "On an empty prompt: leave the window. The session keeps running in the agent service." },
  { action: "newline", does: "A new line in the prompt (so does a trailing \\ before enter)." },
  { action: "palette", does: "The command palette: commands, and pickers for model, posture, effort, mode, resume." },
  { keys: "ctrl+c", does: "End the turn; else clear what is typed; else, twice, leave the window." },
  { keys: "y n A k R e a ?", does: "With a question open: the answers the session offers (? explains each). N refuses with a reason." },
  { action: "agents", does: "On an empty prompt: the agent view -- every session, grouped by what it needs." },
  { action: "agentsOpen", does: "In the agent view: open the selected session in this window (alt+1..9: the n-th)." },
  { action: "agentsPeek", does: "In the agent view: peek at the selected session; type and enter to reply to it." },
  { action: "agentsStop", does: "In the agent view: stop the selected session (twice); it ends and is saved." },
  { action: "agentsRename", does: "In the agent view: rename the selected session." },
  { action: "agentsPin", does: "In the agent view: pin the selected session to the top." },
  { keys: "typing + enter", does: "In the agent view: a new session, in this window's folder, with what you typed." },
]

/** The key table as lines, with the keys this keymap binds. */
export function keyTable(km: Keymap): string[] {
  const rows = KEY_HELP.map((row) => [row.keys ?? km.keys(row.action!).map(describe).join(" / "), row.does] as const)
  const width = Math.max(...rows.map(([keys]) => keys.length))
  return rows.map(([keys, does]) => `  ${keys.padEnd(width)}  ${does}`)
}
