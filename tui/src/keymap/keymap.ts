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
  | "agentsPeek" | "agentsBack" | "agentsStop"

export const DEFAULT_KEYS: Readonly<Record<Action, readonly string[]>> = {
  submit: ["return"],
  newline: ["shift+return", "alt+return", "ctrl+j"],
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
  agentsUp: ["up", "k"],
  agentsDown: ["down", "j"],
  agentsOpen: ["right", "return"],
  agentsPeek: ["space"],
  agentsBack: ["escape", "left"],
  agentsStop: ["ctrl+x"],
}

export interface Binding {
  name: string
  ctrl: boolean
  alt: boolean
  shift: boolean
}

const ALIASES: Record<string, string> = {
  esc: "escape", enter: "return", ret: "return", del: "delete", bs: "backspace",
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
