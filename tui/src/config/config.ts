/**
 * The user's TUI settings: `~/.delfin/tui.json` (DELFIN_TUI_CONFIG moves it).
 *
 *   {
 *     "theme": "auto" | "dark" | "light" | "mono",
 *     "colors": { "accent": "#6290C0", ... },       // theme token overrides
 *     "keys": { "steer": ["ctrl+g"], ... },          // keymap overrides
 *     "replayTurns": 20,                              // turns shown on attach
 *     "showThinking": false,
 *     "screen": "auto" | "full" | "inline"           // see screenModeOf()
 *   }
 *
 * A file that cannot be read or does not parse costs its settings, never
 * the start: the defaults apply and the problem is reported once.
 */
import { homedir } from "node:os"
import { join } from "node:path"

export interface Config {
  theme: "auto" | "dark" | "light" | "mono"
  colors: Record<string, string>
  keys: Record<string, string[]>
  replayTurns: number
  showThinking: boolean
  screen: "auto" | "full" | "inline"
}

export const DEFAULTS: Config = {
  theme: "auto",
  colors: {},
  keys: {},
  replayTurns: 20,
  showThinking: false,
  screen: "auto",
}

export function configPath(env: Record<string, string | undefined> = process.env): string {
  return env.DELFIN_TUI_CONFIG || join(env.HOME || homedir(), ".delfin", "tui.json")
}

export interface Loaded {
  config: Config
  problems: string[]
}

/** Parse a config text; anything unknown or malformed is reported and dropped. */
export function parseConfig(text: string): Loaded {
  const problems: string[] = []
  let raw: unknown
  try {
    raw = JSON.parse(text)
  } catch (error) {
    return { config: { ...DEFAULTS }, problems: [`not JSON: ${(error as Error).message}`] }
  }
  if (!raw || typeof raw !== "object" || Array.isArray(raw)) {
    return { config: { ...DEFAULTS }, problems: ["not a JSON object"] }
  }
  const src = raw as Record<string, unknown>
  const config: Config = { ...DEFAULTS, colors: {}, keys: {} }
  for (const [key, value] of Object.entries(src)) {
    switch (key) {
      case "theme":
        if (value === "auto" || value === "dark" || value === "light" || value === "mono") config.theme = value
        else problems.push(`theme must be "auto", "dark", "light" or "mono"`)
        break
      case "colors":
        if (value && typeof value === "object" && !Array.isArray(value)) {
          for (const [token, hex] of Object.entries(value as Record<string, unknown>)) {
            if (typeof hex === "string" && /^#[0-9a-fA-F]{6}$/.test(hex)) config.colors[token] = hex
            else problems.push(`colors.${token} must be "#rrggbb"`)
          }
        } else problems.push("colors must be an object")
        break
      case "keys":
        if (value && typeof value === "object" && !Array.isArray(value)) {
          for (const [action, keys] of Object.entries(value as Record<string, unknown>)) {
            const list = typeof keys === "string" ? [keys] : keys
            if (Array.isArray(list) && list.every((k) => typeof k === "string")) {
              config.keys[action] = list as string[]
            } else problems.push(`keys.${action} must be a key or a list of keys`)
          }
        } else problems.push("keys must be an object")
        break
      case "replayTurns":
        if (typeof value === "number" && Number.isInteger(value) && value >= 0 && value <= 1000) {
          config.replayTurns = value
        } else problems.push("replayTurns must be a whole number from 0 to 1000")
        break
      case "showThinking":
        if (typeof value === "boolean") config.showThinking = value
        else problems.push("showThinking must be true or false")
        break
      case "screen":
        if (value === "auto" || value === "full" || value === "inline") config.screen = value
        else problems.push(`screen must be "auto", "full" or "inline"`)
        break
      default:
        problems.push(`unknown setting ${JSON.stringify(key)}`)
    }
  }
  return { config, problems }
}

export async function loadConfig(path: string = configPath()): Promise<Loaded> {
  const file = Bun.file(path)
  if (!(await file.exists())) return { config: { ...DEFAULTS }, problems: [] }
  try {
    const loaded = parseConfig(await file.text())
    return { config: loaded.config, problems: loaded.problems.map((p) => `${path}: ${p}`) }
  } catch (error) {
    return { config: { ...DEFAULTS }, problems: [`${path}: ${(error as Error).message}`] }
  }
}

/**
 * Which screen the window takes. "full" (the default): the whole terminal
 * -- its alternate screen, the conversation scrolling inside the window,
 * the composer fixed at the bottom, the shell's lines back as they were
 * when it closes. "inline": under the shell's last line, the conversation
 * in the terminal's own scrollback (`delfin-agent chat`'s way).
 *
 * DELFIN_TUI_SCREEN (which `delfin-agent tui --inline` sets) before the
 * setting; with neither, inline under GNU screen, whose alternate screen
 * is off unless .screenrc turns it on ("altscreen on"), else full.
 */
export function screenModeOf(setting: Config["screen"],
                             env: Record<string, string | undefined> = process.env): "full" | "inline" {
  const asked = String(env.DELFIN_TUI_SCREEN ?? "").trim().toLowerCase()
  if (asked === "full" || asked === "inline") return asked
  if (setting === "full" || setting === "inline") return setting
  return env.STY ? "inline" : "full"
}
