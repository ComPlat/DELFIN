/**
 * The window's own small state across starts, `~/.delfin/tui_state.json`
 * (DELFIN_SCRATCH_STATE moves it, as it moves every ~/.delfin sink; so does
 * DELFIN_TUI_STATE): the pinned sessions, the terminal's last detected
 * theme. Each writer changes only its own keys -- read, merged, written
 * atomically -- so windows and features do not drop each other's. A file
 * that cannot be read is no state; one that cannot be written keeps the
 * state for this window.
 */
import { existsSync, mkdirSync, readFileSync, renameSync, writeFileSync } from "node:fs"
import { homedir } from "node:os"
import { dirname, join } from "node:path"

export function statePath(env: Record<string, string | undefined> = process.env): string {
  const home = env.DELFIN_SCRATCH_STATE || env.HOME || homedir()
  return env.DELFIN_TUI_STATE || join(home, ".delfin", "tui_state.json")
}

export function readState(path: string = statePath()): Record<string, unknown> {
  try {
    if (!existsSync(path)) return {}
    const data = JSON.parse(readFileSync(path, "utf8"))
    return data && typeof data === "object" && !Array.isArray(data) ? data : {}
  } catch {
    return {}
  }
}

/** Set `patch`'s keys, keeping every other key in the file. */
export function updateState(patch: Record<string, unknown>, path: string = statePath()): void {
  try {
    mkdirSync(dirname(path), { recursive: true, mode: 0o700 })
    const tmp = `${path}.${process.pid}.tmp`
    writeFileSync(tmp, JSON.stringify({ ...readState(path), ...patch }, null, 2) + "\n", { mode: 0o600 })
    renameSync(tmp, path)
  } catch {
    // a read-only home keeps the state for this window only
  }
}
