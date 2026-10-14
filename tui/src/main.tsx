/**
 * `delfin-agent tui`: the window. Started by the Python launcher, which has
 * made sure the agent service runs:
 *
 *   DELFIN_TUI_SOCKET  the service's socket
 *   DELFIN_TUI_START   the first session's settings, as JSON (chat's flags)
 *
 * It connects, starts the first session (or attaches to `resume`), and
 * draws it inline: finished output goes into the terminal's own scrollback,
 * only the footer is live. Leaving the window leaves the sessions running.
 */
import { appendFileSync } from "node:fs"
import { dirname } from "node:path"
import { createCliRenderer, getTreeSitterClient } from "@opentui/core"
import { render } from "@opentui/solid"
import { createSignal } from "solid-js"
import { Pins } from "./agents/pins"
import { createAgentScreen } from "./agents/screen"
import { App } from "./app"
import { History } from "./composer/history"
import { loadConfig, screenModeOf } from "./config/config"
import { Keymap } from "./keymap/keymap"
import type { AnyEvent } from "./protocol/generated"
import { startSessionView } from "./session"
import { ensureSession, type Asked } from "./store/state"
import { createAppStore, eventsOf, eventsOfUnknown, listed, withReplay } from "./store/store"
import { readState, updateState } from "./config/state"
import { noColor, themeFor } from "./theme/theme"
import { bounded } from "./scrollback/queue"
import { createScreenItems } from "./scrollback/screenview"
import { ServiceConnection, type ServiceHello } from "./transport/rpc"

// The kitty-graphics probe retitles tmux panes; this TUI draws no images.
process.env.OPENTUI_GRAPHICS ??= "0"

const socket = process.env.DELFIN_TUI_SOCKET ?? ""
if (!socket) {
  console.error("Start the window with `delfin-agent tui`: it needs the agent service.")
  process.exit(2)
}
let start: Record<string, unknown> = {}
try {
  start = JSON.parse(process.env.DELFIN_TUI_START ?? "{}")
} catch {
  console.error("DELFIN_TUI_START is not JSON; starting with defaults")
}

// DELFIN_TUI_DEBUG=<file>: when each startup step was reached (a window
// that is slow to come up says where it waits).
const debugFile = process.env.DELFIN_TUI_DEBUG ?? ""
const t0 = Date.now()
function step(what: string): void {
  if (!debugFile) return
  try {
    appendFileSync(debugFile, `${((Date.now() - t0) / 1000).toFixed(3)}s ${what}\n`)
  } catch {
    // a debug log that cannot be written is no reason to stop
  }
}
step("start")

const loaded = await loadConfig()
const keymap = new Keymap(loaded.config.keys)
const store = createAppStore()

/** What starts the next service when this one vanished (the launcher's
 *  command, DELFIN_TUI_SERVICE); its error, when it could not, is said. */
const serviceCommand: unknown = (() => {
  try {
    return JSON.parse(process.env.DELFIN_TUI_SERVICE ?? "null")
  } catch {
    return null
  }
})()
const startService = Array.isArray(serviceCommand) && serviceCommand.length && serviceCommand.every((a) => typeof a === "string")
  ? async (): Promise<boolean> => {
      // In the runtime directory (the socket's, ours alone): `python -m`
      // imports from the working directory first, and this window's is
      // the project folder -- a json.py there would run outside the cage.
      const proc = Bun.spawn(serviceCommand as string[], { stdin: "ignore", stdout: "ignore", stderr: "pipe",
                                                           cwd: dirname(socket) })
      const timer = setTimeout(() => proc.kill(), 150_000)
      const code = await proc.exited
      clearTimeout(timer)
      if (code !== 0) {
        const said = (await new Response(proc.stderr).text()).trim().split("\n").slice(-3).join(" ")
        view.note(`the agent service could not be started again${said ? `: ${said}` : ""}`)
      }
      return code === 0
    }
  : undefined
const connection = new ServiceConnection(socket, undefined, { startService })
connection.rpc.onEvent((event) => store.push(event))
connection.onState((state) => store.update((d) => { d.connection = state }))

step("config read")
let hello: ServiceHello
try {
  hello = await connection.connect()
} catch (error) {
  console.error(`The agent service does not answer on ${socket}: ${(error as Error).message}`)
  process.exit(3)
}

step("hello answered")
// The whole terminal (the default) or inline under the shell's output
// (`--inline`, tui.json "screen", GNU screen): config.ts screenModeOf.
const full = screenModeOf(loaded.config.screen) === "full"
const renderer = await createCliRenderer({
  screenMode: full ? "alternate-screen" : "split-footer",
  footerHeight: 4,
  // Captured stdout is a split-footer thing (it goes to the scrollback);
  // the window itself never writes to stdout while it runs.
  externalOutputMode: full ? "passthrough" : "capture-stdout",
  // Full screen takes the mouse: the wheel scrolls the conversation, a
  // drag selects and copies (shift+drag is the terminal's own selection).
  // Buttons and drags only (no reports of every move: nothing hovers), and
  // a click never takes the keys from the prompt (its focus is the app's).
  useMouse: full,
  enableMouseMovement: false,
  autoFocus: false,
  exitOnCtrlC: false,
  targetFps: 30,
  maxFps: 60,
  consoleMode: "disabled",
  useKittyKeyboard: {},
})
step("renderer created")
// The highlighter starts its worker and loads the markdown grammar when it
// is first used -- some hundreds of ms in which the first finished block
// waits to be shown: started now, while the session is attached.
void getTreeSitterClient().highlightOnce("# DELFIN\n\n**ready**", "markdown").catch(() => undefined)
// NO_COLOR asks for no colours of the window's own, unless tui.json names
// a theme. "auto" draws at once with the theme the terminal said last time
// and asks it again meanwhile: waiting for the answer before the first
// frame cost 300 ms of every start where none comes (tmux).
const auto = loaded.config.theme === "auto" && !noColor()
const remembered = readState().themeMode
let mode: "dark" | "light" | "mono" = loaded.config.theme !== "auto" ? loaded.config.theme
  : noColor() ? "mono"
  : remembered === "light" ? "light" : "dark"
const [theme, setTheme] = createSignal(themeFor(mode, loaded.config.colors))
if (auto) {
  void renderer.waitForThemeMode(1500).then((said) => {
    if (said !== "dark" && said !== "light") return
    if (said !== remembered) updateState({ themeMode: said })
    if (said === mode) return
    mode = said
    setTheme(themeFor(said, loaded.config.colors))
  }).catch(() => undefined)
}
// `--verbose` shows the thinking, as in chat; tui.json can make it the default.
const showThinking = loaded.config.showThinking || Boolean(start.verbose)
const screen = full ? createScreenItems() : undefined
const view = startSessionView(renderer, store, theme, { showThinking, replayTurns: loaded.config.replayTurns, screen })
store.update((d) => {
  d.mayAct = Boolean(hello.mayAct)
  d.host = String(hello.host ?? "")
})
// Every session the service runs, for the agent view (each summary that
// follows keeps its row up to date).
listed(store, (hello.sessions ?? []) as { sessionId?: string }[])
for (const problem of [...loaded.problems, ...keymap.problems]) view.note(problem)
// Where else this user's sessions are (other login nodes), from the launcher.
try {
  for (const said of JSON.parse(process.env.DELFIN_TUI_NOTES ?? "[]") as unknown[]) {
    if (typeof said === "string" && said) view.note(said)
  }
} catch {
  // none to say
}

/**
 * Which session the window shows is the latest one asked for: opening the
 * first, switching (the resume picker, the agent view) and starting a new
 * one each take a number, and one that comes back after a later request
 * is not put in view -- it runs on in the service, in the agent view.
 */
let wanted = 0

/** The session the newest open asks for, and the id it runs under once
 *  known (a resumed one may run under another). An open it superseded
 *  lets its session go -- unless it is this one: the newer open is about
 *  to show it, and a detach now would leave it in view unfollowed. */
let wantedId = ""
let wantedTarget = ""

/** Whether a superseded open may let ``id`` go: not the session in view,
 *  nor the one the newest open asks for. */
function letGo(id: string): boolean {
  return Boolean(id) && id !== current && id !== wantedId && id !== wantedTarget
}

/** The longest an attach may take -- a wake starts a worker, its MCP
 *  servers and a resume -- before its hold lets the session's events go. */
const ATTACH_TIMEOUT_MS = 240_000

/** What a resume needs of the window's settings: a --fake window resumes
 *  fake sessions (never a real model in one started for tests). */
function resumeSettings(): Record<string, unknown> {
  const out: Record<string, unknown> = {}
  if (start.fake) out.fake = true
  if (start.fakeDelay !== undefined) out.fakeDelay = start.fakeDelay
  return out
}

/** Open the first session: a new one, or the one to resume if it runs here. */
async function openFirst(): Promise<string> {
  const mine = ++wanted
  const running = (hello.sessions ?? []) as { sessionId: string }[]
  const resume = String(start.resume ?? "")
  wantedId = resume
  wantedTarget = ""
  if (resume && running.some((s) => s.sessionId === resume)) {
    let shown = false
    await withReplay(store,
      () => connection.rpc.request("session.attach", { sessionId: resume, turns: loaded.config.replayTurns },
                                   { timeoutMs: ATTACH_TIMEOUT_MS }),
      (attached) => {
        store.update((d) => {
          const s = ensureSession(d, resume, attached.session)
          s.state = attached.state
          s.pending = attached.pending as unknown as Asked[]
          s.ended = null
        })
        if (mine !== wanted) return []
        view.show(resume)
        current = resume
        shown = true
        return attached.events as unknown as AnyEvent[]
      }, eventsOf(resume))
    if (!shown && letGo(resume)) await leave(resume)
    return resume
  }
  const { verbose: _verbose, ...settings } = start
  const created = await withReplay(store,
    () => connection.rpc.request("session.create", { ...settings, cwd: String(start.cwd ?? process.cwd()) }),
    (answer) => {
      const id = answer.sessionId as string
      store.update((d) => {
        const s = ensureSession(d, id, answer.session)
        s.state = answer.state
        s.pending = answer.pending as unknown as Asked[]
      })
      if (mine !== wanted) {
        if (letGo(id)) void leave(id)
        return []
      }
      view.show(id)
      current = id
      // The banner (the session as the service saw it start), then what
      // the session has said since its worker started.
      return [{ type: "session.info", session: id, ts: Date.now() / 1000, seq: 0, ...answer.session } as never,
              ...((answer.events ?? []) as unknown as AnyEvent[])]
    }, eventsOfUnknown(store))
  return created.sessionId as string
}

/**
 * Show another session in this window (the resume picker, the agent
 * view): a divider, its last turns, then its live events. One that runs
 * in the service (or sleeps there) is attached; a saved one is started
 * again from the session store first. The session left keeps running.
 */
async function switchTo(id: string, workspace: string): Promise<void> {
  const mine = ++wanted
  wantedId = id
  wantedTarget = ""
  let target = id
  try {
    const listed = await connection.rpc.request("sessions.list", {})
    const running = (listed.sessions as { sessionId: string }[]).some((s) => s.sessionId === id)
    if (!running) {
      // The id the resumed session runs under is the one to attach to.
      const created = await connection.rpc.request("session.create",
                                                   { cwd: workspace || process.cwd(), resume: id, attach: false,
                                                     ...resumeSettings() })
      target = String(created.sessionId)
    }
    if (mine !== wanted) return
    wantedTarget = target
    let shown = false
    // Shown before in this window: its replay is new rows again.
    store.unsee(target)
    await withReplay(store,
      () => connection.rpc.request("session.attach", { sessionId: target, turns: loaded.config.replayTurns },
                                   { timeoutMs: ATTACH_TIMEOUT_MS }),
      (attached) => {
        store.update((d) => {
          const s = ensureSession(d, target, attached.session)
          s.state = attached.state
          s.pending = attached.pending as unknown as Asked[]
          s.ended = null
        })
        if (mine !== wanted) return []
        const previous = current
        const title = String(attached.session?.title || "") || target.slice(0, 8)
        view.show(target, { divider: `${running ? "now" : "resumed"}: ${title}` })
        current = target
        shown = true
        if (previous && previous !== target) void leave(previous)
        return attached.events as unknown as AnyEvent[]
      }, eventsOf(target))
    if (!shown && letGo(target)) await leave(target)
  } catch (error) {
    if (mine === wanted) {
      // The newest open failed: nothing it asked for stays followed out of
      // view, an older open of the same session included.
      wantedId = ""
      wantedTarget = ""
      if (letGo(target)) void leave(target)
    }
    view.note(`could not open that session: ${(error as Error).message}`)
  }
}

/**
 * Stop following a session that left the view: it runs on in the service,
 * and its row lives on its summaries. What the window knew of its turn is
 * let go -- no turn.ended would ever end it here.
 */
async function leave(id: string): Promise<void> {
  await connection.rpc.request("session.detach", { sessionId: id }).catch(() => undefined)
  store.update((d) => {
    const s = d.sessions[id]
    if (!s || s.state === "ended") return
    s.turn = null
    s.tail = ""
    s.thinking = ""
    s.wait = ""
    s.pending = []
  })
}

/**
 * A new session (the agent view's typing + enter), in this window's
 * folder with the first session's settings, shown after a divider -- if
 * nothing else was asked for meanwhile. The session in view until then
 * stays in view and keeps running in the service. Its id, or null when it
 * did not start (the note says why).
 */
async function newSession(): Promise<string | null> {
  const mine = ++wanted
  wantedId = ""
  wantedTarget = ""
  const { verbose: _verbose, resume: _resume, continue: _continue, fork: _fork, name: _name, ...settings } = start
  let id = ""
  let shown = false
  try {
    await withReplay(store,
      () => connection.rpc.request("session.create", { ...settings, cwd: String(start.cwd ?? process.cwd()) }),
      (answer) => {
        id = answer.sessionId as string
        store.update((d) => {
          const s = ensureSession(d, id, answer.session)
          s.state = answer.state
          s.pending = answer.pending as unknown as Asked[]
        })
        if (mine !== wanted) return []
        const previous = current
        view.show(id, { divider: "a new session" })
        current = id
        shown = true
        if (previous && previous !== id) void leave(previous)
        return [{ type: "session.info", session: id, ts: Date.now() / 1000, seq: 0, ...answer.session } as never,
                ...((answer.events ?? []) as unknown as AnyEvent[])]
      }, eventsOfUnknown(store))
    // Asked for something else meanwhile: it runs on, in the agent view.
    if (!shown && letGo(id)) await leave(id)
    return id
  } catch (error) {
    view.note(`the new session could not start: ${(error as Error).message}`)
    return null
  }
}

let current = ""
let quitting = false

/** Leave the window; the sessions keep running in the agent service. A
 *  second ctrl+d while leaving (a service that does not answer) leaves at
 *  once. */
function quit(): void {
  if (quitting) {
    leaveNow()
    return
  }
  quitting = true
  // The scrollback is held while the agent view has the alternate screen.
  agentScreen?.leave()
  void (async () => {
    if (current) {
      await connection.rpc.request("session.detach", { sessionId: current }, { timeoutMs: 1000 }).catch(() => undefined)
    }
    store.flush()
    await bounded(view.queue.drain(), 3000)
    connection.close()
    view.stop()
    // OpenTUI flushes queued scrollback after destroy; the process ends on
    // its own once that is written (never exit from onDestroy).
    renderer.destroy()
    // Full screen gives the shell its lines back as they were: one line
    // says what became of the session -- written as the process ends, when
    // the terminal is surely the shell's again (OpenTUI restores it after
    // destroy() returns).
    if (full) {
      const note = leavingNote()
      process.once("exit", () => {
        try {
          process.stdout.write(note)
        } catch {
          // a terminal already gone
        }
      })
    }
  })()
}

/** What leaving a full-screen window says: why it closed when the service
 *  went (said in the window too, but that screen is gone), else whether
 *  the session runs on -- and how to come back to it. */
function leavingNote(): string {
  const plain = (text: string) => text.replace(/[\u0000-\u001f\u007f-\u009f]/g, " ")
  if (connection.problem) return `DELFIN: ${plain(connection.problem)}\n`
  const s = current ? store.state.sessions[current] : undefined
  if (!s) return ""
  const name = s.title || s.summary?.title || `session ${s.id.slice(0, 8)}`
  if (s.state === "ended") return `DELFIN: ${plain(name)} has ended${s.ended?.reason ? ` (${plain(s.ended.reason)})` : ""}.\n`
  return `DELFIN: ${plain(name)} keeps running in the agent service -- `
    + `delfin-agent tui --resume ${s.id} opens it again (or ← in any window).\n`
}

function leaveNow(): void {
  connection.close()
  try {
    renderer.destroy()
  } catch {
    // already gone
  }
  setTimeout(() => process.exit(0), 100)
}

// A signal -- a hangup when the SSH connection drops, a SIGTERM -- ends the
// window. OpenTUI destroys the renderer on it but leaves the process to
// the app: with the service's socket open it would never end (and its
// session, still shown "in a window", would never sleep).
let signalled = 0
for (const [name, number] of [["SIGHUP", 1], ["SIGINT", 2], ["SIGQUIT", 3], ["SIGTERM", 15]] as const) {
  // Before OpenTUI's own listener, which destroys the renderer at once.
  process.prependListener(name, () => { signalled ||= number })
}
renderer.on("destroy", () => {
  if (quitting) return
  quitting = true
  // The shell's way: a window ended by a signal exits with 128 + its number.
  if (signalled) process.exitCode = 128 + signalled
  connection.close()
  view.stop()
  ;(setTimeout(() => process.exit(signalled ? 128 + signalled : 0), 1500) as { unref?: () => void }).unref?.()
})

// Made once: a JSX prop is read through a getter, each read running its
// expression again.
const pins = Pins.load()
// The agent view on the alternate screen -- except under GNU screen, whose
// alternate screen is off unless .screenrc turns it on ("altscreen"): there
// the view grows the footer instead.
// Full screen has the view in the window already: no screen to switch to.
const agentScreen = full || (process.env.DELFIN_TUI_AGENTS_SCREEN ?? (process.env.STY ? "footer" : "alternate")) === "footer"
  ? undefined : createAgentScreen(renderer, view.queue)

// The window takes keys from the start: a session that is slow to start
// (MCP servers, a long resume) can still be left with ctrl+d or ctrl+c.
await render(() => (
  <App renderer={renderer} store={store} keymap={keymap} theme={theme} history={History.load()} sentAs={view.sentAs}
       client={connection.rpc} onQuit={quit} onRedraw={() => void view.redraw()}
       onExpand={() => view.expandLast()} details={view.details} canShowDetails={() => view.canShowDetails}
       onDetails={() => view.toggleDetails()}
       onResume={(id, ws) => void switchTo(id, ws)}
       onNewSession={newSession} pins={pins} agentScreen={agentScreen} screen={screen}
       onThinking={(on) => view.setShowThinking(on)} thinking={() => view.showThinking}
       onLines={(title, lines) => view.lines(title, lines)} />
), renderer)

step("window rendered")
store.update((d) => { d.starting = true })
try {
  await openFirst()
} catch (error) {
  if (quitting) {
    // Left while it started: nothing to say (the session, if it started,
    // runs on in the agent service).
    leaveNow()
  } else {
    // Nothing to show: say why on a terminal that is the user's again.
    connection.close()
    renderer.destroy()
    await Bun.sleep(50)
    process.stderr.write(`The session could not start: ${(error as Error).message}\n`)
    process.exit(3)
  }
} finally {
  store.update((d) => { d.starting = false })
  step(current ? "session open" : "no session")
}

// After a reconnect: the session is attached again (a new service wakes it
// from where the old one left it); what this window has already shown is
// not shown twice (the store drops seen sequence numbers), and what it
// knew of a turn that ended meanwhile is let go.
connection.onUp(async (again) => {
  // What this (maybe new) service runs is what the window knows now.
  listed(store, (again.sessions ?? []) as { sessionId?: string }[])
  if (!current) return
  const id = current
  try {
    await withReplay(store,
      () => connection.rpc.request("session.attach", { sessionId: id, turns: 0 }, { timeoutMs: ATTACH_TIMEOUT_MS }),
      (attached) => {
        store.update((d) => {
          const s = ensureSession(d, id, attached.session)
          s.state = attached.state
          s.pending = attached.pending as unknown as Asked[]
          s.ended = null
          if (attached.state === "idle") {
            s.turn = null
            s.tail = ""
            s.thinking = ""
            s.wait = ""
          }
        })
        return attached.events as unknown as AnyEvent[]
      }, eventsOf(id))
  } catch (error) {
    view.note(`the session is not in the agent service any more: ${(error as Error).message}`)
  }
})
let saidClosed = false
connection.onState((state) => {
  if (state !== "closed" || !connection.problem || saidClosed) return
  saidClosed = true
  view.note(connection.problem)
  // Stopped on purpose: there is nothing to come back to here.
  if (connection.stoppedOnPurpose) quit()
})
