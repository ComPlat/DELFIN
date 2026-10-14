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
import { createCliRenderer } from "@opentui/core"
import { render } from "@opentui/solid"
import { createSignal } from "solid-js"
import { Pins } from "./agents/pins"
import { createAgentScreen } from "./agents/screen"
import { App } from "./app"
import { History } from "./composer/history"
import { loadConfig } from "./config/config"
import { Keymap } from "./keymap/keymap"
import type { AnyEvent } from "./protocol/generated"
import { startSessionView } from "./session"
import { ensureSession, type Asked } from "./store/state"
import { createAppStore, withReplay } from "./store/store"
import { themeFor } from "./theme/theme"
import { bounded } from "./scrollback/queue"
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
const connection = new ServiceConnection(socket)
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
const renderer = await createCliRenderer({
  screenMode: "split-footer",
  footerHeight: 4,
  externalOutputMode: "capture-stdout",
  useMouse: false,
  exitOnCtrlC: false,
  targetFps: 30,
  maxFps: 60,
  consoleMode: "disabled",
  useKittyKeyboard: {},
})
step("renderer created")
const mode = loaded.config.theme === "auto"
  ? await renderer.waitForThemeMode(300).catch(() => null)
  : loaded.config.theme
const [theme] = createSignal(themeFor(mode, loaded.config.colors))
// `--verbose` shows the thinking, as in chat; tui.json can make it the default.
const showThinking = loaded.config.showThinking || Boolean(start.verbose)
const view = startSessionView(renderer, store, theme, { showThinking, replayTurns: loaded.config.replayTurns })
store.update((d) => {
  d.mayAct = Boolean(hello.mayAct)
  d.host = String(hello.host ?? "")
})
// Every session the service runs, for the agent view (each summary that
// follows keeps its row up to date).
for (const summary of (hello.sessions ?? []) as { sessionId?: string }[]) {
  if (summary?.sessionId) {
    store.push({ type: "session.summary", session: summary.sessionId, ts: Date.now() / 1000, summary } as never)
  }
}
for (const problem of [...loaded.problems, ...keymap.problems]) view.note(problem)

/** Open the first session: a new one, or the one to resume if it runs here. */
async function openFirst(): Promise<string> {
  const running = (hello.sessions ?? []) as { sessionId: string }[]
  const resume = String(start.resume ?? "")
  if (resume && running.some((s) => s.sessionId === resume)) {
    await withReplay(store,
      () => connection.rpc.request("session.attach", { sessionId: resume, turns: loaded.config.replayTurns }),
      (attached) => {
        store.update((d) => { ensureSession(d, resume, attached.session) })
        view.show(resume)
        return attached.events as unknown as AnyEvent[]
      })
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
      view.show(id)
      // The banner (the session as the service saw it start), then what
      // the session has said since its worker started.
      return [{ type: "session.info", session: id, ts: Date.now() / 1000, seq: 0, ...answer.session } as never,
              ...((answer.events ?? []) as unknown as AnyEvent[])]
    })
  return created.sessionId as string
}

/**
 * Show another session in this window (the resume picker): a divider, its
 * last turns, then its live events. One that runs in the service is
 * attached; a saved one is started again from the session store first.
 * The session left keeps running in the service.
 */
async function switchTo(id: string, workspace: string): Promise<void> {
  try {
    const listed = await connection.rpc.request("sessions.list", {})
    const running = (listed.sessions as { sessionId: string }[]).some((s) => s.sessionId === id)
    let target = id
    if (!running) {
      // The id the resumed session runs under is the one to attach to.
      const created = await connection.rpc.request("session.create",
                                                   { cwd: workspace || process.cwd(), resume: id, attach: false })
      target = String(created.sessionId)
    }
    const previous = current
    // Shown before in this window: its replay is new rows again.
    store.unsee(target)
    await withReplay(store,
      () => connection.rpc.request("session.attach", { sessionId: target, turns: loaded.config.replayTurns }),
      (attached) => {
        store.update((d) => {
          const s = ensureSession(d, target, attached.session)
          s.state = attached.state
          s.pending = attached.pending as unknown as Asked[]
        })
        const title = String(attached.session?.title || "") || target.slice(0, 8)
        view.show(target, { divider: `${running ? "now" : "resumed"}: ${title}` })
        current = target
        return attached.events as unknown as AnyEvent[]
      })
    if (previous && previous !== target) await leave(previous)
  } catch (error) {
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
 * folder with the first session's settings, shown after a divider. The
 * session in view until now keeps running in the service. Prompts typed
 * meanwhile wait for it (`starting`).
 */
async function newSession(): Promise<boolean> {
  const { verbose: _verbose, resume: _resume, continue: _continue, fork: _fork, name: _name, ...settings } = start
  const previous = current
  store.update((d) => { d.starting = true })
  try {
    await withReplay(store,
      () => connection.rpc.request("session.create", { ...settings, cwd: String(start.cwd ?? process.cwd()) }),
      (answer) => {
        const id = answer.sessionId as string
        store.update((d) => {
          const s = ensureSession(d, id, answer.session)
          s.state = answer.state
          s.pending = answer.pending as unknown as Asked[]
        })
        view.show(id, { divider: "a new session" })
        current = id
        return [{ type: "session.info", session: id, ts: Date.now() / 1000, seq: 0, ...answer.session } as never,
                ...((answer.events ?? []) as unknown as AnyEvent[])]
      })
    if (previous && previous !== current) await leave(previous)
    return true
  } catch (error) {
    view.note(`the new session could not start: ${(error as Error).message}`)
    return false
  } finally {
    store.update((d) => { d.starting = false })
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
  })()
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

// Made once: a JSX prop is read through a getter, each read running its
// expression again.
const pins = Pins.load()
// The agent view on the alternate screen -- except under GNU screen, whose
// alternate screen is off unless .screenrc turns it on ("altscreen"): there
// the view grows the footer instead.
const agentScreen = (process.env.DELFIN_TUI_AGENTS_SCREEN ?? (process.env.STY ? "footer" : "alternate")) === "footer"
  ? undefined : createAgentScreen(renderer, view.queue)

// The window takes keys from the start: a session that is slow to start
// (MCP servers, a long resume) can still be left with ctrl+d or ctrl+c.
await render(() => (
  <App renderer={renderer} store={store} keymap={keymap} theme={theme} history={History.load()} sentAs={view.sentAs}
       client={connection.rpc} onQuit={quit} onRedraw={() => void view.redraw()}
       onExpand={() => view.expandLast()} onResume={(id, ws) => void switchTo(id, ws)}
       onNewSession={newSession} pins={pins} agentScreen={agentScreen}
       onThinking={(on) => view.setShowThinking(on)} thinking={() => view.showThinking}
       onLines={(title, lines) => view.lines(title, lines)} />
), renderer)

step("window rendered")
store.update((d) => { d.starting = true })
try {
  current = await openFirst()
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

// After a reconnect: the session is attached again; what this window has
// already shown is not shown twice (the store drops seen sequence numbers).
connection.onUp(async () => {
  if (!current) return
  try {
    await withReplay(store,
      () => connection.rpc.request("session.attach", { sessionId: current, turns: 0 }),
      (attached) => attached.events as unknown as AnyEvent[])
  } catch (error) {
    view.note(`the session is not in the agent service any more: ${(error as Error).message}`)
  }
})
connection.onState((state) => {
  if (state === "closed" && connection.problem) view.note(connection.problem)
})
