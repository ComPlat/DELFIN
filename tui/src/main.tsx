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
import { createCliRenderer } from "@opentui/core"
import { render } from "@opentui/solid"
import { createSignal } from "solid-js"
import { App } from "./app"
import { History } from "./composer/history"
import { loadConfig } from "./config/config"
import { Keymap } from "./keymap/keymap"
import type { AnyEvent } from "./protocol/generated"
import { startSessionView } from "./session"
import { ensureSession, type Asked } from "./store/state"
import { createAppStore, withReplay } from "./store/store"
import { themeFor } from "./theme/theme"
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

const loaded = await loadConfig()
const keymap = new Keymap(loaded.config.keys)
const store = createAppStore()
const connection = new ServiceConnection(socket)
connection.rpc.onEvent((event) => store.push(event))
connection.onState((state) => store.update((d) => { d.connection = state }))

let hello: ServiceHello
try {
  hello = await connection.connect()
} catch (error) {
  console.error(`The agent service does not answer on ${socket}: ${(error as Error).message}`)
  process.exit(3)
}

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
    await withReplay(store,
      () => connection.rpc.request("session.attach", { sessionId: target, turns: loaded.config.replayTurns }),
      (attached) => {
        store.update((d) => { ensureSession(d, target, attached.session) })
        const title = String(attached.session?.title || "") || target.slice(0, 8)
        view.show(target, { divider: `${running ? "now" : "resumed"}: ${title}` })
        current = target
        return attached.events as unknown as AnyEvent[]
      })
    if (previous && previous !== target) {
      await connection.rpc.request("session.detach", { sessionId: previous }).catch(() => undefined)
    }
  } catch (error) {
    view.note(`could not open that session: ${(error as Error).message}`)
  }
}

let current = ""
try {
  current = await openFirst()
} catch (error) {
  // Nothing to show: say why on a terminal that is the user's again.
  connection.close()
  renderer.destroy()
  await Bun.sleep(50)
  process.stderr.write(`The session could not start: ${(error as Error).message}\n`)
  process.exit(3)
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

let quitting = false
function quit(): void {
  if (quitting) return
  quitting = true
  void (async () => {
    if (current) {
      await connection.rpc.request("session.detach", { sessionId: current }).catch(() => undefined)
    }
    store.flush()
    await view.queue.drain()
    connection.close()
    view.stop()
    // OpenTUI flushes queued scrollback after destroy; the process ends on
    // its own once that is written (never exit from onDestroy).
    renderer.destroy()
  })()
}

if (current) {
  await render(() => (
    <App renderer={renderer} store={store} keymap={keymap} theme={theme} history={History.load()} sentAs={view.sentAs}
         client={connection.rpc} onQuit={quit} onRedraw={() => void view.redraw()}
         onExpand={() => view.expandLast()} onResume={(id, ws) => void switchTo(id, ws)}
         onThinking={(on) => view.setShowThinking(on)} thinking={() => view.showThinking}
         onLines={(title, lines) => view.lines(title, lines)} />
  ), renderer)
}
