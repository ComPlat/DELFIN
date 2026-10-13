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

const REPLAY_LIMIT = 400

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
const view = startSessionView(renderer, store, theme,
                              { showThinking: loaded.config.showThinking, replayLimit: REPLAY_LIMIT })
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
  const created = await withReplay(store,
    () => connection.rpc.request("session.create", { ...start, cwd: String(start.cwd ?? process.cwd()) }),
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
         client={connection.rpc} onQuit={quit} onRedraw={() => void view.redraw()} />
  ), renderer)
}
