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
import { briefOf, openingOf, shortPath, suggestedPrompt } from "./agents/attention"
import { Pins } from "./agents/pins"
import { createAgentScreen } from "./agents/screen"
import { App } from "./app"
import { History } from "./composer/history"
import { loadConfig, screenModeOf } from "./config/config"
import { Keymap } from "./keymap/keymap"
import type { AnyEvent, AttentionItem } from "./protocol/generated"
import { leavingNote } from "./leaving"
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

/**
 * What needs the user, from every DELFIN process (a service of version 3
 * and later; one before it has nothing to say): the list now, then its
 * changes as attention.updated. Into the store's queue like those, so a
 * change that came before the answer cannot overwrite it afterwards.
 */
async function loadAttention(service: ServiceHello): Promise<void> {
  const event = (items: unknown, version?: number) =>
    ({ type: "attention.updated", session: "", ts: Date.now() / 1000, items, version }) as never
  if (Number(service.serviceVersion ?? 0) < 3) {
    store.update((d) => {
      d.attention = []
      d.attentionVersion = -1
    })
    return
  }
  try {
    // Versioned (from the service's start on): an update that came before
    // this answer, though read after it, is not overwritten by it.
    const got = await connection.rpc.request("attention.list", {})
    store.push(event(got?.items ?? [], typeof got?.version === "number" ? got.version : undefined))
  } catch {
    // The next attention.updated brings it.
  }
}
void loadAttention(hello)
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

interface SwitchOptions {
  /** A saved session is opened as a copy, under a new id: another
   *  process has it open (two writing one conversation lose saves). */
  fork?: boolean
  /** The divider's words, given the session's title. */
  divider?: (title: string, resumed: boolean) => string
  /** What went wrong is the caller's to say. */
  quiet?: boolean
}

/** Shown; or not: superseded by a newer open (no error), or why not. */
type Switched = { shown: true; id: string } | { shown: false; error?: string }

/**
 * Show another session in this window (the resume picker, the agent
 * view): a divider, its last turns, then its live events. One that runs
 * in the service (or sleeps there) is attached; a saved one is started
 * again from the session store first (or a copy of it). The session left
 * keeps running.
 */
async function switchTo(id: string, workspace: string, options: SwitchOptions = {}): Promise<Switched> {
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
                                                     ...(options.fork ? { fork: true } : {}), ...resumeSettings() })
      target = String(created.sessionId)
    }
    if (mine !== wanted) return { shown: false }
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
        view.show(target, { divider: options.divider?.(title, !running) ?? `${running ? "now" : "resumed"}: ${title}` })
        current = target
        shown = true
        if (previous && previous !== target) void leave(previous)
        return attached.events as unknown as AnyEvent[]
      }, eventsOf(target))
    if (!shown && letGo(target)) await leave(target)
    return shown ? { shown: true, id: target } : { shown: false }
  } catch (error) {
    if (mine === wanted) {
      // The newest open failed: nothing it asked for stays followed out of
      // view, an older open of the same session included.
      wantedId = ""
      wantedTarget = ""
      if (letGo(target)) void leave(target)
    }
    const said = (error as Error).message
    if (!options.quiet) view.note(`could not open that session: ${said}`)
    return mine === wanted ? { shown: false, error: said || "it did not start" } : { shown: false }
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

interface NewOptions {
  /** Its folder (this window's by default). */
  cwd?: string
  /** What it is to look into (an attention item): its first row, and in
   *  front of its first prompt. */
  brief?: string
  briefTitle?: string
  /** What it is called (the agent view, the prompt box's chip). */
  title?: string
  divider?: string
  /** What went wrong is the caller's to say. */
  quiet?: boolean
}

/** A new session asked for: its id (null when it did not start, and
 *  why); superseded when something else was opened meanwhile -- then it
 *  runs on in the agent view, out of view. */
interface Started {
  id: string | null
  superseded: boolean
  error: string
}

/**
 * A new session (the agent view's typing + enter), in this window's
 * folder with the first session's settings, shown after a divider -- if
 * nothing else was asked for meanwhile. The session in view until then
 * stays in view and keeps running in the service. Its id, or null when it
 * did not start (the note says why).
 */
async function newSession(options: NewOptions = {}): Promise<string | null> {
  return (await startSession(options)).id
}

async function startSession(options: NewOptions = {}): Promise<Started> {
  const mine = ++wanted
  wantedId = ""
  wantedTarget = ""
  const { verbose: _verbose, resume: _resume, continue: _continue, fork: _fork, name: _name, ...settings } = start
  const brief = options.brief ? { brief: options.brief, briefTitle: options.briefTitle ?? "" } : {}
  const named = options.title ? { title: options.title } : {}
  let id = ""
  let shown = false
  try {
    await withReplay(store,
      () => connection.rpc.request("session.create", { ...settings, ...brief, ...named,
                                                       cwd: options.cwd || String(start.cwd ?? process.cwd()) }),
      (answer) => {
        id = answer.sessionId as string
        store.update((d) => {
          const s = ensureSession(d, id, answer.session)
          s.state = answer.state
          s.pending = answer.pending as unknown as Asked[]
        })
        if (mine !== wanted) return []
        const previous = current
        view.show(id, { divider: options.divider ?? "a new session" })
        current = id
        shown = true
        if (previous && previous !== id) void leave(previous)
        return [{ type: "session.info", session: id, ts: Date.now() / 1000, seq: 0, ...answer.session } as never,
                ...((answer.events ?? []) as unknown as AnyEvent[])]
      }, eventsOfUnknown(store))
    // Asked for something else meanwhile: it runs on, in the agent view.
    if (!shown && letGo(id)) await leave(id)
    return { id, superseded: !shown, error: "" }
  } catch (error) {
    const said = (error as Error).message
    if (!options.quiet) view.note(`the new session could not start: ${said}`)
    return { id: null, superseded: mine !== wanted, error: said || "it did not start" }
  }
}

/** What a session asked that is answered there now (a parked question,
 *  a confirmation, a plan): said under its replay, where the answer goes. */
function sayWhatItAsked(item: AttentionItem): void {
  const lines = item.kind === "question"
    ? [item.detail || item.text.replace(/^question:\s*/, ""), ...item.options.map((o, i) => `${i + 1}. ${o}`)]
    : (item.detail || item.text).split("\n")
  const what = item.kind === "question" ? "it asked, and nobody was there to answer"
    : item.kind === "plan" ? "it put up a plan, and nobody was there to approve it"
    : "it asked for a confirmation, and nobody was there to give it"
  // After the replay: the store's queue is presented first.
  store.flush()
  view.lines(`${what} -- answer it here`, lines.filter((l) => l.trim()))
}

/**
 * → on an item of what needs the user (attention.ts): its session, resumed
 * in this window -- or, while another process has it open, a copy; a
 * failed job's diagnosis session; else (or when that cannot be opened) a
 * new session about it, in its folder, with the item as its brief. The
 * prompt a new session suggests (for the prompt box: nothing is sent
 * until the user presses enter), or null.
 */
async function openItem(item: AttentionItem): Promise<string | null> {
  const opening = openingOf(item)
  const folder = item.workspace || String(start.cwd ?? process.cwd())
  const first = (text: string) => text.split("\n")[0]
  // What the window says about how it opened it: once that is in view --
  // said before, full screen would clear it with the screen.
  const said: string[] = []
  const sayNow = () => {
    store.flush()
    for (const text of said.splice(0)) view.note(text)
  }
  if (opening.how !== "new") {
    const where = item.host && item.host !== store.state.host ? ` on ${item.host}` : " elsewhere"
    // Named as the list named it: a copy's own title says "(fork)".
    const copy = () => `a copy of ${opening.title} (still open${where})`
    const plain = (title: string) => opening.how === "diagnosis" ? `${opening.title}: ${title}` : `resumed: ${title}`
    let copying = opening.how === "session" && opening.copy
    let got = await switchTo(opening.id, folder, { fork: copying, quiet: true, divider: copying ? copy : plain })
    if (!got.shown && got.error && /being written by another process/i.test(got.error)) {
      // Open elsewhere after all (no presence said so yet): a copy.
      copying = true
      got = await switchTo(opening.id, folder, { fork: true, quiet: true, divider: copy })
    }
    if (got.shown) {
      if (copying) said.push(`this is a copy of ${opening.title}, which is still open${where}: what is done here does not reach it`)
      sayNow()
      if (item.blocking) {
        // Answered there, as an ordinary message: the inbox's entry is
        // closed, so no answer to it goes out a second time.
        void connection.rpc.request("attention.dismiss", { id: item.id }).catch(() => undefined)
        sayWhatItAsked(item)
      }
      return null
    }
    if (got.error === undefined) return null          // something newer was asked for
    said.push(`${opening.title} could not be opened (${first(got.error)}): a new session about it instead`)
  }
  // Named after it: found again in the agent view by what it is about.
  const about = { brief: briefOf(item), briefTitle: item.name, title: item.name, divider: `about: ${item.name}`,
                  quiet: true }
  let made = await startSession({ ...about, cwd: folder })
  if (made.superseded) return null                    // something newer was asked for
  const here = String(start.cwd ?? process.cwd())
  if (!made.id && folder !== here) {
    // Its folder takes no session (gone, on another machine, refused):
    // this window's folder then, with the item as its brief all the same.
    said.push(`${shortPath(folder)} could not take a session (${first(made.error)}): started in ${shortPath(here)}`)
    made = await startSession({ ...about, cwd: here })
    if (made.superseded) return null
  }
  if (!made.id) {
    sayNow()
    view.note(`the session about ${item.name} could not start: ${first(made.error)}`)
    return null
  }
  sayNow()
  // Its prompt only for it: not when something else was opened meanwhile.
  return current === made.id ? suggestedPrompt(item) : null
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
    printWhenGone()
  })()
}

/** Once: what became of the session and the command that opens it again,
 *  written as the process ends -- when the terminal is surely the shell's
 *  again (OpenTUI restores it after destroy() returns), in both screens:
 *  full screen gives the shell its lines back, inline leaves the
 *  conversation above, and the line lands under either. */
let printed = false
function printWhenGone(): void {
  if (printed) return
  printed = true
  const s = current ? store.state.sessions[current] : undefined
  const note = leavingNote({
    problem: connection.problem || undefined,
    session: s ? { id: s.id, title: s.title || s.summary?.title || "", state: s.state, ended: s.ended } : undefined,
  })
  if (!note) return
  process.once("exit", () => {
    try {
      process.stdout.write(note)
    } catch {
      // a terminal already gone
    }
  })
}

function leaveNow(): void {
  connection.close()
  printWhenGone()
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
  // A hangup has no terminal left to print to; any other end says the way back.
  if (signalled !== 1) printWhenGone()
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
       onNewSession={() => newSession()} onOpenItem={openItem} pins={pins} agentScreen={agentScreen} screen={screen}
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
  void loadAttention(again)
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
