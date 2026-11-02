/**
 * The window: one session in view, its scrollback, and the footer.
 *
 * Keys go to the agent view while it is open (its arrows and ctrl keys;
 * letters are a new session's prompt). With a question open the prompt
 * keeps the keys: the question takes its own keys from an empty prompt
 * (approval.tsx), everything else is typing, and a line sent is the
 * question's typed answer -- or, for an approval, a message; the window's
 * other actions wait. Else the composer (completion, history, search),
 * then the keymap's actions. ctrl+c and ctrl+d work in all of them.
 * Leaving the window (ctrl+d, or ctrl+c twice) never ends a session: it
 * keeps running in the agent service, and the next window attaches to it
 * again.
 */
import type { CliRenderer, KeyEvent, ScrollBoxRenderable, Selection } from "@opentui/core"
import { useKeyboard } from "@opentui/solid"
import { createEffect, createSignal, onCleanup, untrack, type JSX } from "solid-js"
import { Pins } from "./agents/pins"
import { WINDOW_COMMANDS, nameFromPrompt, renameOf, resumeOf, shellLead } from "./commands"
import { createSwarmView } from "./swarm/view"
import type { ArgumentChoices } from "./composer/complete"
import { CLICK_PAUSE_MS } from "./mouse"
import { arrivalNote, attentionCounts, markOf } from "./agents/attention"
import type { AgentScreen } from "./agents/screen"
import { createAgentView, needsInput, peekLines, sessionName, type AgentRow } from "./agents/view"
import { createComposer } from "./composer/composer"
import { History } from "./composer/history"
import { ANSWER_PAUSE_MS, createApprovalUi, setFullPreviewKey } from "./footer/approval"
import { Footer } from "./footer/footer"
import { createPicker, type Choice } from "./footer/picker"
import { editPrompt } from "./composer/editor"
import { keyTable, type Keymap } from "./keymap/keymap"
import type { AttentionItem } from "./protocol/generated"
import type { Asked, SessionState } from "./store/state"
import { listed, type AppStore } from "./store/store"
import type { ScreenItems } from "./scrollback/screenview"
import type { Theme } from "./theme/theme"

export interface Requester {
  request(method: string, params?: Record<string, unknown>): Promise<any>
}

export interface AppProps {
  renderer: CliRenderer
  store: AppStore
  keymap: Keymap
  theme: () => Theme
  client: Requester
  /** The prompt history (in memory only when not given). */
  history?: History
  /** A prompt sent with its pastes put back -> the prompt as typed, so the
   *  scrollback shows the placeholder rather than the whole paste. */
  sentAs?: Map<string, string>
  onQuit: () => void
  onRedraw: () => void
  /** The last tool result in full (the palette; the session view does it). */
  onExpand?: () => void
  /** Every tool call's details shown (ctrl+o toggles them). */
  details?: () => boolean
  /** Whether they can be (not under GNU screen, inline). */
  canShowDetails?: () => boolean
  onDetails?: () => void
  /** Show another session in this window -- one that runs in the agent
   *  service, or a saved one (the resume picker, the agent view). */
  onResume?: (sessionId: string, workspace: string) => void
  /** Start a new session and show it (the agent view): its id, or null
   *  when it did not start (the reason is noted). The session in view stays
   *  in view, and takes what is typed, until the new one is up. */
  onNewSession?: () => Promise<string | null>
  /** → on an attention item: opened in this window (its session, a copy,
   *  the diagnosis, a new session about it); the prompt a new one
   *  suggests, for the prompt box -- or null. */
  onOpenItem?: (item: AttentionItem) => Promise<string | null>
  /** Pinned sessions (the agent view's); in memory only when not given. */
  pins?: Pins
  /** The agent view on the alternate screen (agents/screen.ts); without
   *  it the view grows the footer instead. */
  agentScreen?: AgentScreen
  /** Thinking shown or hidden from here on (the palette). */
  onThinking?: (show: boolean) => void
  thinking?: () => boolean
  /** Lines of the window's own, into the scrollback (/keys). */
  onLines?: (title: string, lines: string[]) => void
  /** The pause before a key answers a question (approval.tsx); tests set 0. */
  answerPauseMs?: number
  /** The pause before a click answers a question or picks from a list
   *  (mouse.ts CLICK_PAUSE_MS); tests set 0. */
  clickPauseMs?: number
  /** How long a new item names itself at the bottom right, and for how
   *  long after the window opens nothing is news (tests shorten them). */
  arrival?: { shownMs: number; quietMs: number }
  /** Full screen: the session's rows in the window (session.ts), scrolled
   *  by the wheel and the keys; a selection is copied when the drag ends. */
  screen?: ScreenItems
}

/** How long a new item names itself at the bottom right. */
export const ARRIVAL_SHOWN_MS = 5000
/** Items older than this when they arrive are no news (an inbox read
 *  after a reconnect). */
const ARRIVAL_WINDOW_S = 180
/** What the window finds when it opens is no news. */
const ARRIVAL_QUIET_MS = 2500

/** shift+tab's order, as in the REPL (delfin/agent/posture_cycle.py):
 *  bypassPermissions only where the session says the ladder may stop
 *  there (the managed policy allows it, isolation holds each command). */
export const POSTURES = ["default", "acceptEdits", "plan", "bypassPermissions"] as const

/** The posture after `current`; from one outside the cycle, default.
 *  `bypass` false steps over bypassPermissions; `turnRunning` steps over
 *  plan, as the REPL does (delfin/agent/mode_cycle.py): every step takes
 *  effect at once, and a tool call landing while the cycle stood on plan
 *  was refused. */
export function nextPosture(current: string, bypass = false, turnRunning = false): string {
  const order = POSTURES.filter((m) => (bypass || m !== "bypassPermissions")
    && !(turnRunning && m === "plan"))
  const at = order.indexOf(current as typeof POSTURES[number])
  if (at >= 0) return order[(at + 1) % order.length]!
  // From plan itself while a turn runs: the step after plan on the full ladder.
  const full = POSTURES.filter((m) => bypass || m !== "bypassPermissions")
  const i = full.indexOf(current as typeof POSTURES[number])
  return i < 0 ? "default" : full[(i + 1) % full.length]!
}

/** What the confirm row says while the first step onto bypass waits;
 *  BYPASS_YES is its clickable part (full screen). */
export const BYPASS_YES = "[y] switch"
export function bypassOffer(from: string): string {
  return `⏵⏵ bypass permissions asks nothing from here on: ${BYPASS_YES} · any other key stays in ${from}`
}

export function App(props: AppProps): JSX.Element {
  const store = props.store
  const session = (): SessionState | null => {
    const id = store.state.current
    return id ? store.state.sessions[id] ?? null : null
  }
  const [message, setMessage] = createSignal("")
  /** The task list folded to one line (ctrl+t), for every session in view. */
  const [tasksFolded, setTasksFolded] = createSignal(false)
  let quitArmed = 0

  const say = (text: string) => {
    setMessage(text)
    // Unref'd: a note fading out must not keep a window that left alive.
    if (text) (setTimeout(() => message() === text && setMessage(""), 6000) as { unref?: () => void }).unref?.()
  }
  const call = async (method: string, params: Record<string, unknown>) => {
    try {
      return await props.client.request(method, params)
    } catch (error) {
      say((error as Error).message)
      return undefined
    }
  }
  const sid = () => store.state.current ?? ""

  // -- full screen: scrolling and copying ----------------------------------------
  let scrollBox: ScrollBoxRenderable | null = null
  /** The reader scrolled (keys, the wheel): the next rest repaints. Rows
   *  that only follow the session do not count -- a stream repaints none. */
  let userScrolled = false
  /** PgUp/PgDn a page (less a row, for the eye), ctrl+home/end to the
   *  start or back to the end -- where the box follows the session again. */
  const scroll = (how: "up" | "down" | "top" | "bottom") => {
    const box = scrollBox
    if (!box) return
    userScrolled = true
    const page = Math.max(1, (box.viewport?.height ?? 10) - 1)
    if (how === "top") box.scrollTo(0)
    else if (how === "bottom") box.scrollTo(box.scrollHeight)
    else box.scrollBy(how === "up" ? -page : page)
  }
  /** What a drag selected, to the clipboard -- OSC 52, which reaches the
   *  terminal's clipboard also over SSH; a terminal without it is said. */
  const copy = (selection: Selection | null): boolean => {
    const text = selection?.getSelectedText?.() ?? ""
    if (!text.trim()) return false
    const done = props.renderer.copyToClipboardOSC52(text)
    say(done ? `copied ${text.length} character${text.length === 1 ? "" : "s"}`
      : "this terminal takes no OSC 52 copy: shift+drag selects with its own")
    return true
  }
  if (props.screen) {
    // On release: copied, and the selection let go a moment later (shown
    // that long) -- one left behind would hold the next ctrl+c, and a
    // plain click leaves an empty one.
    let letGo: ReturnType<typeof setTimeout> | null = null
    const onSelection = (selection: Selection) => {
      if (selection?.isDragging) return
      const copied = copy(selection)
      if (letGo) clearTimeout(letGo)
      letGo = setTimeout(() => {
        letGo = null
        if (!props.renderer.isDestroyed) props.renderer.clearSelection()
      }, copied ? 300 : 0)
      ;(letGo as { unref?: () => void }).unref?.()
    }
    props.renderer.on("selection", onSelection)
    // Once the conversation has come to rest after a scroll, every cell is
    // written once more: some terminals (xterm) leave pixels of a block
    // glyph beside where it stood. Only after the reader scrolled.
    let lastTop = -1
    let rest: ReturnType<typeof setTimeout> | null = null
    const onFrame = () => {
      const top = scrollBox?.scrollTop ?? -1
      if (top === lastTop) return
      lastTop = top
      if (!userScrolled) return
      if (rest) clearTimeout(rest)
      rest = setTimeout(() => {
        rest = null
        userScrolled = false
        if (!props.renderer.isDestroyed) props.onRedraw()
      }, 150)
      ;(rest as { unref?: () => void }).unref?.()
    }
    props.renderer.on("frame", onFrame)
    onCleanup(() => {
      props.renderer.off("selection", onSelection)
      props.renderer.off("frame", onFrame)
      if (rest) clearTimeout(rest)
      if (letGo) clearTimeout(letGo)
    })
  }

  /** /rename: what a session is called -- the prompt box, the agent view
   *  and the terminal's title show it; without a name, after what it was
   *  first asked. Done at once, also while a turn runs. */
  const rename = async (id: string, asked: string): Promise<boolean> => {
    const s = store.state.sessions[id]
    const title = asked || nameFromPrompt(String(s?.summary?.firstPrompt ?? ""))
    if (!title) {
      say("nothing to name it after yet: /rename <name>")
      return true
    }
    const result = await call("session.setTitle", { sessionId: id, title })
    if (!result) return false                   // not asked (the note says why): the line stays
    say(`renamed: ${result.title}`)
    return true
  }

  const picker = createPicker({ clickPauseMs: props.clickPauseMs })

  // Prompts sent while a session starts (the first one, or a new one from
  // the agent view): they go out to it once it is up.
  const early: { text: string; shown?: string }[] = []
  const EARLY_NOTE = "the session is starting: this goes out as soon as it is up"
  createEffect(() => {
    if (!sid() || store.state.starting || !early.length) return
    // Gone out: the note saying it would is not true any more.
    if (untrack(message) === EARLY_NOTE) setMessage("")
    for (const { text, shown } of early.splice(0)) void deliver(text, undefined, shown)
  })

  /** What the composer sends: with the agent view open, it is the view's
   *  (a new session, a reply, a name). */
  const send = async (text: string, how?: "steer", shown?: string): Promise<boolean> => {
    if (agents.state.open) {
      typingAt = 0
      // /rename here names the session peeked at, else the one selected --
      // it is never a new session's first prompt (nor a model's turn).
      const named = renameOf(text)
      if (named !== null && !agents.state.renaming) {
        const list = agents.rows()
        const row = agents.state.peek ? list.find((r) => r.id === agents.state.peek!.id) : list[agents.selected(list)]
        if (!row) {
          say("select a session to rename (↑ ↓), then /rename <name>")
          return false
        }
        if (row.kind === "elsewhere") {
          say(`${row.name} runs elsewhere: it is renamed where it runs`)
          return true
        }
        if (row.kind === "attention") {
          say(`this is no session: ${props.keymap.label("agentsOpen")} opens it, ${props.keymap.label("agentsStop")} dismisses it`)
          return true
        }
        return rename(row.id, named)
      }
      return agents.submit(text)
    }
    if (swarmOpen()) {
      typingAt = 0
      // A typed line in the swarm view: the answer for a held question,
      // else a new swarm's goal (validate -> dry run -> held -> y).
      const body = text.trim()
      if (swarm.state.pendingAnswer) return Boolean(await swarm.confirmAnswer(body))
      if (body) return Boolean(await swarm.newSwarm(body))
      return true
    }
    return deliver(text, how, shown)
  }

  /** A prompt for the session in view (or the one starting). */
  const deliver = async (text: string, how?: "steer", shown?: string): Promise<boolean> => {
    typingAt = 0
    if (!sid() || store.state.starting) {
      if (!store.state.starting) return false
      early.push({ text, shown })
      say(EARLY_NOTE)
      return true
    }
    // A bare /model, /effort, /mode or /resume opens its picker; with an
    // argument the command goes to the session as typed.
    const bare = /^\/(model|effort|mode|resume)$/.exec(text.trim())
    if (bare && !how) {
      openPicker(bare[1]!)
      return true
    }
    // The window's own: the session has no terminal to open an editor in,
    // and the keys are the window's.
    const editor = /^\/editor(?:\s+([\s\S]*))?$/.exec(text.trim())
    if (editor && !how) {
      const outcome = await editPrompt(props.renderer, editor[1] ?? "")
      say(outcome.note)
      return outcome.text ? deliver(outcome.text) : true
    }
    const resumed = resumeOf(text)
    if (resumed !== null && !how) {
      if (resumed === sid()) say("that session is the one in view")
      else props.onResume?.(resumed, workspaces.get(resumed) ?? "")
      return true
    }
    if (/^\/keys$/.test(text.trim()) && !how) {
      props.onLines?.("keys in this window", keyTable(props.keymap))
      return true
    }
    // The swarm view (SW2b): a bare /swarm is the window's, never the
    // session's -- the view manages swarms over the service's swarm.*
    // requests. With an argument the line goes to the session as typed.
    if (/^\/swarm$/.test(text.trim()) && !how) {
      void swarm.open()
      return true
    }
    // Also with the key that steers: a rename is the window's, never the turn's.
    const named = renameOf(text)
    if (named !== null) return rename(sid(), named)
    // The session strips the text it shows, as a typed line is stripped.
    if (shown) props.sentAs?.set(text.trim(), shown)
    const result = await call("session.prompt", how ? { sessionId: sid(), text, how } : { sessionId: sid(), text })
    return result !== undefined
  }

  /** A line from the prompt: with a question of the session in view open,
   *  its typed answer (or the reason of N) -- an approval never takes a
   *  line, which then goes out as a message, as it would during the turn.
   *  A shell line ("!") is a command, never an answer. */
  const sendTyped = async (text: string, how?: "steer", shown?: string): Promise<boolean> => {
    const asked = agents.state.open ? undefined : session()?.pending[0]
    // Session routing must remain switchable while a question is waiting.
    // /sync is a local control, never the text answer to that question.
    if (asked && !how && shellLead(text) < 0 && !/^\/sync(?:\s|$)/i.test(text.trim())) {
      const took = approval.submitText(asked, text)
      if (took === "reflex") {
        say("sent as the question opened: not taken as its answer -- enter again to answer with it")
        return false                       // the line comes back into the prompt
      }
      if (took === "sent") {
        typingAt = 0
        return true
      }
    }
    return send(text, how, shown)
  }

  // The agent view is made after the composer: read through this, set
  // once it exists (a "!" in its prompt is a name, an answer or a message).
  let viewOpen = () => false
  const composer = createComposer({
    shellAllowed: () => !viewOpen(),
    request: (method, params) => props.client.request(method, params),
    sessionId: sid,
    history: props.history ?? new History(),
    keymap: props.keymap,
    say,
    send: sendTyped,
    argumentChoices: (command) => argumentChoices(command),
    clickPauseMs: props.clickPauseMs,
    beforePaste: () => approval.releaseHeld(),
  })

  // -- the agent view ----------------------------------------------------------
  const pins = props.pins ?? new Pins()
  const [pinned, setPinned] = createSignal(0)
  const sendTo = async (row: AgentRow, text: string) => {
    const named = renameOf(text)
    if (named !== null) {
      // Not mailed: a session elsewhere is renamed where it runs.
      if (row.kind === "elsewhere") say(`${row.name} runs elsewhere: it is renamed where it runs`)
      else await rename(row.id, named)
      return
    }
    if (row.kind === "elsewhere") {
      if (await call("session.message", { key: row.id, text })) say(`sent to ${row.name}`)
      return
    }
    const result = await call("session.prompt", { sessionId: row.id, text })
    if (!result) return
    say(result.queued ? `queued for ${row.name}: it goes out when its turn ends` : `sent to ${row.name}`)
  }
  const startNew = async (text: string) => {
    if (store.state.starting) {
      // The first session is still starting: this one follows it (the
      // window then shows this one, the latest asked for).
      say("a new session starts as soon as the first one is up…")
      const until = Date.now() + 180_000
      while (store.state.starting && Date.now() < until) await Bun.sleep(100)
    }
    say("starting a new session…")
    const id = await (props.onNewSession?.() ?? Promise.resolve(null))
    if (id) {
      // Its own prompt, to it -- not to whatever is in view by then.
      await call("session.prompt", { sessionId: id, text })
      return
    }
    // Not started (the note says why): back in the view, so enter tries
    // again -- never sent to the session in view; what was typed meanwhile
    // is kept.
    agents.open()
    composer.setText(composer.value() ? `${composer.value()}\n${text}` : text)
  }
  /** An attention item answered or dismissed: true when the inbox took
   *  it; else the note says why. */
  const settle = async (method: string, params: Record<string, unknown>, what: string): Promise<boolean> => {
    const result = await call(method, params)
    if (!result) return false                 // not asked: the note says why
    if (!result.ok) say(`${what}: ${result.error || "it is not waiting any more"}`)
    return Boolean(result.ok)
  }
  /** → on an item: shown in this window; a new session about it gets its
   *  suggested prompt in the box -- sent only when enter is pressed. */
  const openItem = async (item: AttentionItem) => {
    const suggested = await (props.onOpenItem?.(item) ?? Promise.resolve(null))
    if (suggested && !composer.value().trim() && !agents.state.open) composer.setText(suggested)
  }
  const agents = createAgentView({
    state: () => store.state,
    pins: () => {
      pinned()
      return pins.all
    },
    keymap: props.keymap,
    togglePin: (id) => {
      say(pins.toggle(id) ? "pinned: it stays at the top" : "unpinned")
      setPinned((n) => n + 1)
    },
    open: (row) => {
      if (row.status === "sleeping") say(`waking ${row.name}…`)
      props.onResume?.(row.id, row.workspace)
    },
    peek: async (id) => {
      const got = await props.client.request("session.replay", { sessionId: id, turns: 1 })
      return peekLines((got?.events ?? []) as Record<string, any>[])
    },
    stop: (id) => void call("session.shutdown", { sessionId: id }),
    rename: (id, title) => {
      void call("session.setTitle", { sessionId: id, title }).then((r) => r && say(`renamed: ${r.title}`))
    },
    message: (row, text) => void sendTo(row, text),
    newSession: (text) => void startNew(text),
    openItem: (item) => void openItem(item),
    dismiss: (item) => settle("attention.dismiss", { id: item.id }, "not dismissed"),
    answer: (item, text) => settle("attention.answer", { id: item.id, text }, "not answered"),
    draft: () => composer.value(),
    setDraft: (text) => (text ? composer.setText(text) : composer.clear()),
    say,
  })
  viewOpen = () => agents.state.open
  // The swarm view (SW2b): the window's own, over the service's swarm.*
  // requests. It has the footer while it is open (the agents view's
  // pattern), and its keys there too.
  const swarm = createSwarmView({
    keymap: props.keymap,
    now: () => Date.now() / 1000,
    askStatus: () => props.client.request("swarm.status", {}) as never,
    askValidate: (specYaml) => props.client.request("swarm.validate", { specYaml }) as never,
    askStart: (specYaml, confirm) => props.client.request("swarm.start",
      { specYaml, dryRun: !confirm, confirm }) as never,
    askTriage: () => props.client.request("swarm.triage.list", {}).then((r) => (r as any)?.items ?? []) as never,
    askAnswer: (item, answer) => props.client.request("swarm.triage.answer", { item, answer }) as never,
    note: say,
    draft: () => composer.value(),
  })
  const swarmOpen = () => swarm.state.open
  // The alternate screen, while the view is open.
  createEffect(() => {
    if (!props.agentScreen) return
    if (agents.state.open) void props.agentScreen.enter()
    else props.agentScreen.leave()
  })
  // A peeked session that moves on is peeked again (at most twice a second).
  let peekTimer: ReturnType<typeof setTimeout> | undefined
  createEffect(() => {
    const id = agents.state.peek?.id
    const s = id ? store.state.sessions[id] : undefined
    void s?.summary?.lastActivity
    void s?.state
    if (!s || peekTimer) return
    peekTimer = setTimeout(() => {
      peekTimer = undefined
      agents.refreshPeek()
    }, 500)
    ;(peekTimer as { unref?: () => void }).unref?.()
  })
  // While the view is open: the other DELFIN sessions (dashboards, other
  // terminals) and any session this window has not heard of, every 5 s.
  createEffect(() => {
    if (!agents.state.open) return
    let gone = false
    const poll = async () => {
      const got = await props.client.request("sessions.list", {}).catch(() => null)
      if (gone || !got) return
      // Every listed session as the service has it; one no longer listed
      // has ended (a summary missed while disconnected, a restart).
      listed(store, (got.sessions ?? []) as { sessionId?: string }[])
      store.update((d) => {
        d.elsewhere = ((got.elsewhere ?? []) as any[]).filter((e) => e?.key).map((e) => ({
          key: String(e.key), sessionId: e.sessionId ?? null, title: e.title ?? null,
          workspace: String(e.workspace ?? ""), host: String(e.host ?? ""),
        }))
      })
    }
    void poll()
    const timer = setInterval(() => void poll(), 5000)
    ;(timer as { unref?: () => void }).unref?.()
    onCleanup(() => {
      gone = true
      clearInterval(timer)
    })
  })
  /** Sessions other than the one in view that wait for the user. */
  const othersWaiting = () => Object.values(store.state.sessions)
    .filter((s) => s.id !== store.state.current && needsInput(s)).length
  /** The other sessions the service runs (a worker each: none asleep or
   *  ended), and how many of them wait -- for the hints line. Live: every
   *  window hears every session's summary as it changes. */
  const others = () => {
    const running = Object.values(store.state.sessions)
      .filter((s) => s.id !== store.state.current && s.state && s.state !== "ended" && s.state !== "sleeping")
    return { running: running.length, waiting: running.filter(needsInput).length }
  }
  // -- what needs the user, seen from anywhere in the window ----------------------
  // Quiet, at the bottom right of the hints line: what waits for an answer
  // (other sessions here, questions parked anywhere) and what failed. A new
  // one names itself there for a few seconds instead of the count.
  const [arrival, setArrival] = createSignal("")
  let arrivalTimer: ReturnType<typeof setTimeout> | undefined
  const arrive = (text: string) => {
    setArrival(text)
    if (arrivalTimer) clearTimeout(arrivalTimer)
    arrivalTimer = setTimeout(() => setArrival(""), props.arrival?.shownMs ?? ARRIVAL_SHOWN_MS)
    ;(arrivalTimer as { unref?: () => void }).unref?.()
  }
  onCleanup(() => arrivalTimer && clearTimeout(arrivalTimer))
  // What the window starts with is no news: only what comes after.
  const openedAt = Date.now()
  const news = () => Date.now() - openedAt > (props.arrival?.quietMs ?? ARRIVAL_QUIET_MS)
  const known = new Set<string>()
  createEffect(() => {
    const now = Date.now() / 1000
    const fresh: typeof store.state.attention = []
    for (const item of store.state.attention) {
      if (known.has(item.id)) continue
      known.add(item.id)
      // A notice is quiet; an old item the window had not heard of (a
      // reconnect, a service started again) is no news either.
      if (markOf(item) === "notice" || !news() || (item.createdAt && now - item.createdAt > ARRIVAL_WINDOW_S)) continue
      fresh.push(item)
    }
    // Of what came together, the most urgent (the list's order: what
    // waits first), and how many more.
    if (fresh.length) arrive(arrivalNote(fresh[0]!) + (fresh.length > 1 ? ` (+${fresh.length - 1} more)` : ""))
  })
  // Another session that starts to wait (a question, an approval) says so
  // once, in the message row, for as long as an arrival is shown: which
  // one, and the key that opens the agent view (in the view: where it is).
  const waitingBefore = new Set<string>()
  const announce = (text: string) => {
    setMessage(text)
    const off = setTimeout(() => message() === text && setMessage(""), props.arrival?.shownMs ?? ARRIVAL_SHOWN_MS)
    ;(off as { unref?: () => void }).unref?.()
  }
  createEffect(() => {
    for (const s of Object.values(store.state.sessions)) {
      if (!needsInput(s)) {
        waitingBefore.delete(s.id)
        continue
      }
      if (waitingBefore.has(s.id)) continue
      waitingBefore.add(s.id)
      if (s.id === store.state.current || !news()) continue
      const where = untrack(() => agents.state.open) ? "it is under NEEDS YOU"
        : `${props.keymap.label("agents")} to open`
      announce(`◉ session '${sessionName(s)}' waits for you — ${where}`)
    }
  })
  const attention = () => ({ ...attentionCounts(store.state, othersWaiting()), note: arrival() })
  // The terminal's title says when a session waits for the user (the
  // only notification every terminal shows, also from another tab), and
  // how much needs the user: (3) DELFIN · name.
  let shownTitle = ""
  createEffect(() => {
    const waiting = Object.values(store.state.sessions).filter(needsInput).length
    const counts = attentionCounts(store.state, waiting)
    const n = counts.asks + counts.failed
    // Named (/rename, -n): the tab says which session it is, as agent CLIs do.
    const name = session()?.title || session()?.summary?.title || ""
    const title = n ? `(${n}) DELFIN${name ? ` · ${name}` : ""}`
      : name ? `DELFIN · ${name}` : shownTitle ? "DELFIN" : ""
    if (title === shownTitle) return
    shownTitle = title
    try {
      props.renderer.setTerminalTitle(title)
    } catch {
      // a terminal without titles
    }
  })

  // -- the palette and the pickers ---------------------------------------------
  const choicesOf = async (what: string, command: string): Promise<Choice[]> => {
    const got = await props.client.request("session.choices", { sessionId: sid(), what })
    return ((got?.choices ?? []) as { value: string; detail?: string }[])
      .filter((c) => c.value)
      .map((c) => ({ label: c.value, detail: c.detail, value: `${command} ${c.value}`, current: c.value === got.current }))
  }
  const age = (secs: number) => {
    const d = Math.max(0, Date.now() / 1000 - secs)
    return d < 3600 ? `${Math.round(d / 60)}m ago` : d < 86400 ? `${Math.round(d / 3600)}h ago` : `${Math.round(d / 86400)}d ago`
  }
  /** Where each listed session works (a saved one is started there). */
  const workspaces = new Map<string, string>()
  /** The sessions running in the service first, then the saved ones of
   *  this workspace (the resume picker, `/resume `'s menu). */
  const resumeRows = (s: SessionState) => Promise.all([
    props.client.request("sessions.list", {}).catch(() => ({ sessions: [] })),
    props.client.request("sessions.past", { workspace: s.workspace, limit: 50 }),
  ]).then(([listed, past]) => {
    const running = ((listed?.sessions ?? []) as any[]).map((r) => ({
      id: String(r.sessionId), workspace: String(r.workspace ?? ""), title: String(r.title || "(untitled)"),
      detail: ["running here", r.state, r.model, String(r.sessionId).slice(0, 8)].filter(Boolean).join(" · "),
      current: r.sessionId === s.id,
    }))
    const saved = ((past?.sessions ?? []) as any[]).filter((r) => !r.running).map((r) => ({
      id: String(r.sessionId), workspace: String(r.workspace ?? ""), title: String(r.title || "(untitled)"),
      detail: [age(r.updatedAt), r.model, `${r.messages} messages`, String(r.sessionId).slice(0, 8)]
        .filter(Boolean).join(" · "),
      current: false,
    }))
    for (const r of [...running, ...saved]) workspaces.set(r.id, r.workspace)
    return [...running, ...saved]
  })
  /** What may follow a command on the prompt: the window's own commands
   *  here, the session's from session.choices ("args"). */
  const argumentChoices = async (command: string): Promise<ArgumentChoices> => {
    const s = session()
    if (command === "/rename") return { current: "", choices: [], usage: "/rename [name]  (no name: after its first prompt)" }
    if (command === "/resume") {
      if (!s) return { current: "", choices: [], usage: "" }
      const rows = await resumeRows(s)
      return { current: s.id, usage: "/resume <session id>",
               choices: rows.map((r) => ({ value: r.id, label: r.title, detail: r.detail })) }
    }
    const got = await props.client.request("session.choices", { sessionId: sid(), what: "args", command })
    return {
      current: String(got?.current ?? ""),
      usage: String(got?.usage ?? ""),
      choices: ((got?.choices ?? []) as any[]).filter((c) => c?.value)
        .map((c) => ({ value: String(c.value), detail: c.detail ? String(c.detail) : "" })),
    }
  }
  const openPicker = (what: string) => {
    const s = session()
    if (!s) return
    switch (what) {
      case "model":
      case "effort":
      case "mode":
        picker.open(what === "mode" ? "agent mode" : what, choicesOf(what, `/${what}`), (c) => void deliver(c.value))
        return
      case "permissions":
        picker.open("approval posture", [
          { label: "plan", detail: "read and plan; every write and command asked", value: "plan" },
          { label: "default", detail: "ask before each write and each command", value: "default" },
          { label: "acceptEdits", detail: "writes go through; commands still asked", value: "acceptEdits" },
        ].map((c) => ({ ...c, current: c.value === s.approval })), (c) => {
          void call("session.setMode", { sessionId: sid(), mode: c.value })
        })
        return
      case "resume":
        picker.open("resume a session", resumeRows(s).then((rows) => rows.map((r) => ({
          label: r.title, detail: r.detail, value: `${r.id}\t${r.workspace}`, current: r.current,
        }))), (c) => {
          const [id, workspace] = c.value.split("\t")
          if (id && id !== sid()) props.onResume?.(id, workspace ?? "")
        })
        return
    }
  }
  const openPalette = () => {
    const actions: Choice[] = [
      { label: "model…", detail: "switch the model", value: "pick:model" },
      { label: "permissions…", detail: "the approval posture (also shift+tab)", value: "pick:permissions" },
      { label: "effort…", detail: "reasoning effort", value: "pick:effort" },
      { label: "mode…", detail: "the agent mode", value: "pick:mode" },
      { label: "resume…", detail: "a saved session, in this window", value: "pick:resume" },
      { label: "rename…", detail: "name this session (/rename)", value: "type:/rename" },
      { label: props.thinking?.() ? "hide thinking" : "show thinking", detail: "the model's thinking, dimmed", value: "do:thinking" },
      ...(props.canShowDetails?.() === false ? [] : [{ label: props.details?.() ? "hide tool details" : "show tool details",
        detail: props.keymap.label("details"), value: "do:details" }]),
      { label: "expand the last result", detail: props.keymap.label("expand") || "the whole output, into the conversation",
        value: "do:expand" },
      { label: tasksFolded() ? "open the task list" : "fold the task list", detail: props.keymap.label("tasks"),
        value: "do:tasks" },
      { label: "redraw", detail: props.keymap.label("redraw"), value: "do:redraw" },
      { label: "leave the window", detail: `${props.keymap.label("quit")} — the session keeps running`, value: "do:quit" },
    ]
    // The window's own command (rename…, above) wins over one of the
    // session's of the same name, as the REPL's builtins win over a user's.
    const windowNames = new Set(WINDOW_COMMANDS.map((c) => c.name))
    const commands = props.client.request("session.commands", { sessionId: sid() })
      .then((got) => ((got?.commands ?? []) as any[]).filter((c) => !windowNames.has(c.name)).map((c) => ({
        label: c.name, detail: c.summary, value: `${c.takesArgs ? "type" : "run"}:${c.name}`,
      })))
    picker.open("commands", actions, (c) => {
      const [kind, what] = [c.value.slice(0, c.value.indexOf(":")), c.value.slice(c.value.indexOf(":") + 1)]
      if (kind === "pick") openPicker(what)
      else if (kind === "run") void deliver(what)
      // A command with arguments goes in front of what is typed: a draft
      // is never thrown away.
      else if (kind === "type") composer.setText(composer.value().trim() ? `${what} ${composer.value()}` : `${what} `)
      else if (what === "thinking") props.onThinking?.(!props.thinking?.())
      else if (what === "expand") props.onExpand?.()
      else if (what === "details") props.onDetails?.()
      else if (what === "tasks") setTasksFolded(!tasksFolded())
      else if (what === "redraw") props.onRedraw()
      else if (what === "quit") props.onQuit()
    }, commands)
  }

  const reply = async (asked: Asked, params: Record<string, unknown>): Promise<boolean> => {
    const result = await call("permission.reply", { sessionId: sid(), requestId: asked.requestId, ...params })
    const typed = typeof params.text === "string" ? params.text : null
    if (result && result.applied === false) say("that question had already been answered elsewhere")
    else if (typed !== null && result?.message === "dismissed") {
      // A service from before typed answers reads the line as esc: the
      // question is dismissed, and the line goes out as a message.
      say("this agent service takes no typed answer: the question was dismissed, the line is sent as a message")
      void deliver(typed)
    } else if (result?.message && /once|persist|stop asking|readable|reject/.test(result.message)) say(result.message)
    return result !== undefined
  }
  // When the user last typed into the prompt -- and has not sent it: a
  // question that opens then gets the rest of that typing, not answers.
  // Sending a prompt ends the typing (a question its turn asks a moment
  // later is answered by the next key).
  let typingAt = 0
  setFullPreviewKey(props.keymap.label("details") || props.keymap.label("expand"))
  const approval = createApprovalUi({ reply, say, lastKey: () => typingAt,
                                      pauseMs: props.answerPauseMs, clickPauseMs: props.clickPauseMs,
                                      release: (text) => composer.insert(text),
                                      prefill: (text) => composer.insert(text) })
  // The question's own state starts when it appears (its typeahead pause
  // counts from then), and ends with it.
  createEffect(() => {
    const asked = session()?.pending[0]
    approval.sync(asked)
    if (asked) {
      picker.close()                      // a question takes the footer from a picker
      composer.cancelSearch()             // ... and the keys from a ctrl+r search
    }
  })
  // /exit (or /quit) ends the session, and the window with it, as in the
  // REPL: when it happens -- a session shown that had ended with /exit
  // before (resumed from the agent view) is no reason to leave. A session
  // that ended another way stays in view, saying so.
  let shownId = ""
  let exitedBefore = false
  createEffect(() => {
    const s = session()
    const exited = s?.ended?.reason === "exit"
    if ((s?.id ?? "") !== shownId) {
      shownId = s?.id ?? ""
      exitedBefore = exited
      return
    }
    if (exited && !exitedBefore) props.onQuit()
    exitedBefore = exited
  })

  /** The first step onto bypass in a session waits for y (or a click on
   *  its row): the posture it was on, and when it was offered -- a key
   *  sooner than the answer pause is typing, not an answer. */
  const [offer, setOffer] = createSignal<{ id: string; from: string; at: number } | null>(null)
  // A question takes the footer: the offer is not left behind it.
  createEffect(() => {
    if (session()?.pending[0] || session()?.state === "ended") setOffer(null)
  })
  const cycleMode = () => {
    const s = session()
    if (!s || s.state === "ended") return
    // An offer left open: the key moves on past bypass, from it.
    const from = offer()?.id === s.id ? "bypassPermissions" : s.approval
    setOffer(null)
    const next = nextPosture(from, s.bypassBlocked === "", !!s.turn && s.state !== "ended")
    if (next === "bypassPermissions" && !s.bypassConfirmed) {
      setOffer({ id: s.id, from: s.approval, at: Date.now() })
      return
    }
    void call("session.setMode", { sessionId: sid(), mode: next })
  }
  /** A click on the offer's [y] switch: as y, after CLICK_PAUSE_MS from
   *  the offer (a double shift+tab cannot be a click; a hand already on
   *  its way to where the row appeared is not an answer). */
  const clickOffer = (): boolean => {
    const o = offer()
    if (!o || o.id !== sid()) return false
    if (Date.now() - o.at < (props.clickPauseMs ?? CLICK_PAUSE_MS)) return false
    setOffer(null)
    void call("session.setMode", { sessionId: sid(), mode: "bypassPermissions", confirmed: true })
    return true
  }
  /** y on the offer: the switch, marked as confirmed by a person. */
  const takeOffer = (key: KeyEvent): boolean => {
    const o = offer()
    if (!o) return false
    // shift+tab again: cycleMode moves on past bypass from the offer.
    if (o.id === sid() && props.keymap.matches("cycleMode", key)) return false
    setOffer(null)
    if (o.id !== sid()) return false
    const deliberate = Date.now() - o.at >= (props.answerPauseMs ?? ANSWER_PAUSE_MS)
    if (deliberate && !key.ctrl && !key.meta && (key.name === "y" || key.sequence === "y")) {
      void call("session.setMode", { sessionId: sid(), mode: "bypassPermissions", confirmed: true })
      return true
    }
    // Any other key: the posture stays, and the key does what it does.
    say(`still ${o.from}`)
    return false
  }

  const draft = () => composer.value()

  /** ctrl+c: end the turn; else clear what is typed; else, twice, leave. */
  const ladder = (s: SessionState | null) => {
    if (s?.turn && s.state !== "ended") {
      void call("session.cancel", { sessionId: sid() })
      return
    }
    if (draft()) {
      composer.clear()
      return
    }
    const now = Date.now()
    if (now - quitArmed < 2000) {
      props.onQuit()
      return
    }
    quitArmed = now
    say("ctrl+c again to leave — the sessions keep running in the agent service")
  }

  /** Whether this key reached the open question (approval.handleKey). */
  let reachedQuestion = false
  const onKey = (key: KeyEvent): void => {
    reachedQuestion = false
    try {
      route(key)
    } finally {
      // ctrl+c, the agent view, ctrl+d: a key the question held answers nothing.
      if (!reachedQuestion) approval.dropHeld()
    }
  }
  const route = (key: KeyEvent): void => {
    const s = session()
    const km = props.keymap
    const asked = s?.pending[0]
    if (props.screen) {
      // The conversation's own scroll keys (a drag copies on release, so
      // ctrl+c keeps its ladder).
      if (!agents.state.open) {
        for (const [action, how] of [["scrollUp", "up"], ["scrollDown", "down"], ["scrollTop", "top"],
                                     ["scrollBottom", "bottom"]] as const) {
          if (km.matches(action, key)) {
            key.preventDefault()
            scroll(how)
            return
          }
        }
      }
    }
    const typing = key.sequence && key.sequence.length >= 1 && !key.ctrl && !key.meta
      && !km.matches("submit", key) && key.name !== "escape"
    if (!asked && typing) typingAt = Date.now()
    if (agents.state.open) {
      // The view has the footer: its keys first, then the composer
      // (letters are a new session's prompt). A question of the session
      // in view waits; its row says so.
      if (key.ctrl && key.name === "c" && !composer.search()) {
        key.preventDefault()
        if (agents.cancelRename()) return
        if (draft()) composer.clear()
        else agents.close()
        return
      }
      if (km.matches("quit", key) && !draft()) {
        key.preventDefault()
        props.onQuit()
        return
      }
      if (agents.handleKey(key) || composer.handleKey(key)) key.preventDefault()
      return
    }
    if (swarmOpen()) {
      // The swarm view has the footer: its keys, then the composer
      // (letters type the answer for a held question, or a new goal).
      if (key.ctrl && key.name === "c" && !composer.search()) {
        key.preventDefault()
        if (draft()) composer.clear()
        else swarm.close()
        return
      }
      if (km.matches("quit", key) && !draft()) {
        key.preventDefault()
        props.onQuit()
        return
      }
      if (swarm.handleKey(key) || composer.handleKey(key)) key.preventDefault()
      return
    }
    if (offer() && !asked && takeOffer(key)) {
      key.preventDefault()
      return
    }
    if (km.matches("agents", key) && !composer.value() && !composer.search() && !picker.state.open
        && !(asked && approval.state.mode === "reason")) {
      key.preventDefault()
      agents.open()
      return
    }
    if (key.ctrl && key.name === "c" && !composer.search()) {
      key.preventDefault()
      // An open picker or palette is what ctrl+c leaves first.
      if (picker.state.open) {
        picker.close()
        return
      }
      ladder(s)
      return
    }
    if (asked) {
      picker.close()                       // the question holds the turn: no picker beside it
      // ctrl+d on an empty prompt leaves the window; the question stays
      // open for any window.
      if (km.matches("quit", key) && !draft()) {
        key.preventDefault()
        props.onQuit()
        return
      }
      // ctrl+o: the whole preview, into the scrollback (the panel shows 24
      // lines, fewer on a short terminal).
      if (km.matches("details", key) || km.matches("expand", key)) {
        key.preventDefault()
        const what = asked.kind === "plan" ? "the plan" : asked.kind === "ask" ? "the question" : `${asked.tool}, in full`
        const body = asked.kind === "ask"
          ? [String(asked.question ?? ""), ...(asked.choices ?? []).map((c, i) => `${i + 1}. ${c}`)]
          : String(asked.preview ?? "").split("\n")
        props.onLines?.(what, body)
        return
      }
      // The prompt keeps the keys: the question's own from an empty prompt
      // are its answer; the rest is typing (approval.tsx).
      reachedQuestion = true
      if (approval.handleKey(asked, key, !draft() && !composer.search())) {
        key.preventDefault()
        return
      }
      if (typing) typingAt = Date.now()
      if (composer.handleKey(key)) {
        key.preventDefault()
        return
      }
      // The window's own actions wait for the answer: esc never ends the
      // turn behind a question, shift+tab never moves the posture under it.
      if ((["interrupt", "steer", "cycleMode", "palette", "redraw", "tasks"] as const)
        .some((action) => km.matches(action, key))) key.preventDefault()
      return
    }
    if (picker.state.open) {
      key.preventDefault()
      picker.handleKey(key)
      return
    }
    if (km.matches("palette", key)) {
      key.preventDefault()
      if (s && s.state !== "ended") openPalette()
      return
    }
    if (composer.handleKey(key)) {
      key.preventDefault()
      return
    }
    if (km.matches("quit", key) && !draft()) {
      key.preventDefault()
      props.onQuit()
      return
    }
    if (km.matches("interrupt", key) && s?.turn && s.state !== "ended") {
      key.preventDefault()
      void call("session.cancel", { sessionId: sid() })
      return
    }
    if (km.matches("cycleMode", key)) {
      key.preventDefault()
      cycleMode()
      return
    }
    if (km.matches("redraw", key)) {
      key.preventDefault()
      props.onRedraw()
      return
    }
    if (km.matches("details", key)) {
      key.preventDefault()
      props.onDetails?.()
      return
    }
    if (km.matches("expand", key)) {
      key.preventDefault()
      props.onExpand?.()
      return
    }
    if (km.matches("tasks", key)) {
      key.preventDefault()
      setTasksFolded(!tasksFolded())
      return
    }
    if (km.matches("steer", key)) {
      key.preventDefault()
      // Into the running turn; with no turn running it is just sent.
      if (draft().trim()) composer.submit(s?.turn ? "steer" : undefined)
    }
  }
  useKeyboard(onKey)

  return (
    <Footer
      renderer={props.renderer}
      theme={props.theme}
      keymap={props.keymap}
      session={session}
      connection={() => store.state.connection}
      starting={() => store.state.starting}
      message={message}
      composer={composer}
      approval={approval}
      picker={picker}
      tasksFolded={tasksFolded}
      offer={() => offer()?.id === sid() ? bypassOffer(offer()!.from) : ""}
      offerYes={BYPASS_YES}
      onOfferClick={() => void clickOffer()}
      agents={agents}
      swarm={swarm}
      agentsOnOwnScreen={Boolean(props.agentScreen)}
      attention={attention}
      others={others}
      onOthersClick={() => {
        if (!agents.state.open) agents.open()
      }}
      focused={() => !picker.state.open}
      screen={props.screen}
      scrollRef={(box) => { scrollBox = box }}
      onUserScroll={() => { userScrolled = true }}
    />
  )
}
