/**
 * The window: one session in view, its scrollback, and the footer.
 *
 * Keys go to the open question first (it holds the turn; approval.tsx),
 * then to the composer (completion, history, search), then to the keymap's
 * actions. ctrl+c and ctrl+d work in all of them. Leaving the window
 * (ctrl+d, or ctrl+c twice) never ends a session: it keeps running in the
 * agent service, and the next window attaches to it again.
 */
import type { CliRenderer, KeyEvent } from "@opentui/core"
import { useKeyboard } from "@opentui/solid"
import { createEffect, createSignal, type JSX } from "solid-js"
import { createComposer } from "./composer/composer"
import { History } from "./composer/history"
import { createApprovalUi } from "./footer/approval"
import { Footer } from "./footer/footer"
import { createPicker, type Choice } from "./footer/picker"
import { editPrompt } from "./composer/editor"
import { keyTable, type Keymap } from "./keymap/keymap"
import type { Asked, SessionState } from "./store/state"
import type { AppStore } from "./store/store"
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
  /** ctrl+o: the last tool result in full (the session view does it). */
  onExpand?: () => void
  /** Show a saved session in this window (the resume picker). */
  onResume?: (sessionId: string, workspace: string) => void
  /** Thinking shown or hidden from here on (the palette). */
  onThinking?: (show: boolean) => void
  thinking?: () => boolean
  /** Lines of the window's own, into the scrollback (/keys). */
  onLines?: (title: string, lines: string[]) => void
  /** The pause before a key answers a question (approval.tsx); tests set 0. */
  answerPauseMs?: number
}

/** shift+tab's order, as in the REPL (never onto bypass). */
export const POSTURES = ["plan", "default", "acceptEdits"] as const

/** The posture after `current`; from one outside the cycle, default. */
export function nextPosture(current: string): string {
  const at = POSTURES.indexOf(current as typeof POSTURES[number])
  return at < 0 ? "default" : POSTURES[(at + 1) % POSTURES.length]!
}

export function App(props: AppProps): JSX.Element {
  const store = props.store
  const session = (): SessionState | null => {
    const id = store.state.current
    return id ? store.state.sessions[id] ?? null : null
  }
  const [message, setMessage] = createSignal("")
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

  const picker = createPicker()

  const send = async (text: string, how?: "steer", shown?: string): Promise<boolean> => {
    if (!sid()) return false
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
      return outcome.text ? send(outcome.text) : true
    }
    if (/^\/keys$/.test(text.trim()) && !how) {
      props.onLines?.("keys in this window", keyTable(props.keymap))
      return true
    }
    // The session strips the text it shows, as a typed line is stripped.
    if (shown) props.sentAs?.set(text.trim(), shown)
    const result = await call("session.prompt", how ? { sessionId: sid(), text, how } : { sessionId: sid(), text })
    return result !== undefined
  }

  const composer = createComposer({
    request: (method, params) => props.client.request(method, params),
    sessionId: sid,
    history: props.history ?? new History(),
    keymap: props.keymap,
    say,
    send,
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
  const openPicker = (what: string) => {
    const s = session()
    if (!s) return
    switch (what) {
      case "model":
      case "effort":
      case "mode":
        picker.open(what === "mode" ? "agent mode" : what, choicesOf(what, `/${what}`), (c) => void send(c.value))
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
        // The sessions running in the service first, then the saved ones of
        // this workspace (the agent view will show the running ones too).
        picker.open("resume a session", Promise.all([
          props.client.request("sessions.list", {}).catch(() => ({ sessions: [] })),
          props.client.request("sessions.past", { workspace: s.workspace, limit: 50 }),
        ]).then(([listed, past]) => {
          const running = ((listed?.sessions ?? []) as any[]).map((r) => ({
            label: r.title || "(untitled)",
            detail: ["running here", r.state, r.model, String(r.sessionId).slice(0, 8)].filter(Boolean).join(" · "),
            value: `${r.sessionId}\t${r.workspace}`,
            current: r.sessionId === s.id,
          }))
          const saved = ((past?.sessions ?? []) as any[]).filter((r) => !r.running).map((r) => ({
            label: r.title || "(untitled)",
            detail: [age(r.updatedAt), r.model, `${r.messages} messages`, String(r.sessionId).slice(0, 8)]
              .filter(Boolean).join(" · "),
            value: `${r.sessionId}\t${r.workspace}`,
          }))
          return [...running, ...saved]
        }), (c) => {
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
      { label: props.thinking?.() ? "hide thinking" : "show thinking", detail: "the model's thinking, dimmed", value: "do:thinking" },
      { label: "expand the last result", detail: props.keymap.label("expand"), value: "do:expand" },
      { label: "task list on/off", detail: props.keymap.label("tasks"), value: "do:tasks" },
      { label: "redraw", detail: props.keymap.label("redraw"), value: "do:redraw" },
      { label: "leave the window", detail: `${props.keymap.label("quit")} — the session keeps running`, value: "do:quit" },
    ]
    const commands = props.client.request("session.commands", { sessionId: sid() })
      .then((got) => ((got?.commands ?? []) as any[]).map((c) => ({
        label: c.name, detail: c.summary, value: `${c.takesArgs ? "type" : "run"}:${c.name}`,
      })))
    picker.open("commands", actions, (c) => {
      const [kind, what] = [c.value.slice(0, c.value.indexOf(":")), c.value.slice(c.value.indexOf(":") + 1)]
      if (kind === "pick") openPicker(what)
      else if (kind === "run") void send(what)
      else if (kind === "type") composer.setText(`${what} `)
      else if (what === "thinking") props.onThinking?.(!props.thinking?.())
      else if (what === "expand") props.onExpand?.()
      else if (what === "tasks") void call("session.toggleTasks", { sessionId: sid() })
      else if (what === "redraw") props.onRedraw()
      else if (what === "quit") props.onQuit()
    }, commands)
  }

  const reply = async (asked: Asked, params: Record<string, unknown>): Promise<boolean> => {
    const result = await call("permission.reply", { sessionId: sid(), requestId: asked.requestId, ...params })
    if (result && result.applied === false) say("that question had already been answered elsewhere")
    else if (result?.message && /once|persist|stop asking|readable|reject/.test(result.message)) say(result.message)
    return result !== undefined
  }
  const approval = createApprovalUi({ reply, say,
                                      pauseMs: props.answerPauseMs })
  // The question's own state starts when it appears (its typeahead pause
  // counts from then), and ends with it.
  createEffect(() => {
    const asked = session()?.pending[0]
    approval.sync(asked)
    if (asked) picker.close()             // a question takes the footer from a picker
  })
  // /exit (or /quit) ends the session, and the window with it, as in the
  // REPL. A session that ended another way stays in view, saying so.
  createEffect(() => {
    if (session()?.ended?.reason === "exit") props.onQuit()
  })

  const cycleMode = () => {
    const s = session()
    if (!s || s.state === "ended") return
    void call("session.setMode", { sessionId: sid(), mode: nextPosture(s.approval) })
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

  const onKey = (key: KeyEvent): void => {
    const s = session()
    const km = props.keymap
    const asked = s?.pending[0]
    if (key.ctrl && key.name === "c" && !composer.search()) {
      key.preventDefault()
      ladder(s)
      return
    }
    if (asked) {
      key.preventDefault()
      picker.close()                       // the question holds the turn: it has the keys
      // ctrl+d leaves the window; the question stays open for any window.
      if (km.matches("quit", key)) {
        props.onQuit()
        return
      }
      approval.handleKey(asked, key)
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
    if (km.matches("expand", key)) {
      key.preventDefault()
      props.onExpand?.()
      return
    }
    if (km.matches("tasks", key)) {
      key.preventDefault()
      if (s && s.state !== "ended") void call("session.toggleTasks", { sessionId: sid() })
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
      message={message}
      composer={composer}
      approval={approval}
      picker={picker}
      focused={() => !picker.state.open}
    />
  )
}
