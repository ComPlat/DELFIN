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
import type { Keymap } from "./keymap/keymap"
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

  const composer = createComposer({
    request: (method, params) => props.client.request(method, params),
    sessionId: sid,
    history: props.history ?? new History(),
    keymap: props.keymap,
    say,
    send: async (text, how, shown) => {
      if (!sid()) return false
      // The session strips the text it shows, as a typed line is stripped.
      if (shown) props.sentAs?.set(text.trim(), shown)
      const result = await call("session.prompt", how ? { sessionId: sid(), text, how } : { sessionId: sid(), text })
      return result !== undefined
    },
  })

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
  createEffect(() => approval.sync(session()?.pending[0]))

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
      // ctrl+d leaves the window; the question stays open for any window.
      if (km.matches("quit", key)) {
        props.onQuit()
        return
      }
      approval.handleKey(asked, key)
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
      focused={() => true}
    />
  )
}
