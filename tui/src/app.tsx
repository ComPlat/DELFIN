/**
 * The window: one session in view, its scrollback, and the footer.
 *
 * Keys go to the open question first (it holds the turn), then to the
 * keymap's actions; everything else is the composer's. Leaving the window
 * (ctrl+d, or ctrl+c twice) never ends a session: it keeps running in the
 * agent service, and the next window attaches to it again.
 */
import type { CliRenderer, InputRenderable, KeyEvent } from "@opentui/core"
import { useKeyboard } from "@opentui/solid"
import { createSignal, type JSX } from "solid-js"
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
  onQuit: () => void
  onRedraw: () => void
}

const POSTURES = ["plan", "default", "acceptEdits"] as const

function keyChar(key: KeyEvent): string {
  return key.sequence && key.sequence.length === 1 ? key.sequence : key.name
}

export function App(props: AppProps): JSX.Element {
  const store = props.store
  const session = (): SessionState | null => {
    const id = store.state.current
    return id ? store.state.sessions[id] ?? null : null
  }
  const [message, setMessage] = createSignal("")
  const [persistFor, setPersistFor] = createSignal<{ asked: Asked; choice: string } | null>(null)
  const [picked, setPicked] = createSignal<number[]>([])
  let input: InputRenderable | undefined
  let quitArmed = 0

  const say = (text: string) => {
    setMessage(text)
    if (text) setTimeout(() => message() === text && setMessage(""), 6000)
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

  const reply = async (asked: Asked, params: Record<string, unknown>) => {
    const result = await call("permission.reply", { sessionId: sid(), requestId: asked.requestId, ...params })
    if (result && result.applied === false) say("that question had already been answered elsewhere")
    else if (result?.message && /once|persist|stop asking|readable/.test(result.message)) say(result.message)
  }

  const onApprovalKey = (asked: Asked, key: KeyEvent): void => {
    const pending = persistFor()
    if (pending) {
      setPersistFor(null)
      void reply(pending.asked, { choice: pending.choice, persist: keyChar(key) === "y" })
      return
    }
    if (asked.kind === "ask") {
      if (key.name === "escape") {
        setPicked([])
        void reply(asked, { answers: [] })
        return
      }
      const n = Number(keyChar(key))
      const choices = asked.choices ?? []
      if (Number.isInteger(n) && n >= 1 && n <= choices.length) {
        if (!asked.multiSelect) {
          void reply(asked, { answers: [choices[n - 1]] })
          return
        }
        setPicked((p) => (p.includes(n) ? p.filter((x) => x !== n) : [...p, n].sort()))
        say(`picked: ${picked().map((i) => choices[i - 1]).join(", ") || "nothing"} (enter sends)`)
        return
      }
      if (asked.multiSelect && key.name === "return") {
        const answers = picked().map((i) => choices[i - 1])
        setPicked([])
        void reply(asked, { answers })
      }
      return
    }
    const char = keyChar(key)
    const option = asked.options.find((o) => o.key === char)
    if (!option) return
    if (char === "?") {
      void props.client.request("permission.reply", { sessionId: sid(), requestId: asked.requestId, choice: "?" })
        .catch((error) => say(String((error as Error).message).split("\n").slice(0, 3).join(" ")))
      return
    }
    if (char === "A" || char === "k") {
      // A rule that is global and permanent takes a second yes, as in the REPL.
      setPersistFor({ asked, choice: char })
      const pattern = char === "k" ? asked.suggestion : "this exact command"
      say(`persist ${pattern ?? ""} for every future session? [y] persist  [any other key] just this once`)
      return
    }
    void reply(asked, { choice: char })
  }

  const cycleMode = () => {
    const s = session()
    if (!s) return
    const next = POSTURES[(POSTURES.indexOf(s.approval as typeof POSTURES[number]) + 1) % POSTURES.length]!
    void call("session.setMode", { sessionId: sid(), mode: next })
  }

  const draft = () => input?.value ?? ""

  const onKey = (key: KeyEvent): void => {
    const s = session()
    const km = props.keymap
    const asked = s?.pending[0]
    if (asked) {
      key.preventDefault()
      onApprovalKey(asked, key)
      return
    }
    if (key.ctrl && key.name === "c") {
      key.preventDefault()
      if (s?.turn) {
        void call("session.cancel", { sessionId: sid() })
        return
      }
      const now = Date.now()
      if (now - quitArmed < 2000) {
        props.onQuit()
        return
      }
      quitArmed = now
      say("ctrl+c again to leave — the sessions keep running in the agent service")
      return
    }
    if (km.matches("quit", key) && !draft()) {
      key.preventDefault()
      props.onQuit()
      return
    }
    if (km.matches("interrupt", key) && s?.turn) {
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
    if (km.matches("steer", key) && draft().trim()) {
      key.preventDefault()
      const text = draft()
      if (input) input.value = ""
      void call("session.prompt", { sessionId: sid(), text, how: "steer" })
    }
  }
  useKeyboard(onKey)

  const submit = (value: string) => {
    const text = value.trim()
    if (input) input.value = ""
    if (!text || !sid()) return
    void call("session.prompt", { sessionId: sid(), text })
  }

  return (
    <Footer
      renderer={props.renderer}
      theme={props.theme}
      keymap={props.keymap}
      session={session}
      connection={() => store.state.connection}
      message={message}
      onSubmit={submit}
      inputRef={(el) => (input = el)}
      focused={() => true}
    />
  )
}
