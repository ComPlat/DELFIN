/**
 * The composer: the prompt the user types into, at the bottom of the
 * footer.
 *
 *   enter          send (with a completion open: take it)
 *   shift/alt+enter, ctrl+j, a trailing "\" + enter    a new line
 *   ↑ / ↓          the history, from the first / last line of the prompt
 *   ctrl+r         search the history (again: older; enter takes it)
 *   tab            complete a /command, its argument or an @path
 *   click          on a row of the list: take it, as enter does (full
 *                  screen; not within CLICK_PAUSE_MS of the list changing)
 *   esc            close the completion list
 *   !              at the start of the prompt: shell mode (below)
 *
 * Shell mode: a line that starts with "!" is not a prompt for the model --
 * the session runs it as a command, through the same deny-list, secret
 * scan, allow-list and approval as the model's own bash calls
 * (delfin/agent/repl.py `_handle_line`, `SHELL_PREFIX`). The composer
 * takes the "!" out of the text and shows it as the prompt marker
 * instead, with no completion (nothing of a /command or an @path applies
 * to a shell line); backspace at the start of the line, or esc, puts the
 * line back to a prompt. What is sent is "!" + the line, as typed, so the
 * session sees exactly what it saw before the mode existed. A shell line
 * is never steered into a running turn: the model would read it as text.
 *
 * Large pastes stand in as `[Pasted text #1 +42 lines]` (paste.ts) and go
 * out whole. `createComposer` holds the state and the key handling, so the
 * app can route keys to it in the right order (an open question first);
 * `ComposerView` draws it.
 */
import { decodePasteBytes, stripAnsiSequences, TextAttributes, type ColorInput, type KeyEvent,
         type TextareaRenderable } from "@opentui/core"
import { usePaste } from "@opentui/solid"
import { For, Show, createEffect, createSignal, on, type JSX } from "solid-js"
import type { Keymap } from "../keymap/keymap"
import { fit, width } from "../text"
import type { Theme } from "../theme/theme"
import {
  applyCandidate, argumentAt, matchArguments, matchCommands, pathCandidates, tokenAt,
  type ArgumentChoices, type Candidate, type Command, type Token,
} from "./complete"
import { CLICK_PAUSE_MS, clickable } from "../mouse"
import { HistoryWalk, type Entry, type History } from "./history"
import { Pastes, normalisePaste, shouldCollapse } from "./paste"
import { SHELL_PREFIX, WINDOW_COMMANDS, shellLead } from "../commands"

/** The prompt grows to this many rows, then scrolls. */
export const COMPOSER_MAX_ROWS = 8
/** Completion rows shown at once. The window's own commands grew from 2
 *  to 3 (/swarm, SW2b): a session's 4 common commands would otherwise
 *  lose the 7th row off the visible menu. */
export const POPUP_ROWS = 7

export interface ComposerDeps {
  request(method: string, params?: Record<string, unknown>): Promise<any>
  sessionId(): string
  history: History
  keymap: Keymap
  /** A prompt to send; `how: "steer"` sends it into the running turn.
   *  `shown` is the prompt as typed, when pastes were put back into it.
   *  False (or a promise of it) when it could not be sent: the prompt then
   *  comes back into the composer, if nothing else was typed meanwhile. */
  send(text: string, how?: "steer", shown?: string): Promise<boolean> | boolean | void
  say(text: string): void
  /** What may follow `command` (the app's: the window's own commands, else
   *  the session's session.choices "args"). Without it no argument list. */
  argumentChoices?(command: string): Promise<ArgumentChoices>
  /** The pause before a click on the list is taken (tests: 0). */
  clickPauseMs?: number
  /** Whether a leading "!" means shell mode here (not while the prompt is
   *  for something else: the agent view's names, answers and messages).
   *  Always, when not given. */
  shellAllowed?: () => boolean
  /** Just before a paste goes into the prompt (a key held by an open
   *  question goes in first: approval.tsx HOLD_MS). */
  beforePaste?: () => void
}

interface Popup {
  token: Token
  candidates: Candidate[]
  selected: number
}

interface Search {
  query: string
  at: number            // the history index shown, -1 for none
  draft: Entry
}

export interface Composer {
  /** The prompt as drawn (a signal: it follows the textarea's change
   *  events, which arrive just after the edit). */
  readonly text: () => string
  /** The prompt right now, for decisions made in a key handler. */
  value(): string
  readonly rows: () => number
  readonly popup: () => Popup | null
  readonly search: () => Search | null
  /** The line is a shell command ("!"): the marker, the border and the
   *  missing completion say so. */
  readonly shell: () => boolean
  /** Rows the composer takes in the footer, completion and search included. */
  readonly height: () => number
  attach(textarea: TextareaRenderable): void
  /** A key, before the textarea sees it; true when the composer took it. */
  handleKey(key: KeyEvent): boolean
  /** Paste bytes, before the textarea sees them; true when collapsed. */
  handlePaste(text: string): boolean
  submit(how?: "steer"): boolean
  clear(): void
  setText(text: string, pastes?: Record<string, string>): void
  /** Text typed at the cursor, as if by keys (a key an open question held). */
  insert(text: string): void
  /** Re-measure after a width change. */
  remeasure(): void
  /** Leave a ctrl+r search, the draft put back (a question opened). */
  cancelSearch(): void
  /** A click on the list's row `index` (of the candidates): as enter on
   *  it; false when not taken (too soon, a usage line, no list). */
  clickAt(index: number): boolean
}

const isPlain = (key: KeyEvent) => !key.ctrl && !key.meta && !key.option && !(key as { super?: boolean }).super

export function createComposer(deps: ComposerDeps): Composer {
  let area: TextareaRenderable | undefined
  const [text, setTextSignal] = createSignal("")
  const [shell, setShell] = createSignal(false)
  const shellAllowed = () => deps.shellAllowed?.() ?? true
  const [rows, setRows] = createSignal(1)
  const [popup, setPopupSignal] = createSignal<Popup | null>(null)
  // When the list's rows last changed: a click that soon is not taken.
  let popupShownAt = 0
  const setPopup = (next: Popup | null) => {
    const before = popup()
    if (!next || !before || next.candidates !== before.candidates) popupShownAt = Date.now()
    setPopupSignal(next)
  }
  const argsCache = new Map<string, { at: number; got: ArgumentChoices }>()
  let argAsk = 0
  const [search, setSearch] = createSignal<Search | null>(null)
  const pastes = new Pastes()
  let walk = new HistoryWalk(deps.history)
  // The text the composer put in itself (a recalled prompt, a completion).
  // The textarea reports changes a moment after they happen, so a change
  // is the composer's own when it leaves exactly this text; anything else
  // is the user's, and ends the history walk and refreshes completion.
  let quiet: string | null = null
  let dismissedAt = -1             // esc closed the list for the token starting here
  let commands: Command[] | null = null
  let commandsAt = 0
  let commandsFor = ""
  let pathAsk = 0

  // Wrapped lines of the whole prompt (virtualLineCount counts only the
  // rows in view), and the cursor's row among them.
  const alive = () => Boolean(area && !area.isDestroyed)
  const totalRows = () => (alive() ? Math.max(1, area!.editorView.getTotalVirtualLineCount() || 1) : 1)
  const cursorRow = () => (alive() ? area!.scrollY + area!.visualCursor.visualRow : 0)

  const measure = () => {
    if (!alive()) return                   // a remeasure after the window closed
    setRows(Math.min(COMPOSER_MAX_ROWS, totalRows()))
  }

  /** The line as it is sent: the "!" of shell mode in front. */
  const full = () => (shell() ? SHELL_PREFIX : "") + (area?.plainText ?? "")

  /** Put `value` (a whole line, "!" and all) in the textarea, the cursor
   *  at `cursor` of it (the end when not given): a leading "!" becomes
   *  shell mode and leaves the text. */
  const place = (value: string, cursor?: number) => {
    if (!area) return
    const lead = shellAllowed() ? shellLead(value) : -1
    const inner = lead >= 0 ? value.slice(lead) : value
    setShell(lead >= 0)
    quiet = inner
    area.setText(inner)
    if (cursor === undefined) area.gotoBufferEnd()
    else area.cursorOffset = Math.max(0, Math.min(inner.length, cursor - Math.max(0, lead)))
    setTextSignal(full())
    measure()
  }
  const setAreaText = (value: string) => place(value)

  /** Out of shell mode: the "!" is gone, the rest stays as a prompt (one
   *  that begins with "!" again is a shell line again). */
  const leaveShell = () => {
    if (!area) return
    setShell(false)
    place(area.plainText, area.cursorOffset)
    if (!shell()) refreshCompletion()
  }

  // The prompt became something else's (the agent view opened): the "!"
  // goes back into the text; back to a prompt, a leading "!" is shell
  // mode again.
  createEffect(on(shellAllowed, () => {
    if (!area || area.isDestroyed) return
    const cursor = area.cursorOffset + (shell() ? SHELL_PREFIX.length : 0)
    place(full(), cursor)
  }, { defer: true }))

  const loadCommands = async (): Promise<Command[]> => {
    const sid = deps.sessionId()
    if (commands && commandsFor === sid && Date.now() - commandsAt < 60_000) return commands
    try {
      const got = await deps.request("session.commands", { sessionId: sid })
      const theirs = (got?.commands ?? []) as Command[]
      // The window's own (/rename) win over the session's of the same name,
      // as the REPL's builtins win over a user's: what enter runs is listed.
      const mine = new Set(WINDOW_COMMANDS.map((c) => c.name))
      commands = [...(WINDOW_COMMANDS as Command[]), ...theirs.filter((c) => !mine.has(c.name))]
      commandsAt = Date.now()
      commandsFor = sid
    } catch {
      commands = commands ?? (WINDOW_COMMANDS as Command[])
    }
    return commands
  }

  /** What may follow `command`, if it takes anything: asked again after
   *  10 s (memories, staged diffs and jobs change while a session runs). */
  const loadArguments = async (command: string): Promise<ArgumentChoices | null> => {
    if (!deps.argumentChoices) return null
    const list = await loadCommands()
    const known = list.find((c) => c.name.toLowerCase() === command.toLowerCase())
    if (!known?.takesArgs) return null
    const key = `${deps.sessionId()}\t${known.name}`
    const hit = argsCache.get(key)
    if (hit && Date.now() - hit.at < 10_000) return hit.got
    try {
      const got = await deps.argumentChoices(known.name)
      argsCache.set(key, { at: Date.now(), got })
      return got
    } catch {
      return null                      // a service that cannot say: no list
    }
  }

  const refreshCompletion = () => {
    if (!area) return
    if (shell()) {
      // A shell line has no /commands or @paths of DELFIN's.
      setPopup(null)
      return
    }
    const value = area.plainText
    const token = tokenAt(value, area.cursorOffset) ?? argumentAt(value, area.cursorOffset)
    if (token?.kind === "argument") {
      if (token.start === dismissedAt) {
        setPopup(null)
        return
      }
      const ask = ++argAsk
      void loadArguments(token.command!).then((got) => {
        if (ask !== argAsk || !area) return
        const now = argumentAt(area.plainText, area.cursorOffset)
        if (!got || !now || now.command !== token.command || now.text !== token.text) {
          if (!got) setPopup(null)
          return
        }
        const candidates = matchArguments(got, now.text)
        // Only exactly what is typed: no list.
        const useful = candidates.some((c) => c.hint || c.insert !== now.text)
        setPopup(useful ? { token: now, candidates, selected: 0 } : null)
      })
      return
    }
    if (!token || token.start === dismissedAt) {
      if (!token) dismissedAt = -1
      setPopup(null)
      return
    }
    if (token.kind === "command") {
      void loadCommands().then((list) => {
        const now = area ? tokenAt(area.plainText, area.cursorOffset) : null
        if (!now || now.kind !== "command" || now.text !== token.text) return
        const candidates = matchCommands(list, token.text)
        setPopup(candidates.length ? { token: now, candidates, selected: 0 } : null)
      })
      return
    }
    const ask = ++pathAsk
    void deps.request("session.completePath", { sessionId: deps.sessionId(), fragment: token.text.slice(1) })
      .then((got) => {
        if (ask !== pathAsk || !area) return
        const now = tokenAt(area.plainText, area.cursorOffset)
        if (!now || now.kind !== "path" || now.text !== token.text) return
        const candidates = pathCandidates((got?.paths ?? []) as string[])
        // Nothing to offer, or only exactly what is typed: no list.
        const useful = candidates.filter((c) => c.insert !== token.text)
        setPopup(useful.length ? { token: now, candidates, selected: 0 } : null)
      })
      .catch(() => setPopup(null))
  }

  const onContent = () => {
    if (!area) return
    const value = area.plainText
    if (!shell() && shellAllowed() && shellLead(value) >= 0 && !search()) {
      // A "!" in front by another way than the key (a paste, an edit at
      // the start): shell mode, as the session will read the line.
      place(value, area.cursorOffset)
      setPopup(null)
      return
    }
    setTextSignal(full())
    measure()
    if (value === quiet) return
    quiet = null
    if (search()) return
    // An edit of a recalled prompt does not end the walk: ↓ still brings
    // back the draft that was being typed (as readline does).
    refreshCompletion()
  }

  const take = (candidate: Candidate) => {
    const p = popup()
    if (!area || !p) return
    // The list may be older than what is typed now (its answer was on its
    // way): it applies only to the token at the cursor, and only while the
    // candidate still extends what is typed.
    const now = p.token.kind === "argument" ? argumentAt(area.plainText, area.cursorOffset)
      : tokenAt(area.plainText, area.cursorOffset)
    // An argument's value is matched by its words, not as a prefix.
    const extends_ = p.token.kind === "argument" || candidate.insert.startsWith(now?.text ?? "")
    if (!now || now.kind !== p.token.kind || now.start !== p.token.start || !extends_) {
      setPopup(null)
      refreshCompletion()
      return
    }
    const done = applyCandidate(area.plainText, now, candidate.insert)
    quiet = done.text
    area.setText(done.text)
    area.cursorOffset = done.cursor
    setTextSignal(done.text)
    measure()
    setPopup(null)
    dismissedAt = -1
    // A folder: the list goes on inside it; a command that takes an
    // argument: the list of what it takes.
    if (candidate.descends || (p.token.kind === "command" && candidate.insert.endsWith(" "))) refreshCompletion()
  }

  /** Enter (or a click) on a candidate: taken, and a value or a command
   *  without arguments runs. A usage line is not taken. */
  const choose = (chosen: Candidate): boolean => {
    if (chosen.hint) return false
    take(chosen)
    if (chosen.runs) submit()
    return true
  }

  const clear = () => {
    pastes.clear()
    walk.reset()
    setPopup(null)
    setSearch(null)
    dismissedAt = -1
    setAreaText("")
  }

  const submit = (how?: "steer"): boolean => {
    if (!area) return false
    // Stripped as the REPL strips a typed line (the session does too).
    const raw = full().trim()
    if (!raw) return false
    // A "!" with nothing after it: no command to run, and no prompt either.
    if (shell() && shellLead(raw) === raw.length) return false
    if (shell() && how === "steer") {
      // Steered, the model would read the line as text: it runs after the turn.
      deps.say("a shell line is not steered: it runs when the turn ends")
      how = undefined
    }
    const expanded = pastes.expand(raw)
    const used = pastes.used(raw)
    deps.history.add({ text: raw, pastes: used })
    walk = new HistoryWalk(deps.history)
    clear()
    void Promise.resolve(deps.send(expanded, how, expanded === raw ? undefined : raw)).then((sent) => {
      // Not sent (the service is away): the prompt is not lost.
      if (sent === false && area && !area.plainText) {
        pastes.restore(used)
        setAreaText(raw)
      }
    })
    return true
  }

  const recall = (entry: Entry | null) => {
    if (!entry) return
    pastes.restore(entry.pastes)
    setAreaText(entry.text)
    setPopup(null)
  }

  const current = (): Entry => {
    const raw = full()
    return { text: raw, pastes: pastes.used(raw) }
  }

  // -- ctrl+r ------------------------------------------------------------------
  const showMatch = (s: Search) => {
    const entry = s.at >= 0 ? deps.history.entries[s.at] : null
    setAreaText(entry ? entry.text : s.draft.text)
  }

  const searchKey = (s: Search, key: KeyEvent): boolean => {
    const name = key.name
    if (name === "escape" || (key.ctrl && (name === "c" || name === "g"))) {
      setSearch(null)
      pastes.restore(s.draft.pastes)
      setAreaText(s.draft.text)
      return true
    }
    if (key.ctrl && name === "r") {
      const older = deps.history.search(s.query, s.at >= 0 ? s.at : deps.history.entries.length)
      if (older >= 0) {
        const next = { ...s, at: older }
        setSearch(next)
        showMatch(next)
      } else deps.say("no older match")
      return true
    }
    if (name === "return" || name === "tab" || name === "right" || name === "left") {
      setSearch(null)
      const entry = s.at >= 0 ? deps.history.entries[s.at] : null
      if (entry) recall(entry)
      else {
        pastes.restore(s.draft.pastes)
        setAreaText(s.draft.text)
      }
      return true
    }
    let query = s.query
    if (name === "backspace") query = query.slice(0, -1)
    else if (isPlain(key) && key.sequence && key.sequence.length >= 1 && !/[\x00-\x1f\x7f]/.test(key.sequence)) {
      query += key.sequence
    } else if (name === "space") query += " "
    else return true                    // nothing else edits the search
    const at = query ? deps.history.search(query) : -1
    const next = { ...s, query, at }
    setSearch(next)
    showMatch(next)
    return true
  }

  // -- keys ----------------------------------------------------------------------
  const handleKey = (key: KeyEvent): boolean => {
    if (!area) return false
    const s = search()
    if (s) return searchKey(s, key)
    if (shell() && isPlain(key)) {
      // Out of shell mode: esc, or backspace on the "!" (the start of the line).
      if (key.name === "escape" || (key.name === "backspace" && area.cursorOffset === 0)) {
        leaveShell()
        return true
      }
    } else if (!shell() && shellAllowed() && isPlain(key) && key.sequence === SHELL_PREFIX
               && area.cursorOffset === 0) {
      // "!" at the start: the line becomes a shell line, the "!" its marker.
      setShell(true)
      setPopup(null)
      setTextSignal(full())
      return true
    }
    const km = deps.keymap
    const p = popup()
    if (p) {
      if (key.name === "up" || key.name === "down") {
        const n = p.candidates.length
        const selected = (p.selected + (key.name === "up" ? n - 1 : 1)) % n
        setPopup({ ...p, selected })
        return true
      }
      if (key.name === "tab" && isPlain(key)) {
        const chosen = p.candidates[p.selected]!
        if (!chosen.hint) take(chosen)
        return true
      }
      if (key.name === "escape") {
        dismissedAt = p.token.start
        setPopup(null)
        return true
      }
      if (km.matches("submit", key)) {
        const chosen = p.candidates[p.selected]!
        // Already typed out in full, or only the usage shown: Enter sends,
        // as without the list.
        if (!chosen.hint && chosen.insert.trimEnd() !== p.token.text) {
          choose(chosen)
          return true
        }
        setPopup(null)
      }
    }
    if (km.matches("historySearch", key)) {
      const next: Search = { query: "", at: -1, draft: current() }
      setSearch(next)
      setPopup(null)
      return true
    }
    if (km.matches("newline", key)) {
      area.newLine()
      return true
    }
    if (km.matches("submit", key)) {
      const before = area.plainText.slice(0, area.cursorOffset)
      if (before.endsWith("\\") && area.cursorOffset === area.plainText.length) {
        // A trailing backslash continues the prompt on a new line.
        area.deleteCharBackward()
        area.newLine()
        return true
      }
      submit()
      return true
    }
    if (km.matches("historyPrev", key)) {
      if (cursorRow() > 0) return false             // up within the prompt
      recall(walk.older(current()))
      return true
    }
    if (km.matches("historyNext", key)) {
      if (cursorRow() < totalRows() - 1) return false
      if (!walk.walking) return false
      recall(walk.newer())
      return true
    }
    if (key.name === "tab" && isPlain(key) && !key.shift) {
      const token = tokenAt(area.plainText, area.cursorOffset) ?? argumentAt(area.plainText, area.cursorOffset)
      if (token) {
        dismissedAt = -1
        refreshCompletion()
      }
      return true                         // a tab never goes into the prompt
    }
    return false
  }

  const handlePaste = (raw: string): boolean => {
    deps.beforePaste?.()
    if (!area || search()) return false
    const text = normalisePaste(raw)
    if (!shouldCollapse(text)) return false
    area.insertText(pastes.add(text))
    return true
  }

  return {
    text,
    value: () => full(),
    rows,
    popup,
    search,
    shell,
    height: () => rows() + (search() ? 1 : 0) + Math.min(POPUP_ROWS, popup()?.candidates.length ?? 0),
    attach(textarea) {
      area = textarea
      textarea.onContentChange = onContent
      textarea.onCursorChange = () => {
        // The cursor moved off the token the list is for: the list follows.
        const p = popup()
        if (!p || !area) return
        const now = p.token.kind === "argument" ? argumentAt(area.plainText, area.cursorOffset)
          : tokenAt(area.plainText, area.cursorOffset)
        if (!now || now.start !== p.token.start) setPopup(null)
      }
      measure()
    },
    handleKey,
    handlePaste,
    submit,
    clear,
    setText(value, p) {
      pastes.restore(p)
      setAreaText(value)
    },
    insert(value) {
      if (alive()) area!.insertText(value)
    },
    remeasure: measure,
    clickAt(index) {
      const p = popup()
      const chosen = p?.candidates[index]
      if (!p || !chosen) return false
      const pause = deps.clickPauseMs ?? CLICK_PAUSE_MS
      if (pause > 0 && Date.now() - popupShownAt < pause) return false
      return choose(chosen)
    },
    cancelSearch() {
      const s = search()
      if (!s) return
      setSearch(null)
      pastes.restore(s.draft.pastes)
      setAreaText(s.draft.text)
    },
  }
}

export interface ComposerViewProps {
  composer: Composer
  theme: () => Theme
  /** The terminal's width: the list and the search line are cut to it. */
  width?: () => number
  placeholder: () => string
  focused: () => boolean
  /** The paste handler is active only while the prompt has the keys. */
  pasteActive: () => boolean
  /** The prompt marker's colour (the posture's); the accent when not given. */
  mark?: () => ColorInput
}

/** The completion list, the search line and the prompt itself. */
export function ComposerView(props: ComposerViewProps): JSX.Element {
  const t = props.theme
  const c = props.composer
  usePaste((event) => {
    if (!props.pasteActive()) return
    const text = stripAnsiSequences(decodePasteBytes(event.bytes))
    if (c.handlePaste(text)) event.preventDefault()
  })
  const shown = () => {
    const p = c.popup()
    if (!p) return []
    const first = Math.max(0, Math.min(p.selected - POPUP_ROWS + 1, p.candidates.length - POPUP_ROWS))
    return p.candidates.slice(first, first + POPUP_ROWS).map((cand, i) => ({ cand, selected: first + i === p.selected,
                                                                             index: first + i }))
  }
  return (
    <box flexShrink={0} flexDirection="column">
      <For each={shown()}>
        {(row) => {
          // One row each, whatever the width: the footer counts them so.
          const room = Math.max(8, (props.width?.() ?? 100) - 5)
          const label = fit(row.cand.label, Math.max(4, Math.floor(room * 0.5)))
          const detail = row.cand.detail ? fit(`  ${row.cand.detail}`, Math.max(0, room - width(label))) : ""
          // A row is the candidate: a click takes it, as enter does.
          const hit = row.cand.hint ? undefined : clickable(() => "take", () => void c.clickAt(row.index))
          if (row.cand.hint) {
            return <text wrapMode="none" fg={t().faint}>{"    " + fit(`usage: ${row.cand.label}`, room)}</text>
          }
          return (
            <text wrapMode="none" ref={hit?.ref} onMouseDown={hit?.onMouseDown} onMouseUp={hit?.onMouseUp}
                  onMouseDrag={hit?.onMouseDrag}>
              <span style={{ fg: row.selected ? t().accent : t().faint }}>{row.selected ? "  › " : "    "}</span>
              <span style={{ fg: row.selected ? t().bright : t().text, bold: Boolean(row.selected) }}>
                {label}
              </span>
              <span style={{ fg: t().faint }}>{detail}</span>
            </text>
          )
        }}
      </For>
      <Show when={c.search()}>
        {(s: () => Search) => {
          const room = () => Math.max(8, (props.width?.() ?? 100) - 1)
          // The hint goes first when the width runs short; the query stays.
          const hint = () => {
            const full = s().query && s().at < 0 ? "  (no match)" : "  ctrl+r older · enter take · esc back"
            return room() - 18 - width(full) >= Math.min(12, width(s().query) + 1) ? full : ""
          }
          return (
            <text wrapMode="none">
              <span style={{ fg: t().accentSoft }}>{"  search history: "}</span>
              <span style={{ fg: t().bright }}>{fit(s().query, Math.max(4, room() - 18 - width(hint())))}</span>
              <span style={{ fg: t().faint }}>{hint()}</span>
            </text>
          )
        }}
      </Show>
      <box flexDirection="row">
        {/* Shell mode: "!" in the shell's colour -- the mark is the text,
            so it reads the same without colour (NO_COLOR). */}
        <text flexShrink={0} width={2} fg={c.shell() ? t().shell : props.mark?.() ?? t().accent}
              attributes={TextAttributes.BOLD}>{c.shell() ? SHELL_PREFIX + " " : "❯ "}</text>
        <textarea
          flexBasis={0}
          flexShrink={1}
          minWidth={4}
          ref={(el: TextareaRenderable) => c.attach(el)}
          flexGrow={1}
          height={c.rows()}
          wrapMode="word"
          focused={props.focused()}
          placeholder={props.placeholder()}
          placeholderColor={t().faint}
          textColor={t().bright}
          focusedTextColor={t().bright}
          backgroundColor="transparent"
          focusedBackgroundColor="transparent"
          cursorColor={c.shell() ? t().shell : t().accent}
          keyBindings={[{ name: "return", action: "submit" }]}
          onSubmit={() => void c.submit()}
        />
      </box>
    </box>
  )
}
