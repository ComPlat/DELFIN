/**
 * DELFIN's clusters and what runs in the background: two screens of the
 * window's own, from the agent service's clusters.list and jobs.list
 * (serviceVersion 5 and later).
 *
 *   /clusters   every cluster of remote.clusters: logged in, logged out or
 *               refused, its target, and while logged in the session
 *               budget and the active jobs; the selected one says what is
 *               wrong (a refused one's problems, an error)
 *                 ↑ ↓      select           enter / l   log in (login.tsx)
 *                 o        log out          r           ask again
 *                 space    details          esc         back
 *   /jobs       the agents' background work (a group per session), the
 *               jobs DELFIN's job monitor watches, and the active jobs on
 *               each logged-in cluster; a source that could not be read is
 *               said at the end, never shown as empty. Read-only.
 *                 ↑ ↓      select a job     r           ask again
 *                 space    details          esc         back
 *
 * A screen takes the agent view's place (the whole terminal but the
 * prompt) and every key while it is open, and asks again every 5 s. The
 * service looks at the clusters over ssh (with a cache): an answer may
 * take some seconds, and "asking…" is said until the first one is in.
 *
 * Space shows the selected entry's details under its row (cluster.detail,
 * job.detail: serviceVersion 5): a cluster's session, DELFIN's installs
 * there and its host key; a job's Slurm facts, where it runs and its
 * output's last lines. One entry's at a time, asked again with every poll
 * while they are open -- the last answer stays until the next is in. What
 * does not fit under the row is counted, never drawn over it.
 */
import type { KeyEvent } from "@opentui/core"
import { For, batch, createEffect, onCleanup, type JSX } from "solid-js"
import { createStore } from "solid-js/store"
import { describe, type Keymap } from "../keymap/keymap"
import { oneLine } from "../store/state"
import { fit, width } from "../text"
import type { Theme } from "../theme/theme"

export type Screen = "clusters" | "jobs"

/** The first service that answers clusters.list, cluster.login/logout and jobs.list. */
export const CLUSTERS_SERVICE_VERSION = 5

/** The first that answers cluster.detail and job.detail (space). */
export const DETAILS_SERVICE_VERSION = 6

/** How often an open screen asks the service again. */
export const CLUSTERS_POLL_MS = 5000

/** A cluster as clusters.list has it (delfin/agent/service/clusters_feed.py). */
export interface ClusterRow {
  name: string
  target: string
  state: "refused" | "logged_out" | "logged_in"
  /** Why DELFIN does not use it (a refused one). */
  problems: string[]
  /** What a login asks for: the masked field's fields. */
  requiresPassword: boolean
  requiresOtp: boolean
  /** The session budget while logged in, when the cluster says one. */
  budget: { used: number; limit: number; unit: string } | null
  /** The user's active jobs there while logged in. */
  jobs: number | null
  error: string | null
}

export interface JobRow {
  id: string
  kind: string
  label: string
  state: string
  detail: string
  /** An agent's: shells, watches, agents or wakeups. */
  group: string
}

export interface JobGroup {
  source: "agent" | "monitor" | "cluster"
  /** The session's (agent), "DELFIN's job monitor", the cluster's. */
  title: string
  /** A cluster's name, an agent's session id: what job.detail asks with. */
  cluster: string
  session: string
  rows: JobRow[]
  /** A cluster's: how many of its jobs ended recently. */
  finished: number
}

/** A part of an entry's details (cluster.detail, job.detail): its title,
 *  label and value pairs, and lines (an output tail, an indented list). */
export interface DetailSection {
  title: string
  pairs: [string, string][]
  lines: string[]
}

/** The details open under an entry's row. */
export interface Details {
  /** The entry's: a cluster's name, a job's jobKey. */
  key: string
  /** Its first answer is not in yet. */
  loading: boolean
  sections: DetailSection[]
  error: string
}

/** A field as one plain line, whatever the service sent. */
const clean = (value: unknown, max: number) =>
  oneLine(typeof value === "string" ? value : value == null ? "" : String(value), max)

/** A line shown as it came (an output tail, an indented list): what
 *  `clean` takes out taken out word by word, the spacing kept. */
const cleanLine = (value: unknown, max: number) =>
  (typeof value === "string" ? value : value == null ? "" : String(value))
    .replace(/[^ ]+/g, (word) => clean(word, max)).trimEnd().slice(0, max)

const STATES = new Set(["refused", "logged_out", "logged_in"])

/** The rows clusters.list sent, each field as its type says -- a row
 *  without a name is left out, one of an unknown state is refused. */
export function clusterRows(raw: unknown): ClusterRow[] {
  if (!Array.isArray(raw)) return []
  const out: ClusterRow[] = []
  for (const r of raw as Record<string, any>[]) {
    if (!r || typeof r !== "object" || !clean(r.name, 64)) continue
    const b = r.budget && typeof r.budget === "object" ? r.budget : null
    out.push({
      name: clean(r.name, 64), target: clean(r.target, 120),
      state: STATES.has(r.state) ? r.state : "refused",
      problems: (Array.isArray(r.problems) ? r.problems : []).map((p: unknown) => clean(p, 300)).filter(Boolean),
      requiresPassword: Boolean(r.requiresPassword), requiresOtp: Boolean(r.requiresOtp),
      budget: b && Number.isFinite(Number(b.used)) && Number.isFinite(Number(b.limit))
        ? { used: Number(b.used), limit: Number(b.limit), unit: clean(b.unit, 16) } : null,
      jobs: typeof r.jobs === "number" && Number.isFinite(r.jobs) ? r.jobs : null,
      error: clean(r.error, 300) || null,
    })
  }
  return out
}

const SOURCES = new Set(["agent", "monitor", "cluster"])

/** The groups jobs.list sent, each field as its type says. */
export function jobGroups(raw: unknown): JobGroup[] {
  if (!Array.isArray(raw)) return []
  return (raw as Record<string, any>[]).filter((g) => g && typeof g === "object").map((g) => ({
    source: SOURCES.has(g.source) ? g.source : "agent",
    title: clean(g.title, 120) || clean(g.cluster, 64) || "(untitled)",
    cluster: clean(g.cluster, 64), session: clean(g.session, 64),
    rows: (Array.isArray(g.rows) ? g.rows : []).filter((r: unknown) => r && typeof r === "object")
      .map((r: Record<string, unknown>) => ({ id: clean(r.id, 64), kind: clean(r.kind, 40), label: clean(r.label, 120),
                                              state: clean(r.state, 24), detail: clean(r.detail, 160),
                                              group: clean(r.group, 16) })),
    finished: Math.max(0, Number(g.finished) || 0),
  }))
}

/** A job's identity from one answer to the next: its source, its cluster
 *  or session, an agent's group (a shell's id and a watch's may be the
 *  same) and its id. */
export const jobKey = (g: JobGroup, r: JobRow): string =>
  JSON.stringify([g.source, g.cluster || g.session, r.group, r.id])

/** What job.detail is asked for a row: what its group says it is. */
const jobAsk = (g: JobGroup, r: JobRow): Record<string, unknown> =>
  g.source === "cluster" ? { source: "cluster", cluster: g.cluster, id: r.id }
    : g.source === "agent" ? { source: "agent", session: g.session, group: r.group, id: r.id }
    : { source: "monitor", id: r.id }

/** The sections cluster.detail or job.detail sent, each field as its type
 *  says; a line keeps its spacing. */
export function detailSections(raw: unknown): DetailSection[] {
  if (!Array.isArray(raw)) return []
  return (raw as Record<string, any>[]).filter((s) => s && typeof s === "object").map((s) => ({
    title: clean(s.title, 60),
    pairs: (Array.isArray(s.pairs) ? s.pairs : []).filter(Array.isArray)
      .map((p: unknown[]) => [clean(p[0], 40), clean(p[1], 300)] as [string, string]),
    lines: (Array.isArray(s.lines) ? s.lines : []).map((line: unknown) => cleanLine(line, 300)),
  }))
}

/** A line of the open details: a title, a pair, a line -- or a note
 *  (loading…, what does not fit) or the error. */
export interface DetailLine {
  text: string
  look?: "head" | "note" | "error"
}

/** The open details as lines: each section's title, its pairs (the labels
 *  padded to the section's longest), its lines as they came. */
export function detailLines(details: Details): DetailLine[] {
  if (details.loading) return [{ text: "loading…", look: "note" }]
  const out: DetailLine[] = details.error ? [{ text: details.error, look: "error" }] : []
  for (const s of details.sections) {
    if (s.title) out.push({ text: s.title, look: "head" })
    const pad = Math.max(0, ...s.pairs.map(([label]) => width(label)))
    for (const [label, value] of s.pairs) out.push({ text: `${label}${" ".repeat(pad - width(label))}  ${value}` })
    for (const text of s.lines) out.push({ text })
  }
  return out.length ? out : [{ text: "no details", look: "note" }]
}

/** Why a refused cluster cannot be logged into: its first problem -- and
 *  ctun's command when that is its host key (the service's words name it). */
export function refusal(row: ClusterRow): string {
  const why = row.problems[0] || row.error || "DELFIN does not use it yet"
  const pin = `ctun -t ${row.name} hostkey --pin`
  return `${row.name} cannot be logged into: ${why}${/host ?key/i.test(why) && !why.includes(pin) ? ` (${pin})` : ""}`
}

const num = (n: number) => Number.isInteger(n) ? String(n) : n.toFixed(1)

export function budgetText(b: ClusterRow["budget"]): string {
  return b ? `${num(b.used)}/${num(b.limit)}${b.unit ? ` ${b.unit}` : ""}` : ""
}

export interface ClustersViewState {
  open: Screen | null
  /** The first answer is not in yet. */
  loading: boolean
  /** Why the service did not answer. */
  error: string
  clusters: ClusterRow[]
  groups: JobGroup[]
  /** The sources of /jobs that could not be read. */
  errors: string[]
  /** The selected entry: a cluster (/clusters), a job of any group (/jobs). */
  selected: number
  /** The entry whose details are open (space), and what the service said. */
  details: Details | null
}

export interface ClustersViewDeps {
  request(method: string, params: Record<string, unknown>): Promise<any>
  keymap: Keymap
  say(text: string): void
  /** Enter on a logged-out cluster: the masked login field for it. */
  login(row: ClusterRow): void
  /** The service's version: details from DETAILS_SERVICE_VERSION on. */
  serviceVersion(): number
  /** How often an open screen asks again (CLUSTERS_POLL_MS; tests shorten it). */
  pollMs?: number
}

/** A line of /jobs: a group's title, a job, a group without jobs, a source not read. */
export interface JobLine {
  kind: "title" | "row" | "empty" | "error"
  text: string
  row?: JobRow
  /** A job's: its place among every group's jobs, and its jobKey. */
  index?: number
  key?: string
}

export interface ClustersView {
  state: ClustersViewState
  open(screen: Screen): void
  close(): void
  /** Ask the service again (r, the poll, a login that went through), and
   *  for the open details. */
  refresh(): void
  /** A key while a screen is open; every key is its own. */
  handleKey(key: KeyEvent): boolean
  selectedRow(): ClusterRow | undefined
  /** The selected entry's index: a cluster's, a job's among every group's
   *  (-1: there is none). */
  selectedAt(): number
  jobLines(): JobLine[]
  /** The keys that go back (esc, ←), as the hints line names them. */
  backKeys(): string
}

/** An entry whose details space shows: its key (a cluster's name, a job's
 *  jobKey) and the request that asks for them. */
interface Entry {
  key: string
  method: string
  params: Record<string, unknown>
}

export function createClustersView(deps: ClustersViewDeps): ClustersView {
  const [state, setState] = createStore<ClustersViewState>({
    open: null, loading: false, error: "", clusters: [], groups: [], errors: [], selected: 0, details: null,
  })
  /** The screens with a request on its way: one at a time each. One
   *  asked for meanwhile (r, a login that went through) goes after it; the
   *  poll's waits for its next turn. */
  const asking = new Set<Screen>()
  const again = new Set<Screen>()

  const open = (screen: Screen) => {
    if (state.open === screen) return refresh()
    setState({ open: screen, loading: true, error: "", clusters: [], groups: [], errors: [], selected: 0, details: null })
  }
  const close = () => setState({ open: null, details: null })

  /** The open screen's entries, in the order shown. */
  const entries = (): Entry[] => state.open === "jobs"
    ? state.groups.flatMap((g) => g.rows.map((r) => ({ key: jobKey(g, r), method: "job.detail", params: jobAsk(g, r) })))
    : state.clusters.map((row) => ({ key: row.name, method: "cluster.detail", params: { name: row.name } }))
  const count = () => state.open === "jobs" ? state.groups.reduce((n, g) => n + g.rows.length, 0) : state.clusters.length
  const selectedAt = () => Math.min(state.selected, count() - 1)

  const ask = (queue: boolean) => {
    const screen = state.open
    if (!screen) return
    if (asking.has(screen)) {
      if (queue) again.add(screen)
      return
    }
    asking.add(screen)
    void deps.request(screen === "clusters" ? "clusters.list" : "jobs.list", {}).then((got) => {
      if (state.open !== screen) return
      const was = entries()[selectedAt()]?.key
      batch(() => {
        if (screen === "clusters") setState({ clusters: clusterRows(got?.clusters), loading: false, error: "" })
        else {
          setState({ groups: jobGroups(got?.groups), loading: false, error: "",
                     errors: (Array.isArray(got?.errors) ? got.errors : []).map((e: unknown) => clean(e, 300)).filter(Boolean) })
        }
        // The selection stays on its entry (else as near to where it was
        // as the list allows); details whose entry is gone close.
        const keys = entries().map((e) => e.key)
        const at = was === undefined ? -1 : keys.indexOf(was)
        setState({ selected: at >= 0 ? at : Math.max(0, Math.min(state.selected, keys.length - 1)) })
        if (state.details && !keys.includes(state.details.key)) setState({ details: null })
      })
    }, (error) => {
      if (state.open === screen) setState({ loading: false, error: clean((error as Error)?.message ?? error, 300) })
    }).finally(() => {
      asking.delete(screen)
      if (again.delete(screen) && state.open === screen) refresh()
    })
  }

  /** The open details' request on its way (its entry's key): one at a
   *  time. The last answer stays on the screen until the next is in; one
   *  for an entry no longer open is let go. */
  let detailsAsking: { key: string } | null = null
  const askDetails = () => {
    const key = state.details?.key
    const entry = key === undefined ? undefined : entries().find((e) => e.key === key)
    if (!entry || detailsAsking?.key === entry.key) return
    const mine = { key: entry.key }
    detailsAsking = mine
    const answer = (sections: DetailSection[], error: string) => {
      if (state.details?.key === entry.key) setState("details", { loading: false, sections, error })
    }
    void deps.request(entry.method, entry.params).then(
      (got) => answer(detailSections(got?.sections), clean(got?.error, 300)),
      (error) => answer([], clean((error as Error)?.message ?? error, 300)),
    ).finally(() => {
      if (detailsAsking === mine) detailsAsking = null
    })
  }
  /** Space: the selected entry's details under its row, asked for now --
   *  or closed, when they are its already. */
  const toggleDetails = () => {
    if (deps.serviceVersion() < DETAILS_SERVICE_VERSION) return deps.say("restart the agent service to see details")
    const entry = entries()[selectedAt()]
    if (!entry) return
    if (state.details?.key === entry.key) return setState({ details: null })
    setState({ details: { key: entry.key, loading: true, sections: [], error: "" } })
    askDetails()
  }

  const refresh = () => {
    ask(true)
    askDetails()
  }
  // While a screen is open: asked now, then every 5 s (and its open
  // details with it).
  createEffect(() => {
    if (!state.open) return
    ask(false)
    const timer = setInterval(() => {
      ask(false)
      askDetails()
    }, deps.pollMs ?? CLUSTERS_POLL_MS)
    ;(timer as { unref?: () => void }).unref?.()
    onCleanup(() => clearInterval(timer))
  })

  const selectedRow = () => state.clusters[Math.min(state.selected, state.clusters.length - 1)]

  const jobLines = (): JobLine[] => {
    const out: JobLine[] = []
    let index = 0
    for (const g of state.groups) {
      const finished = g.source === "cluster" && g.finished > 0 ? ` · ${g.finished} finished recently` : ""
      out.push({ kind: "title", text: g.title + finished })
      for (const row of g.rows) out.push({ kind: "row", text: row.label, row, index: index++, key: jobKey(g, row) })
      if (!g.rows.length) out.push({ kind: "empty", text: g.source === "cluster" ? "no active jobs" : "nothing running" })
    }
    if (!state.groups.length && !state.errors.length) out.push({ kind: "empty", text: "nothing runs in the background" })
    if (state.errors.length) out.push({ kind: "title", text: "could not be read" })
    for (const e of state.errors) out.push({ kind: "error", text: e })
    return out
  }

  const logout = async (row: ClusterRow) => {
    deps.say(`logging out of ${row.name}…`)
    try {
      const got = await deps.request("cluster.logout", { name: row.name })
      deps.say(!got?.ok ? `not logged out of ${row.name}: ${got?.error || got?.code || "the service said no"}`
        : got.closed === false ? `${row.name} was not logged in` : `logged out of ${row.name}`)
    } catch (error) {
      deps.say(`not logged out of ${row.name}: ${clean((error as Error)?.message ?? error, 300)}`)
    }
    refresh()
  }

  const handleKey = (key: KeyEvent): boolean => {
    if (!state.open) return false
    const plain = !key.ctrl && !key.meta
    const letter = plain ? (key.sequence ?? "").toLowerCase() : ""
    if (key.name === "escape" || deps.keymap.matches("clustersBack", key)) {
      close()
      return true
    }
    if (key.name === "up" || key.name === "down") {
      // Entry by entry, around at the ends: /jobs steps over the groups'
      // titles and what could not be read.
      const step = key.name === "up" ? -1 : 1
      const n = count()
      if (n) setState({ selected: (Math.min(state.selected, n - 1) + step + n) % n })
      return true
    }
    if (key.name === "space" && plain) {
      toggleDetails()
      return true
    }
    if (letter === "r") {
      refresh()
      return true
    }
    const row = state.open === "clusters" ? selectedRow() : undefined
    if (row && (deps.keymap.matches("submit", key) || letter === "l")) {
      if (row.state === "logged_out") deps.login(row)
      else if (row.state === "logged_in") deps.say(`already logged in to ${row.name}: o logs out`)
      else deps.say(refusal(row))
      return true
    }
    if (row && letter === "o") {
      if (row.state === "logged_in") void logout(row)
      else deps.say(`${row.name} is not logged in`)
      return true
    }
    return true                           // nothing else leaves the screen
  }

  const backKeys = () => deps.keymap.keys("clustersBack").map(describe).join("/") || "esc"

  return { state, open, close, refresh, handleKey, selectedRow, selectedAt, jobLines, backKeys }
}

/** The hints line while a screen is open: what its keys do now. */
export function clustersHelp(view: ClustersView, room = Infinity): string {
  const row = view.state.open === "clusters" ? view.selectedRow() : undefined
  const back = `${view.backKeys()} back`
  const parts = view.state.open === "jobs" ? ["↑↓ select", "space details", "r refresh", back]
    : ["↑↓ select", row?.state === "logged_out" ? "enter log in" : row?.state === "logged_in" ? "o log out" : "",
       "space details", "r refresh", back].filter(Boolean)
  const kept: string[] = []
  for (const part of parts) {
    if (width([...kept, part].join(" · ")) > room) break
    kept.push(part)
  }
  return kept.join(" · ")
}

const LOOK: Record<ClusterRow["state"], { glyph: string; word: string }> = {
  logged_in: { glyph: "●", word: "logged in" },
  logged_out: { glyph: "○", word: "logged out" },
  refused: { glyph: "✗", word: "refused" },
}

/** What the selected cluster shows under its row: a refused one's
 *  problems, another's error. */
function detailOf(row: ClusterRow): string[] {
  if (row.state === "refused") return row.problems.length ? row.problems : [refusal(row)]
  return row.error ? [row.error] : []
}

interface Line {
  kind: "title" | "row" | "detail" | "empty" | "error" | "job"
  text: string
  row?: ClusterRow
  job?: JobRow
  index?: number
  look?: DetailLine["look"]
}

export function ClustersPanel(props: { view: ClustersView; theme: () => Theme; width: () => number;
                                       height: () => number }): JSX.Element {
  const t = props.theme
  const v = props.view
  const inner = () => Math.max(20, props.width() - 2)
  const room = () => Math.max(1, props.height() - 1)
  const lines = (): Line[] => {
    const out: Line[] = []
    const add = (kind: Line["kind"], text: string, look?: Line["look"]) => out.push({ kind, text, look })
    /** Under an entry's row: what the selected one says (`said`), then its
     *  details while they are open -- as many lines as fit under the row
     *  at the top of the screen, the last of them saying how many more. */
    const under = (key: string, said: string[]) => {
      const details = v.state.details?.key === key ? detailLines(v.state.details) : []
      const shown: DetailLine[] = [...said.map((text) => ({ text })), ...details]
      const space = Math.max(1, room() - 1)
      const kept = shown.length > space ? space - 1 : shown.length
      for (const l of shown.slice(0, kept)) add("detail", l.text, l.look)
      const more = shown.length - kept
      if (more) add("detail", `… ${more} more line${more === 1 ? "" : "s"}`, "note")
    }
    if (v.state.loading) add("empty", "asking the agent service… (it looks at each cluster over ssh)")
    else if (v.state.error) add("error", `the agent service did not answer: ${v.state.error}`)
    const at = v.selectedAt()
    if (v.state.open === "jobs") {
      if (!v.state.loading && !v.state.error) {
        for (const l of v.jobLines()) {
          if (l.kind !== "row") add(l.kind, l.text)
          else {
            out.push({ kind: "job", text: l.text, job: l.row, index: l.index })
            under(l.key!, [])
          }
        }
      }
    } else {
      if (!v.state.loading && !v.state.error && !v.state.clusters.length) {
        add("empty", "no clusters: remote.clusters in DELFIN's settings names none")
      }
      v.state.clusters.forEach((row, index) => {
        out.push({ kind: "row", text: row.name, row, index })
        // The short lines step aside while the full details are open (a
        // refused cluster's would say its problems twice).
        under(row.name, index === at && v.state.details?.key !== row.name ? detailOf(row) : [])
      })
    }
    // The selection and what is under it in view -- the row itself always.
    const sel = out.findIndex((l) => l.index === at)
    let end = sel + 1
    while (out[end]?.kind === "detail") end++
    let first = Math.max(0, Math.min(sel - Math.floor(room() / 3), out.length - room()))
    if (end - first > room()) first = Math.min(sel, end - room())
    return out.slice(Math.max(0, first), Math.max(0, first) + room())
  }
  const summary = () => {
    if (v.state.open === "jobs") {
      const n = v.state.groups.reduce((sum, g) => sum + g.rows.length, 0)
      return v.state.loading ? "" : [`${n} job${n === 1 ? "" : "s"}`,
                                     v.state.errors.length ? `${v.state.errors.length} not read` : ""].filter(Boolean).join(" · ")
    }
    const list = v.state.clusters
    const ins = list.filter((r) => r.state === "logged_in").length
    return v.state.loading ? "" : `${list.length} cluster${list.length === 1 ? "" : "s"}${ins ? ` · ${ins} logged in` : ""}`
  }
  const nameWidth = () => Math.min(18, Math.max(8, ...v.state.clusters.map((r) => width(r.name))))
  const tone = (state: ClusterRow["state"]) => state === "logged_in" ? t().ok : state === "refused" ? t().err : t().muted
  const stateTone = (state: string) => /fail|error|cancel|timeout|memory|kill/i.test(state) ? t().err
    : /run|complet/i.test(state) ? t().ok : /pend|queue|wait/i.test(state) ? t().warn : t().muted
  const detailTone = (look: Line["look"]) => look === "head" ? t().text : look === "note" ? t().faint
    : look === "error" ? t().err : t().muted
  return (
    <box flexDirection="column" height={props.height()} flexShrink={0}>
      <text wrapMode="none">
        <span style={{ fg: t().accent, bold: true }}>{v.state.open === "jobs" ? "◆ JOBS  " : "◆ CLUSTERS  "}</span>
        <span style={{ fg: t().muted }}>{summary()}</span>
      </text>
      <For each={lines()}>
        {(line) => {
          if (line.kind === "title") return <text fg={t().accentSoft} wrapMode="none">{fit(` ${line.text}`, inner())}</text>
          if (line.kind === "empty") return <text fg={t().faint} wrapMode="none">{fit(`   ${line.text}`, inner())}</text>
          if (line.kind === "error") return <text fg={t().err} wrapMode="none">{fit(`   ✗ ${line.text}`, inner())}</text>
          if (line.kind === "detail") {
            return (
              <text wrapMode="none">
                <span style={{ fg: t().muted }}>{"      │ "}</span>
                <span style={{ fg: detailTone(line.look), bold: line.look === "head" }}>{fit(line.text, inner() - 8)}</span>
              </text>
            )
          }
          const on = line.index === v.selectedAt()
          if (line.kind === "job") {
            const job = line.job!
            const kind = fit(job.kind, 10).padEnd(10)
            const room = Math.max(10, Math.floor(inner() * 0.36))
            const label = fit(job.label || job.id, room)
            const state = fit(job.state, 12)
            const detail = fit(job.detail, Math.max(0, inner() - 3 - 10 - 2 - room - 2 - 12 - 2))
            return (
              <text wrapMode="none">
                <span style={{ fg: on ? t().accent : t().faint }}>{on ? " › " : "   "}</span>
                <span style={{ fg: t().muted }}>{kind + "  "}</span>
                <span style={{ fg: on ? t().bright : t().text, bold: on }}>{label + " ".repeat(Math.max(0, room - width(label)))}</span>
                <span style={{ fg: stateTone(job.state) }}>{"  " + state.padEnd(12)}</span>
                <span style={{ fg: t().faint }}>{detail ? "  " + detail : ""}</span>
              </text>
            )
          }
          const row = line.row!
          const look = LOOK[row.state]
          const right = [budgetText(row.budget), row.jobs === null ? "" : `${row.jobs} job${row.jobs === 1 ? "" : "s"}`]
            .filter(Boolean).join("  ")
          const target = fit(row.target, Math.max(0, inner() - 5 - nameWidth() - 2 - 10 - 2 - (right ? width(right) + 2 : 0)))
          return (
            <text wrapMode="none">
              <span style={{ fg: on ? t().accent : t().faint }}>{on ? " › " : "   "}</span>
              <span style={{ fg: tone(row.state) }}>{look.glyph + " "}</span>
              <span style={{ fg: on ? t().bright : t().text, bold: on }}>{fit(row.name, nameWidth()).padEnd(nameWidth())}</span>
              <span style={{ fg: tone(row.state) }}>{"  " + look.word.padEnd(10)}</span>
              <span style={{ fg: t().muted }}>{"  " + target}</span>
              <span style={{ fg: t().faint }}>{right ? "  " + right : ""}</span>
            </text>
          )
        }}
      </For>
    </box>
  )
}
