import { describe, expect, test } from "bun:test"
import { CommitQueue } from "../../src/scrollback/queue"
import { replayOnWidthChange } from "../../src/scrollback/resize"
import { Transcript, type Item } from "../../src/scrollback/transcript"
import { isDrawableDiff } from "../../src/scrollback/writers"

const sleep = (ms: number) => new Promise<void>((r) => setTimeout(r, ms))

/** A fake terminal: what was written, in order, including wipes. */
function fakeTerminal() {
  const screen: string[] = []
  const target = {
    requestLive() {},
    dropLive() {},
    resetSplitFooterForReplay() {
      screen.length = 0
    },
  }
  const write = async (item: Item) => {
    await sleep(Number(item.delay ?? 0))
    screen.push(item.id)
  }
  const queue = new CommitQueue(target)
  return { screen, transcript: new Transcript(queue, target, write), queue }
}

const item = (id: string, delay = 0): Item => ({ id, kind: "row", delay })

describe("Transcript", () => {
  test("items are kept in the order the terminal received them", async () => {
    const { screen, transcript } = fakeTerminal()
    await Promise.all([transcript.append(item("a", 20)), transcript.append(item("b")),
                       transcript.append(item("c", 5))])
    expect(screen).toEqual(["a", "b", "c"])
    expect(transcript.items.map((i) => i.id)).toEqual(["a", "b", "c"])
  })

  test("a replay wipes and rewrites, and what is appended meanwhile is written once, after it", async () => {
    const { screen, transcript } = fakeTerminal()
    transcript.append(item("a"))
    transcript.append(item("b", 10))
    const replay = transcript.replay()
    transcript.append(item("c"))
    await replay
    await transcript.append(item("d"))
    expect(screen).toEqual(["a", "b", "c", "d"])
    expect(transcript.replayCount).toBe(1)
  })

  test("a limited replay writes only the newest items", async () => {
    const { screen, transcript } = fakeTerminal()
    for (const id of ["a", "b", "c", "d"]) transcript.append(item(id))
    await transcript.replay(2)
    expect(screen).toEqual(["c", "d"])
  })
})

describe("replayOnWidthChange", () => {
  function fakeSource(width: number) {
    const listeners = new Set<() => void>()
    return {
      terminalWidth: width,
      terminalHeight: 40,
      on: (_e: "resize", l: () => void) => listeners.add(l),
      off: (_e: "resize", l: () => void) => listeners.delete(l),
      resize(w: number) {
        this.terminalWidth = w
        for (const l of listeners) l()
      },
    }
  }

  test("a footer that changes height (same width) never replays", async () => {
    const { transcript } = fakeTerminal()
    const source = fakeSource(100)
    replayOnWidthChange(source, transcript, { settleMs: 10 })
    source.resize(100)
    source.resize(100)
    await sleep(40)
    expect(transcript.replayCount).toBe(0)
  })

  test("a storm of width changes replays once, after it settles", async () => {
    const { transcript } = fakeTerminal()
    const source = fakeSource(100)
    replayOnWidthChange(source, transcript, { settleMs: 20 })
    for (const w of [99, 95, 90, 85, 80]) {
      source.resize(w)
      await sleep(3)
    }
    await sleep(80)
    expect(transcript.replayCount).toBe(1)
  })

  test("going back to the written width before it settles replays once: the window was drawn at the other", async () => {
    const { transcript } = fakeTerminal()
    const source = fakeSource(100)
    replayOnWidthChange(source, transcript, { settleMs: 20 })
    source.resize(90)
    source.resize(100)
    await sleep(60)
    expect(transcript.replayCount).toBe(1)
  })

  test("a storm that settles once at a width and ends at it writes the session again at the end", async () => {
    // tmux's sequence: the storm settles at 100 (that replay is written while
    // the pane already has another width), goes on, and ends at 100 again.
    const { screen, transcript } = fakeTerminal()
    await transcript.append(item("a"))
    const widths: number[] = []
    const source = fakeSource(110)
    replayOnWidthChange(source, transcript, {
      settleMs: 20,
      onSettled: (w) => {
        widths.push(w)
        void transcript.replay()
      },
    })
    source.resize(100)
    await sleep(60)
    for (const w of [120, 135, 100]) source.resize(w)
    await sleep(60)
    expect(widths).toEqual([100, 100])
    expect(transcript.replayCount).toBe(2)
    expect(screen).toEqual(["a"])
  })

  test("a size the terminal reported but the renderer never applied replays once", async () => {
    // OpenTUI applies only the last size reported while rows are queued,
    // and says nothing when that is the size it has: 100 -> 70 -> 100
    // during a replay reaches only the reports.
    const { transcript } = fakeTerminal()
    const source = fakeSource(100)
    let report: ((w: number, h: number) => void) | undefined
    let stopped = false
    const stop = replayOnWidthChange(source, transcript, {
      settleMs: 20,
      reports: (listener) => {
        report = listener
        return () => { stopped = true }
      },
    })
    report!(100, 40)
    await sleep(60)
    expect(transcript.replayCount).toBe(0)
    report!(70, 40)
    report!(100, 40)
    await sleep(60)
    expect(transcript.replayCount).toBe(1)
    stop()
    expect(stopped).toBe(true)
  })

  test("a footer that changes height while a width change settles neither cancels nor repeats it", async () => {
    const { transcript } = fakeTerminal()
    const source = fakeSource(100)
    replayOnWidthChange(source, transcript, { settleMs: 20 })
    source.resize(90)
    source.resize(90)
    await sleep(60)
    source.resize(90)
    await sleep(60)
    expect(transcript.replayCount).toBe(1)
  })
})

describe("isDrawableDiff", () => {
  test("a diff whose hunk header lies about its lines is not drawable", () => {
    const bad = "--- a/x\n+++ b/x\n@@ -41,7 +41,9 @@\n ctx\n-old\n+new1\n+new2\n ctx2\n"
    const good = "--- a/x\n+++ b/x\n@@ -41,3 +41,4 @@\n ctx\n-old\n+new1\n+new2\n ctx2\n"
    expect(isDrawableDiff(bad)).toBe(false)
    expect(isDrawableDiff(good)).toBe(true)
    expect(isDrawableDiff("not a diff")).toBe(false)
  })
})

test("a session that runs again after it ended says its next end too", async () => {
  const { Presenter } = await import("../../src/scrollback/present")
  const appended: { kind: string; text?: unknown }[] = []
  const transcript = { append: async (item: { kind: string; text?: unknown }) => { appended.push(item) } }
  const p = new Presenter("s1", transcript as never, { showThinking: false })
  const draft = { sessions: { s1: { tail: "", thinking: "" } } } as never
  const ev = (type: string, fields: Record<string, unknown> = {}) => ({ type, session: "s1", ts: 0, ...fields }) as never
  p.handle([ev("session.info", { seq: 1, workspace: "/w" }), ev("session.ended", { seq: 2, reason: "shutdown" })], draft)
  // Woken (a new run's session.info), then its real end.
  p.handle([ev("session.info", { seq: 3, workspace: "/w" }), ev("session.ended", { seq: 4, reason: "exit" })], draft)
  // Woken under a replay that leaves the new run's session.info out (the
  // header stands for it): its first message says it goes on.
  p.handle([ev("user.message", { seq: 6, text: "again", source: "user" }), ev("session.removed", { reason: "stopped" })], draft)
  // A removal after the worker's own last word is not said twice.
  p.handle([ev("session.ended", { seq: 7, reason: "exit" })], draft)
  expect(appended.filter((i) => i.kind === "ended").map((i) => i.text)).toEqual(["session ended (shutdown)",
                                                                              "session ended (exit)",
                                                                              "session ended (stopped)"])
})
