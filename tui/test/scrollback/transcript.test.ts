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

  test("going back to the original width before it settles replays nothing", async () => {
    const { transcript } = fakeTerminal()
    const source = fakeSource(100)
    replayOnWidthChange(source, transcript, { settleMs: 20 })
    source.resize(90)
    source.resize(100)
    await sleep(60)
    expect(transcript.replayCount).toBe(0)
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
