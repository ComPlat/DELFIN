import { describe, expect, test } from "bun:test"
import { CommitQueue, bounded } from "../../src/scrollback/queue"

function liveCounter() {
  const calls: string[] = []
  return { calls, target: { requestLive: () => calls.push("live"), dropLive: () => calls.push("drop") } }
}

const sleep = (ms: number) => new Promise<void>((r) => setTimeout(r, ms))

describe("CommitQueue", () => {
  test("tasks finish in the order they were queued, however long each takes", async () => {
    const { target } = liveCounter()
    const q = new CommitQueue(target)
    const done: number[] = []
    const waits = [30, 1, 15, 0, 5]
    await Promise.all(waits.map((ms, i) => q.run(async () => {
      await sleep(ms)
      done.push(i)
    })))
    expect(done).toEqual([0, 1, 2, 3, 4])
  })

  test("a failing task is reported and the next one still runs", async () => {
    const { target } = liveCounter()
    const errors: unknown[] = []
    const q = new CommitQueue(target, (e) => errors.push(e))
    const done: string[] = []
    q.run(() => { throw new Error("boom") })
    q.run(async () => { await sleep(1); throw new Error("async boom") })
    await q.run(() => { done.push("after") })
    expect(done).toEqual(["after"])
    expect(errors.map((e) => (e as Error).message)).toEqual(["boom", "async boom"])
  })

  test("the renderer is held live exactly while work is queued", async () => {
    const { calls, target } = liveCounter()
    const q = new CommitQueue(target)
    q.run(() => sleep(5))
    q.run(() => sleep(5))
    expect(calls).toEqual(["live"])
    await q.drain()
    expect(calls).toEqual(["live", "drop"])
    expect(q.size).toBe(0)
    await q.run(() => {})
    expect(calls).toEqual(["live", "drop", "live", "drop"])
  })

  test("an error reporter that throws does not stall the queue", async () => {
    const { target } = liveCounter()
    const q = new CommitQueue(target, () => { throw new Error("reporter") })
    q.run(() => { throw new Error("task") })
    let ran = false
    await q.run(() => { ran = true })
    expect(ran).toBe(true)
  })
})

describe("bounded", () => {
  test("gives up after the bound", async () => {
    const never = new Promise<string>(() => {})
    const started = Date.now()
    expect(await bounded(never, 30)).toBeUndefined()
    expect(Date.now() - started).toBeLessThan(500)
  })

  test("returns the value when it comes first", async () => {
    expect(await bounded(Promise.resolve("ok"), 1000)).toBe("ok")
  })
})
