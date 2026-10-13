/**
 * The connection to the agent service across a reconnect: a cut inside a
 * character does not garble the next connection, a hello that never comes
 * does not hang the window, failed attempts leave no socket open, and a
 * service of another protocol version ends the retrying with the reason.
 */
import { expect, test } from "bun:test"
import { mkdtempSync } from "node:fs"
import { tmpdir } from "node:os"
import { join } from "node:path"
import { CONNECTION_LOST, RpcClient, RpcError, ServiceConnection, type Wire } from "../../src/transport/rpc"

const HELLO = '{"protocolVersion":1,"serviceVersion":1,"pid":1,"host":"h","mayAct":true,"sessions":[]}'

test("the first line after a reconnect survives a cut mid-character", async () => {
  const sent: string[] = []
  const wire = (): Wire => ({ write: (l) => sent.push(l), close: () => {} })
  const rpc = new RpcClient(wire())
  const line = new TextEncoder().encode('{"jsonrpc":"2.0","method":"event","params":{"type":"notice","session":"s","text":"zz€"}}\n')
  rpc.feed(line.subarray(0, line.indexOf(0xe2) + 1))          // cut inside "€"
  rpc.lost("the agent service closed the connection")
  rpc.attach(wire())
  const hello = rpc.request("hello", { protocolVersion: 1 })
  const id = JSON.parse(sent.at(-1)!).id
  rpc.feed(new TextEncoder().encode(`{"jsonrpc":"2.0","id":${id},"result":{"ok":true}}\n`))
  expect(await Promise.race([hello.then(() => "answered"), Bun.sleep(300).then(() => "still waiting")])).toBe("answered")
})

test("a request with a time limit fails when no answer comes; a lost one says the connection went", async () => {
  const rpc = new RpcClient({ write: () => {}, close: () => {} })
  const err = await rpc.request("hello", {}, { timeoutMs: 50 }).catch((e) => e) as RpcError
  expect(err.message).toContain("no answer within")
  const waiting = rpc.request("sessions.list", {})
  rpc.lost("gone")
  expect(((await waiting.catch((e) => e)) as RpcError).code).toBe(CONNECTION_LOST)
})

function server(path: string, answer: (n: number) => string, after: (socket: any, n: number) => void = () => {}) {
  let opened = 0
  let open = 0
  const s = Bun.listen({
    unix: path,
    socket: {
      open() {
        opened++
        open++
      },
      close() {
        open--
      },
      data(socket, data) {
        for (const line of Buffer.from(data).toString().split("\n").filter(Boolean)) {
          const msg = JSON.parse(line)
          if (msg.method !== "hello") continue
          socket.write(`{"jsonrpc":"2.0","id":${msg.id},${answer(opened)}}\n`)
          after(socket, opened)
        }
      },
    },
  })
  return { stop: () => s.stop(true), opened: () => opened, open: () => open }
}

for (const where of ["mid-character", "between characters"]) {
  test(`a reconnect after the service died ${where} comes back up`, async () => {
    const path = join(mkdtempSync(join(tmpdir(), "tui-rc-")), "s.sock")
    const svc = server(path, () => `"result":${HELLO}`, (socket, n) => {
      if (n !== 1) return
      const ev = new TextEncoder().encode('{"jsonrpc":"2.0","method":"event","params":{"type":"notice","session":"s","seq":1,"text":"zz€"}}\n')
      socket.write(ev.subarray(0, ev.indexOf(0xe2) + (where === "mid-character" ? 1 : 0)))
      setTimeout(() => socket.end(), 50)
    })
    const conn = new ServiceConnection(path, [50])
    await conn.connect()
    await Bun.sleep(1000)
    svc.stop()
    expect(svc.opened()).toBeGreaterThanOrEqual(2)
    expect(conn.state).toBe("up")
    conn.close()
  })
}

test("a service of another protocol version: no sockets left open, and the reason said", async () => {
  const path = join(mkdtempSync(join(tmpdir(), "tui-rc-")), "s.sock")
  const svc = server(path,
    (n) => n === 1 ? `"result":${HELLO}`
      : `"error":{"code":-32001,"message":"the service speaks protocol 2, this window 1"}`,
    (socket, n) => {
      if (n === 1) setTimeout(() => socket.end(), 30)        // the old service goes away
    })
  const conn = new ServiceConnection(path, [20])
  await conn.connect()
  await Bun.sleep(800)
  const stillOpen = svc.open()
  const attempts = svc.opened() - 1
  svc.stop()
  expect(stillOpen).toBe(0)
  expect(attempts).toBe(1)                                  // it stopped trying
  expect(conn.state).toBe("closed")
  expect(conn.problem).toContain("protocol 2")
})
