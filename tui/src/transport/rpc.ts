/**
 * JSON-RPC over the agent service's Unix socket: one JSON object per line.
 *
 * `RpcClient` speaks the framing and matches answers to requests; it is
 * transport-agnostic so a test hands it a fake. `ServiceConnection` opens
 * the socket with Bun, says hello, and reconnects when the service goes
 * away and comes back (an upgrade, a restart): requests fail while it is
 * down, events stop, and `onUp` runs again once it is back, so the app can
 * attach to its sessions again.
 */
import type { AnyEvent, Requests, Results } from "../protocol/generated"
import { ERRORS, PROTOCOL_VERSION } from "../protocol/generated"

export class RpcError extends Error {
  constructor(readonly code: number, message: string) {
    super(message)
  }
}

/** What carries lines both ways. */
export interface Wire {
  write(line: string): void
  close(): void
}

type Pending = { resolve: (value: any) => void; reject: (error: Error) => void; method: string }

export class RpcClient {
  private nextId = 0
  private pending = new Map<number, Pending>()
  private buffer = ""
  private readonly decoder = new TextDecoder()
  private readonly listeners = new Set<(event: AnyEvent) => void>()

  constructor(private wire: Wire | null = null) {}

  attach(wire: Wire): void {
    this.wire = wire
  }

  onEvent(listener: (event: AnyEvent) => void): () => void {
    this.listeners.add(listener)
    return () => this.listeners.delete(listener)
  }

  request<M extends keyof Requests & keyof Results>(method: M, params: Requests[M]): Promise<Results[M]>
  request(method: string, params?: Record<string, unknown>): Promise<any>
  request(method: string, params: Record<string, unknown> = {}): Promise<any> {
    if (!this.wire) return Promise.reject(new RpcError(-1, "not connected"))
    const id = ++this.nextId
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject, method })
      try {
        this.wire!.write(JSON.stringify({ jsonrpc: "2.0", id, method, params }) + "\n")
      } catch (error) {
        this.pending.delete(id)
        reject(new RpcError(-1, `could not send ${method}: ${(error as Error).message}`))
      }
    })
  }

  /** Bytes from the wire, in order; whole lines are handled as they complete. */
  feed(chunk: Uint8Array | string): void {
    this.buffer += typeof chunk === "string" ? chunk : this.decoder.decode(chunk, { stream: true })
    let nl: number
    while ((nl = this.buffer.indexOf("\n")) >= 0) {
      const line = this.buffer.slice(0, nl).trim()
      this.buffer = this.buffer.slice(nl + 1)
      if (line) this.handle(line)
    }
  }

  private handle(line: string): void {
    let message: any
    try {
      message = JSON.parse(line)
    } catch {
      return
    }
    if (!message || typeof message !== "object") return
    if (message.method === "event") {
      for (const listener of this.listeners) {
        try {
          listener(message.params as AnyEvent)
        } catch (error) {
          console.error("event listener failed:", error)
        }
      }
      return
    }
    const waiting = this.pending.get(message.id)
    if (!waiting) return
    this.pending.delete(message.id)
    if (message.error) {
      waiting.reject(new RpcError(Number(message.error.code), String(message.error.message)))
    } else {
      waiting.resolve(message.result)
    }
  }

  /** The wire is gone: every request still waiting fails. */
  lost(reason: string): void {
    this.wire = null
    this.buffer = ""
    for (const [id, waiting] of this.pending) {
      waiting.reject(new RpcError(ERRORS.WORKER_GONE, `${waiting.method}: ${reason}`))
      this.pending.delete(id)
    }
  }
}

export type ConnectionState = "connecting" | "up" | "reconnecting" | "closed"

export interface ServiceHello {
  protocolVersion: number
  serviceVersion: number
  pid: number
  host: string
  mayAct: boolean
  sessions: Record<string, unknown>[]
}

export class ServiceConnection {
  readonly rpc = new RpcClient()
  state: ConnectionState = "connecting"
  hello: ServiceHello | null = null
  private closing = false
  private socket: any = null
  private readonly stateListeners = new Set<(state: ConnectionState) => void>()
  private readonly upListeners = new Set<(hello: ServiceHello) => void>()

  constructor(readonly path: string, private readonly retryMs: number[] = [200, 500, 1000, 2000, 5000]) {}

  onState(listener: (state: ConnectionState) => void): () => void {
    this.stateListeners.add(listener)
    return () => this.stateListeners.delete(listener)
  }

  /** Runs after every successful hello: the first connect and each reconnect. */
  onUp(listener: (hello: ServiceHello) => void): () => void {
    this.upListeners.add(listener)
    return () => this.upListeners.delete(listener)
  }

  private setState(state: ConnectionState): void {
    this.state = state
    for (const listener of this.stateListeners) listener(state)
  }

  /** Connect and say hello; resolves with the hello, rejects only when the
   *  first attempt fails (a window with no service has nothing to show). */
  async connect(): Promise<ServiceHello> {
    await this.open()
    return this.greet()
  }

  private async open(): Promise<void> {
    const rpc = this.rpc
    const encoder = new TextEncoder()
    // Bytes the socket did not take yet; written on drain, in order.
    let unsent = new Uint8Array(0)
    const self = this
    const flush = (socket: any) => {
      while (unsent.length) {
        const wrote = socket.write(unsent)
        if (wrote <= 0) return
        unsent = unsent.subarray(wrote)
      }
    }
    this.socket = await Bun.connect({
      unix: this.path,
      socket: {
        data(_socket, data) {
          rpc.feed(data)
        },
        drain(socket) {
          flush(socket)
        },
        close() {
          self.lost("the agent service closed the connection")
        },
        error(_socket, error) {
          self.lost(`the connection failed: ${error.message}`)
        },
      },
    })
    const socket = this.socket
    rpc.attach({
      write(line: string) {
        const bytes = encoder.encode(line)
        if (unsent.length) {
          const joined = new Uint8Array(unsent.length + bytes.length)
          joined.set(unsent)
          joined.set(bytes, unsent.length)
          unsent = joined
          return
        }
        unsent = bytes
        flush(socket)
      },
      close() {
        socket.end()
      },
    })
  }

  private async greet(): Promise<ServiceHello> {
    const hello = await this.rpc.request("hello", { protocolVersion: PROTOCOL_VERSION, client: "delfin-tui" })
    this.hello = hello as ServiceHello
    this.setState("up")
    for (const listener of this.upListeners) listener(this.hello)
    return this.hello
  }

  private lost(reason: string): void {
    // Whatever waits on this connection fails now -- also a hello sent
    // while reconnecting, which would otherwise wait forever.
    this.rpc.lost(reason)
    this.socket = null
    if (this.state === "reconnecting" || this.state === "closed") return
    if (this.closing) {
      this.setState("closed")
      return
    }
    this.setState("reconnecting")
    void this.reconnect()
  }

  private async reconnect(): Promise<void> {
    for (let attempt = 0; !this.closing; attempt++) {
      const wait = this.retryMs[Math.min(attempt, this.retryMs.length - 1)]!
      await Bun.sleep(wait)
      if (this.closing) return
      try {
        await this.open()
        await this.greet()
        return
      } catch {
        this.rpc.lost("still reconnecting")
        this.socket = null
      }
    }
  }

  close(): void {
    this.closing = true
    this.setState("closed")
    try {
      this.socket?.end()
    } catch {
      // already gone
    }
  }
}
