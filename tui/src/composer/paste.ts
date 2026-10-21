/**
 * Large pastes stand in the composer as one short placeholder,
 * `[Pasted text #1 +42 lines]`, and are put back when the prompt is sent:
 * a 500-line log must not push the prompt out of view, and must still
 * reach the model whole.
 *
 * The placeholder is ordinary text. Deleting it drops the paste; editing
 * it into something else leaves that text as typed (only an exact
 * placeholder is expanded).
 */

/** A paste longer than this many lines, or characters, is collapsed. */
export const COLLAPSE_LINES = 8
export const COLLAPSE_CHARS = 2000

export function normalisePaste(text: string): string {
  return text.replace(/\r\n?/g, "\n")
}

export function shouldCollapse(text: string): boolean {
  return text.length > COLLAPSE_CHARS || text.split("\n").length > COLLAPSE_LINES
}

export class Pastes {
  private next = 1
  private readonly byLabel = new Map<string, string>()

  /** Keep `text`; the placeholder that stands for it. */
  add(text: string): string {
    const lines = text.split("\n").length
    const label = `[Pasted text #${this.next++} +${lines} lines]`
    this.byLabel.set(label, text)
    return label
  }

  /** The prompt with every placeholder still in it put back. */
  expand(text: string): string {
    let out = text
    for (const [label, content] of this.byLabel) {
      const at = out.indexOf(label)
      if (at >= 0) out = out.slice(0, at) + content + out.slice(at + label.length)
    }
    return out
  }

  /** The placeholders `text` still holds, with what they stand for. */
  used(text: string): Record<string, string> {
    const out: Record<string, string> = {}
    for (const [label, content] of this.byLabel) if (text.includes(label)) out[label] = content
    return out
  }

  /** A prompt from the history comes back with its pastes. */
  restore(pastes: Record<string, string> | undefined): void {
    this.clear()
    for (const [label, content] of Object.entries(pastes ?? {})) {
      this.byLabel.set(label, content)
      const n = Number(/#(\d+)/.exec(label)?.[1] ?? 0)
      if (n >= this.next) this.next = n + 1
    }
  }

  clear(): void {
    this.byLabel.clear()
    this.next = 1
  }

  get size(): number {
    return this.byLabel.size
  }
}
