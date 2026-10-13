import { describe, expect, test } from "bun:test"
import { AnswerBlocks, splitMarkdown } from "../../src/scrollback/blocks"

const ANSWER = `## What was wrong

ORCA 6.1 prints the final energy with an extra
column, so the parser picked up the **wrong field**.

- \`delfin/xtb.py\` now reads the value with a regex
- the old layout still parses

| Output | Before | After |
| --- | --- | --- |
| ORCA 6.1.1 | wrong field | ok |

\`\`\`python
energy = float(match.group(1))
\`\`\`

Run the tests to check it.
`

function feed(text: string, sizes: number[]): { committed: string[]; tails: string[] } {
  const committed: string[] = []
  const tails: string[] = []
  const answer = new AnswerBlocks((block) => committed.push(block), (tail) => tails.push(tail))
  let i = 0
  let k = 0
  while (i < text.length) {
    const n = sizes[k++ % sizes.length]!
    answer.append(text.slice(i, i + n))
    i += n
  }
  answer.finish()
  return { committed, tails }
}

describe("splitMarkdown", () => {
  test("every block but the last is final while streaming", () => {
    const { blocks, tail } = splitMarkdown("# Title\n\nFirst para.\n\nSecond pa", false)
    expect(blocks).toEqual(["# Title", "First para."])
    expect(tail).toBe("Second pa")
  })

  test("an unclosed code fence stays in the tail until it closes", () => {
    const open = splitMarkdown("Intro.\n\n```ts\nconst x = 1\n", false)
    expect(open.blocks).toEqual(["Intro."])
    expect(open.tail).toContain("```ts")
    const closed = splitMarkdown("Intro.\n\n```ts\nconst x = 1\n```\n\nAfter", false)
    expect(closed.blocks).toEqual(["Intro.", "```ts\nconst x = 1\n```"])
  })

  test("final commits the rest and leaves no tail", () => {
    const { blocks, tail } = splitMarkdown("A.\n\nB.", true)
    expect(blocks).toEqual(["A.", "B."])
    expect(tail).toBe("")
  })

  test("CRLF arrives as LF", () => {
    expect(splitMarkdown("A.\r\n\r\nB.", true).blocks).toEqual(["A.", "B."])
  })

  test("blank or empty input commits nothing", () => {
    expect(splitMarkdown("", true)).toEqual({ blocks: [], tail: "" })
    expect(splitMarkdown("\n\n", true).blocks).toEqual([])
  })
})

describe("AnswerBlocks", () => {
  const whole = splitMarkdown(ANSWER, true).blocks

  test.each([[[1]], [[7]], [[3, 11, 1, 40]], [[1000]]])(
    "streaming in chunks of %p commits exactly the blocks of the whole answer, once each",
    (sizes) => {
      const { committed, tails } = feed(ANSWER, sizes)
      expect(committed).toEqual(whole)
      expect(tails.at(-1)).toBe("")
    },
  )

  test("what is committed plus the tail is always the text so far", () => {
    const seen: string[] = []
    let tail = ""
    const answer = new AnswerBlocks((b) => seen.push(b), (t) => { tail = t })
    for (const ch of ANSWER) {
      answer.append(ch)
      const rebuilt = [...seen, tail].join("")
      // Blocks lose only the blank lines between them; no character of
      // content may go missing or appear twice.
      expect(rebuilt.replace(/\s+/g, "")).toBe(answer.content.replace(/\s+/g, ""))
    }
  })
})
