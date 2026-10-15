/**
 * R11 (probe): AnswerBlocks commits block i once blocks.length > i. If a
 * later delta merges block i-1 with what followed (a loose list continuing
 * after a blank line), the merged block keeps index i-1 -- "committed" --
 * and the text that merged into it is never committed. Compare what a
 * stream commits, delta by delta, with what the final text splits into.
 */
import { expect, test } from "bun:test"
import { AnswerBlocks, splitMarkdown } from "../../src/scrollback/blocks"

const answers = [
  "- a\n\n- b\n\n- c\n\nDone.",
  "1. one\n\n2. two\n\n3. three\n\nAfter.",
  "Intro.\n\n- x\n\n  more about x\n\n- y\n\nEnd.",
  "* first\n\n* second\n",
  "Para one.\n\nPara two.\n\n```py\nprint(1)\n```\n\nTail.",
]

for (const text of answers) {
  for (const step of [1, 2, 3, 6]) {
    test(`${JSON.stringify(text).slice(0, 30)} in ${step}-char deltas`, () => {
      const committed: string[] = []
      const blocks = new AnswerBlocks((b) => committed.push(b))
      for (let i = 0; i < text.length; i += step) blocks.append(text.slice(i, i + step))
      blocks.finish()
      const final = splitMarkdown(text, true).blocks
      const got = committed.join("\n\n").replace(/\s+/g, " ")
      const want = final.join("\n\n").replace(/\s+/g, " ")
      if (got !== want) console.log(`step ${step}: committed ${JSON.stringify(committed)}\n        final     ${JSON.stringify(final)}`)
      expect(got).toBe(want)
    })
  }
}
