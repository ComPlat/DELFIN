/**
 * Lines of a frame as the tests read them: without the gutter, the column
 * the window keeps free at the left (theme GUTTER). Only the gutter is
 * taken off -- a line with something in it would keep it, and
 * test/render/gutter.test.tsx checks that column 0 stays blank.
 */
import { GUTTER } from "../src/theme/theme"

const PAD = " ".repeat(GUTTER)

export function ungutter(line: string): string {
  return line.startsWith(PAD) ? line.slice(GUTTER) : line
}

export function ungutterFrame(frame: string): string {
  return frame.split("\n").map(ungutter).join("\n")
}
