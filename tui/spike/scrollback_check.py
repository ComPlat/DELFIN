"""Check a terminal's scrollback against what the spike says it wrote.

Every row the spike commits in test mode carries a unique tag ("{t00001}",
"{t00002}", ...), and its report lists, in order, every tag the scrollback
must hold (``expected``). The terminal's history -- however a harness got
it out of the terminal -- passes when its tags are exactly that list:

  * nothing missing (a row deleted instead of pushed into scrollback, the
    failure of OpenTUI issue #1630);
  * nothing twice (a row written again by a footer that moved);
  * nothing out of order.

Optional checks, each for one harness situation:

  --prefill N   N shell lines ("shell line 1".."N") were printed before the
                spike started; they must all still be there, in order, above
                the spike's first row (the spike started mid-screen).
  --replays N   the spike must have replayed exactly N times (a resize).
  --width W     no row may be wider than W cells (rows were re-wrapped).

Tags are found in the history with every row trimmed and the rows joined,
so a tag that a narrow row split over two lines still counts once.

Usage (also importable: ``check(history, report, ...) -> list[str]``):
    python3 -I spike/scrollback_check.py <history.txt> <report.json>
        [--prefill N] [--replays N] [--width W]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata

TAG = re.compile(r"\{t(\d{5})\}")
SHELL_LINE = re.compile(r"^shell line (\d+)$")


def cells(text: str) -> int:
    """Terminal cells ``text`` takes (wide and full-width characters: 2)."""
    width = 0
    for ch in text:
        if unicodedata.combining(ch):
            continue
        width += 2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1
    return width


def tags_in(history: str) -> list[str]:
    joined = "".join(line.strip() for line in history.splitlines())
    return [f"{{t{m.group(1)}}}" for m in TAG.finditer(joined)]


def _first_difference(got: list[str], want: list[str]) -> int:
    for i, (a, b) in enumerate(zip(got, want)):
        if a != b:
            return i
    return min(len(got), len(want))


def check(history: str, report: dict, *, prefill: int = 0,
          replays: int | None = None, width: int = 0) -> list[str]:
    problems: list[str] = []
    for err in report.get("errors", []) or []:
        problems.append(f"spike error: {err}")
    want = list(report.get("expected") or [])
    if not want:
        problems.append("the report lists no expected rows")
    got = tags_in(history)
    if got != want:
        seen = set(got)
        missing = [t for t in want if t not in seen]
        twice = sorted({t for t in got if got.count(t) > 1})
        stray = sorted(set(got) - set(want))
        at = _first_difference(got, want)
        problems.append(
            f"rows differ: {len(got)} found, {len(want)} expected; "
            f"{len(missing)} missing, {len(twice)} twice, {len(stray)} not expected; "
            f"first difference at #{at}: found {got[at:at + 3]} expected {want[at:at + 3]}")
        if missing:
            problems.append(f"missing (first 8): {missing[:8]}")
        if twice:
            problems.append(f"twice (first 8): {twice[:8]}")
        if stray:
            problems.append(f"not expected (first 8): {stray[:8]}")
    if prefill:
        lines = history.splitlines()
        first_tag = next((i for i, line in enumerate(lines) if TAG.search(line)),
                         len(lines))
        numbers = [int(m.group(1)) for line in lines[:first_tag]
                   if (m := SHELL_LINE.match(line.strip()))]
        if numbers != list(range(1, prefill + 1)):
            lost = sorted(set(range(1, prefill + 1)) - set(numbers))
            problems.append(f"shell lines above the spike: {len(numbers)} of {prefill} "
                            f"intact (lost: {lost[:8]})")
    if replays is not None and report.get("replays") != replays:
        problems.append(f"replays: {report.get('replays')} (expected {replays})")
    if width:
        wide = [(i, cells(line.rstrip())) for i, line in enumerate(history.splitlines())
                if cells(line.rstrip()) > width]
        if wide:
            problems.append(f"{len(wide)} rows wider than {width} cells "
                            f"(first: line {wide[0][0]}, {wide[0][1]} cells)")
    return problems


def summary(history: str, report: dict) -> str:
    return (f"{len(tags_in(history))} tagged rows found, "
            f"{len(report.get('expected') or [])} expected, "
            f"replays {report.get('replays')}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("history")
    ap.add_argument("report")
    ap.add_argument("--prefill", type=int, default=0)
    ap.add_argument("--replays", type=int, default=None)
    ap.add_argument("--width", type=int, default=0)
    args = ap.parse_args()
    history = open(args.history, encoding="utf-8", errors="replace").read()
    report = json.load(open(args.report))
    problems = check(history, report, prefill=args.prefill,
                     replays=args.replays, width=args.width)
    print(summary(history, report))
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: every row exactly once, in order")
    return 0


if __name__ == "__main__":
    sys.exit(main())
