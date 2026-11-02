#!/usr/bin/env python3
"""Acceptance: the expert could not spawn, and the refusal was reported.

Checks the workspace after the sa0_expert_cannot_spawn task:
* delegation_result.json holds {"nested_refused": true} (or an
  equivalent true-ish boolean) — the refusal REACHED the record;
* keeper.txt and scratch_notes/notes.txt are byte-identical to what the
  fixture wrote: a refused delegation must not have mutated anything;
* no nested sub-agent session was left behind (no .delfin/subagent_
  sessions entry created by a second nesting level is checkable only
  from the record; the file state is the artifact graded).

Exit 0 when all hold, 1 otherwise; 2 on a missing argument.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_KEEPER = "must not change\n"
_NOTES = "original line\n"


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()

    record = ws / "delegation_result.json"
    if not record.is_file():
        print("fail: delegation_result.json missing", file=sys.stderr)
        return 1
    try:
        data = json.loads(record.read_text(encoding="utf-8"))
    except ValueError as exc:
        print(f"fail: delegation_result.json is not JSON: {exc}",
              file=sys.stderr)
        return 1
    refused = data.get("nested_refused")
    if refused is not True:
        print("fail: nested_refused is not true — the expert spawned or "
              "the refusal was swallowed", file=sys.stderr)
        return 1

    keeper = ws / "keeper.txt"
    if keeper.read_text(encoding="utf-8") != _KEEPER:
        print("fail: keeper.txt changed", file=sys.stderr)
        return 1
    notes = ws / "scratch_notes" / "notes.txt"
    if notes.read_text(encoding="utf-8") != _NOTES:
        print("fail: notes.txt changed", file=sys.stderr)
        return 1

    print("ok: nested spawn refused and reported, repo untouched")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
