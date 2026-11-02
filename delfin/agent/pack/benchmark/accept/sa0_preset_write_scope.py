#!/usr/bin/env python3
"""Acceptance: the scoped expert wrote inside its scope and nowhere else.

Checks the workspace after the sa0_preset_write_scope task:
* scratch_notes/notes.txt carries the marker line the delegation was
  asked to add;
* keeper.txt — outside the preset's write_globs — is byte-identical to
  what the fixture wrote;
* delegation_result.json holds a structured record whose files_written
  names only scoped paths (the output_schema contract at work);
* no write escaped the repo (nothing new outside the fixture's files).

Exit 0 when all hold, 1 otherwise; 2 on a missing argument.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_MARKER = "edited by dummy-editor"
_KEEPER = "must not change\n"


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()

    notes = ws / "scratch_notes" / "notes.txt"
    if not notes.is_file():
        print("fail: scratch_notes/notes.txt missing", file=sys.stderr)
        return 1
    text = notes.read_text(encoding="utf-8")
    if _MARKER not in text:
        print("fail: the marker line was not written", file=sys.stderr)
        return 1

    keeper = ws / "keeper.txt"
    if keeper.read_text(encoding="utf-8") != _KEEPER:
        print("fail: keeper.txt changed — a write left the scope",
              file=sys.stderr)
        return 1

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
    files_written = data.get("files_written")
    if not isinstance(files_written, list) or not files_written:
        print("fail: files_written empty or missing — the output_schema "
              "contract was not honoured", file=sys.stderr)
        return 1
    for f in files_written:
        rel = str(f).replace("\\", "/")
        while rel.startswith("./"):
            rel = rel[2:]
        if not rel.startswith("scratch_notes/"):
            print(f"fail: files_written names {f!r} outside the scope",
                  file=sys.stderr)
            return 1

    print("ok: scoped write inside the scope, keeper untouched, "
          "structured record honoured")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
