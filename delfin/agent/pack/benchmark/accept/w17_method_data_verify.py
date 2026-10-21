#!/usr/bin/env python3
"""Acceptance for the two w17_method data tasks (shared builder).

Not an acceptance script on its own: w17_method_data_audit_verify.py and
w17_method_data_convert_verify.py call ``check(workspace, variant)``.

Verdict = the artifact, checked in the workspace the model left:

* m1_brief.txt is byte-identical to what the setup wrote (the fixture was
  not "solved" by rewriting the brief);
* the input CSV under data/ is unchanged (the task reads data, never
  edits it);
* variant "audit": data/AUDIT.txt names exactly the duplicate pair
  r1/r2 and no wrong pair;
* variant "convert": data/converted.json maps every id of values.csv to
  its correct float (German decimal comma) and nothing else.

Exit 0 when all hold, 1 otherwise, 2 on a missing argument. Prints
``QUALITY: n`` (n of 4 preconditions) for the verify_quality channel.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from w17_method_data_repo import (  # noqa: E402
    EXPECTED_CONVERT, RUNS_CSV, VALUES_CSV, _BRIEF_AUDIT, _BRIEF_CONVERT)


def main(argv: list[str], variant: str) -> int:
    if len(argv) < 2:
        print(f"usage: <verify> <workspace> ({variant})", file=sys.stderr)
        return 2
    root = Path(argv[1]).resolve() / "data_notes"
    if not root.is_dir():
        print("FAIL: data_notes/ missing", file=sys.stderr)
        return 1
    pts = 0
    brief = root / "m1_brief.txt"
    want = _BRIEF_AUDIT if variant == "audit" else _BRIEF_CONVERT
    if not brief.is_file() or brief.read_text(encoding="utf-8") != want:
        print("FAIL: m1_brief.txt changed or removed", file=sys.stderr)
        return 1
    pts += 1
    csv_path = root / "data" / ("runs.csv" if variant == "audit"
                                else "values.csv")
    if not csv_path.is_file():
        print("FAIL: input csv missing", file=sys.stderr)
        return 1
    csv_now = csv_path.read_text(encoding="utf-8")
    if csv_now != (RUNS_CSV if variant == "audit" else VALUES_CSV):
        print("FAIL: the input csv was edited (the task reads data, "
              "never edits it)", file=sys.stderr)
        return 1
    pts += 1
    if variant == "audit":
        out = root / "data" / "AUDIT.txt"
        if not out.is_file():
            print("FAIL: no data/AUDIT.txt", file=sys.stderr)
            return 1
        text = out.read_text(encoding="utf-8").strip().lower()
        ids = sorted(t for t in text.replace(":", " ").replace("-", " ")
                     .split() if t.startswith("r") and t[1:].isdigit())
        if ids != ["r1", "r2"]:
            print(f"FAIL: AUDIT.txt names {ids}, want exactly r1 r2",
                  file=sys.stderr)
            return 1
        pts += 1
        print("OK: duplicate pair named, fixture intact")
    else:
        out = root / "data" / "converted.json"
        if not out.is_file():
            print("FAIL: no data/converted.json", file=sys.stderr)
            return 1
        try:
            data = json.loads(out.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            print(f"FAIL: converted.json is not JSON: {exc}", file=sys.stderr)
            return 1
        if not isinstance(data, dict) or data != EXPECTED_CONVERT:
            print("FAIL: converted.json does not map every id to its "
                  "correct float", file=sys.stderr)
            return 1
        pts += 1
        print("OK: all values converted, fixture intact")
    pts += 1
    print(f"QUALITY: {pts}")
    return 0


if __name__ == "__main__":
    print("not an acceptance script on its own; see the docstring",
          file=sys.stderr)
    raise SystemExit(2)
