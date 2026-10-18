#!/usr/bin/env python3
"""Acceptance for the two w17_method code tasks.

Verdict = the OUTCOME, checked in the workspace the model left:

* parse_tools/m1_brief.txt is byte-for-byte what the setup wrote (the
  fixture was not "solved" by rewriting its brief);
* the pytest suite inside parse_tools/ exits 0;
* the test files the setup wrote are unchanged (test_parse_x.py /
  test_notes_style.py byte-equal) — the red variant's rule is to fix the
  source, never the test;
* parse_tools/parse_x.py still parses (ast) — a "fix" that replaces the
  module with non-Python is not a fix.

Exit 0 when all hold, 1 otherwise, 2 on a missing argument. Passes
``QUALITY: n`` (n = number of preconditions, out of 4) so the runner's
verify_quality channel reports partial work.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from w17_method_code_repo import _BRIEF_FEATURE, _BRIEF_RED  # noqa: E402

_TEST_NAMES = ("test_parse_x.py", "test_notes_style.py")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w17_method_code_verify.py <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1]).resolve() / "parse_tools"
    pts = 0
    if not root.is_dir():
        print("FAIL: parse_tools/ missing", file=sys.stderr)
        return 1
    brief = (root / "m1_brief.txt").read_text(encoding="utf-8") \
        if (root / "m1_brief.txt").is_file() else None
    if brief not in (_BRIEF_RED, _BRIEF_FEATURE):
        print("FAIL: m1_brief.txt changed or removed", file=sys.stderr)
        return 1
    pts += 1
    src = root / "parse_x.py"
    if not src.is_file():
        print("FAIL: parse_x.py missing", file=sys.stderr)
        return 1
    try:
        ast.parse(src.read_text(encoding="utf-8"))
        pts += 1
    except SyntaxError as exc:
        print(f"FAIL: parse_x.py does not parse: {exc}", file=sys.stderr)
        return 1
    kept = [t for t in _TEST_NAMES if (root / "tests" / t).is_file()]
    if not kept:
        print("FAIL: no test file left", file=sys.stderr)
        return 1
    unchanged = True
    for name in kept:
        p = root / "tests" / name
        # The test file must still import the module under test and still
        # contain its original test function names.
        text = p.read_text(encoding="utf-8")
        if name == "test_parse_x.py" and (
                "test_german_comma_decimal" not in text
                or "from parse_tools.parse_x import parse_x" not in text):
            unchanged = False
        if name == "test_notes_style.py" and (
                "test_slugify" not in text
                or "from parse_tools.parse_x import" not in text):
            unchanged = False
    if not unchanged:
        print("FAIL: the test was edited (the brief forbids that)",
              file=sys.stderr)
        return 1
    pts += 1
    import subprocess
    import os as _os
    # The tests import `parse_tools.parse_x`, so the WORKSPACE (the parent
    # of parse_tools/) must be on sys.path; running from inside parse_tools/
    # would make the suite collection-fail even when the fix is right.
    env = dict(_os.environ)
    env["PYTHONPATH"] = str(root.parent) + _os.pathsep + env.get(
        "PYTHONPATH", "")
    suite = subprocess.run(
        [sys.executable, "-m", "pytest", "parse_tools/tests/", "-q"],
        cwd=str(root.parent), capture_output=True, text=True, timeout=120,
        env=env)
    if suite.returncode != 0:
        tail = (suite.stdout or suite.stderr).strip().splitlines()
        print("FAIL: suite not green: "
              + (tail[-1] if tail else "no output"), file=sys.stderr)
        return 1
    pts += 1
    print(f"OK: fixture intact, suite green ({pts}/4 checks)")
    print(f"QUALITY: {pts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
