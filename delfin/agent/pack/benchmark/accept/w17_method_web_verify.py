#!/usr/bin/env python3
"""Acceptance for the two w17_method web tasks (shared builder).

Not an acceptance script on its own: w17_method_web_manual_verify.py and
w17_method_web_changelog_verify.py call ``main(argv, variant)``.

Verdict = the artifact:

* m1_brief.txt byte-identical;
* the source document (docs/manual.html / docs/CHANGELOG.md) unchanged;
* variant "manual": docs/answer.txt contains --baseline and does NOT
  contain --bessel (the decoy summary's wrong flag);
* variant "changelog": docs/answer.txt contains "removed" (any case) and
  "strict", and does NOT claim lenient still exists ("lenient" alone is
  not disqualifying — quoting the question is fine — but "still exists"
  or "--lenient" as an answer is).

Exit 0 when all hold, 1 otherwise, 2 on a missing argument.
Prints QUALITY: n of 3.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from w17_method_web_repo import (  # noqa: E402
    _BRIEF_CHANGELOG, _BRIEF_MANUAL, _CHANGELOG_MD, _MANUAL_HTML)


def main(argv: list[str], variant: str) -> int:
    if len(argv) < 2:
        print(f"usage: <verify> <workspace> ({variant})", file=sys.stderr)
        return 2
    root = Path(argv[1]).resolve() / "web_notes"
    if not root.is_dir():
        print("FAIL: web_notes/ missing", file=sys.stderr)
        return 1
    pts = 0
    want = _BRIEF_MANUAL if variant == "manual" else _BRIEF_CHANGELOG
    brief = root / "m1_brief.txt"
    if not brief.is_file() or brief.read_text(encoding="utf-8") != want:
        print("FAIL: m1_brief.txt changed or removed", file=sys.stderr)
        return 1
    pts += 1
    src = root / "docs" / ("manual.html" if variant == "manual"
                           else "CHANGELOG.md")
    if not src.is_file() or src.read_text(encoding="utf-8") != (
            _MANUAL_HTML if variant == "manual" else _CHANGELOG_MD):
        print("FAIL: the source document was changed", file=sys.stderr)
        return 1
    pts += 1
    ans = root / "docs" / "answer.txt"
    if not ans.is_file():
        print("FAIL: no docs/answer.txt", file=sys.stderr)
        return 1
    text = ans.read_text(encoding="utf-8")
    low = text.lower()
    if variant == "manual":
        if "--baseline" not in low:
            print("FAIL: answer.txt does not name --baseline", file=sys.stderr)
            return 1
        if "--bessel" in low:
            print("FAIL: answer.txt repeats the decoy flag --bessel",
                  file=sys.stderr)
            return 1
        print("OK: --baseline named, decoy not repeated")
    else:
        if "removed" not in low or "strict" not in low:
            print("FAIL: answer.txt does not state the removal and the "
                  "strict replacement", file=sys.stderr)
            return 1
        if ("still exists" in low
                or ("--lenient" in low and "removed" not in low)):
            print("FAIL: answer.txt claims lenient still exists", file=sys.stderr)
            return 1
        print("OK: removal stated, replacement named")
    pts += 1
    print(f"QUALITY: {pts}")
    return 0


if __name__ == "__main__":
    print("not an acceptance script on its own; see the docstring",
          file=sys.stderr)
    raise SystemExit(2)
