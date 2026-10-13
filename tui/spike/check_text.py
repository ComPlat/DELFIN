"""Check a captured scrollback text against a spike report.

Same rules as check_tmux.py and check_xterm.ts, for text captured any other
way (for example GNOME Terminal's select-all + copy):

  * every marker the spike reported appears, in order;
  * all 1,000 replay lines appear exactly once and in sequence.

Usage: python3 -I spike/check_text.py <scrollback.txt> <report.json>
"""

from __future__ import annotations

import json
import re
import sys


def main() -> int:
    history = open(sys.argv[1], encoding="utf-8", errors="replace").read()
    report = json.load(open(sys.argv[2]))
    problems: list[str] = []
    pos = 0
    for marker in report.get("markers", []):
        at = history.find(marker, pos)
        if at < 0:
            problems.append(f"marker missing or out of order: {marker!r}")
        else:
            pos = at
    numbers = [int(n) for n in re.findall(r"replay line (\d+) —", history)]
    if numbers != list(range(1, 1001)):
        seen = set(numbers)
        missing = sorted(set(range(1, 1001)) - seen)
        problems.append(f"replay lines: {len(numbers)} seen, {len(numbers) - len(seen)}"
                        f" duplicated, {len(missing)} missing (first: {missing[:5]})")
    print(f"{history.count(chr(10))} lines captured")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: markers in order, 1000 replay lines exactly once")
    return 0


if __name__ == "__main__":
    sys.exit(main())
