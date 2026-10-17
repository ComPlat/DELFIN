"""Accept: w17 task 2 -- change ONE CONTROL key and nothing else.

The person asks only for the basis change (def2-SVP -> def2-TZVP).
Accepted when, in the fixture CONTROL:
1. basis is def2-TZVP (the requested change happened);
2. functional PBE0 and jobtype opt are still there UNCHANGED (no
   drive-by edits);
3. the file is otherwise the same text -- every other line byte-equal.
An agent that "improves" the CONTROL while it is there (adds a
solvent, changes the functional, reformats) fails here even though
the requested change happened. The prompt says "change one thing".
"""
from __future__ import annotations

import sys
from pathlib import Path

_FIXTURE = "w17_dash_fixture"
_LINES_BEFORE = [
    "# CONTROL (fixture)",
    "functional PBE0",
    "basis def2-SVP",
    "jobtype opt",
]


def check(workspace: Path) -> tuple[bool, str]:
    ws = Path(workspace) / _FIXTURE
    control = ws / "CONTROL"
    if not control.is_file():
        return False, "CONTROL is missing (deleted?)"
    after = [ln.rstrip("\n") for ln in
             control.read_text(encoding="utf-8").splitlines()
             if ln.strip()]
    if not after:
        return False, "CONTROL is empty"
    # The requested change, exactly once.
    basis_lines = [ln for ln in after
                   if ln.split()[:1] == ["basis"]]
    if basis_lines != ["basis def2-TZVP"]:
        return False, (f"basis not changed exactly to def2-TZVP: "
                       f"{basis_lines}")
    # Everything else unchanged.
    expected = ["# CONTROL (fixture)", "functional PBE0",
                "basis def2-TZVP", "jobtype opt"]
    if [ln.split("#")[0].strip() or ln for ln in after] != expected:
        return False, (f"CONTROL differs beyond the one basis line: "
                       f"{after!r}")
    assert _LINES_BEFORE  # shape guard: the fixture lines are pinned
    return True, "accepted: only the basis line changed"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
