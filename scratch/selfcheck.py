"""Does an arm's own branch pass its own gate?

    PYTHONPATH=. python selfcheck.py

Run from the arm's worktree. Calls ``land.check`` on THIS worktree against
origin/main, with the real configured identity and the real agent address --
no stub. A gate that refuses the branch it was written on, or that lands it
while its own PR.md proves nothing, is measured here and nowhere else.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, os.getcwd())
import delfin  # noqa: E402
from delfin.agent.swarm import land as L  # noqa: E402

print(f"# module: {os.path.realpath(delfin.__file__)}")
print(f"# identity: {L.identity_email()!r}")
co = getattr(L, "coauthor_email", lambda: "(none: no such function)")()
print(f"# agent:    {co!r}")

wt = Path.cwd()
for call in (lambda: L.check(wt, "origin/main", run="w54", tid="g1-int"),
             lambda: L.check(wt, "origin/main")):
    try:
        ok, why = call()
    except TypeError:
        continue
    print(f"# own branch: {'PASSES' if ok else 'REFUSED'} -- {why[:400]}")
    break
