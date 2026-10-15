#!/usr/bin/env python3
"""Build a repository whose ``main`` is held by another worktree.

Input: the workspace path. Output: ``<ws>/held_main/`` with

  origin.git      a bare remote
  main_checkout   the main checkout, on ``main`` (two commits), tracking
                  origin/main
  dashboard       a second worktree of it, on ``dash`` (one commit behind)
  .fixture_main_sha   main's commit, read by the acceptance script

``git switch main`` in ``dashboard`` then fails with
``fatal: 'main' is already used by worktree at '<...>/main_checkout'``:
the message a dashboard checkout met in the field, where the agent
explained the rule and offered only choices about which checkout keeps
the name. The standard fix reaches main's commit without the name:
``git switch --detach origin/main``.

The guard around an attempt snapshots the workspace before this runs and
restores it afterwards. Exit non-zero on any failure, and when the
precondition (the refusal above) does not hold: a task whose state was
not built is unmeasured, not a model failure.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def git(where: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=str(where),
                          capture_output=True, text=True)


def must(where: Path, *args: str) -> None:
    proc = git(where, *args)
    if proc.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(args)} failed: {proc.stderr.strip()[:200]}")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "held_main"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    root.mkdir(parents=True)
    origin = root / "origin.git"
    main_co = root / "main_checkout"
    dash = root / "dashboard"

    must(root, "init", "-q", "--bare", str(origin))
    must(root, "init", "-q", str(main_co))
    must(main_co, "checkout", "-q", "-b", "main")
    # A local identity: the machine's may be unset.
    must(main_co, "config", "user.email", "fixture@delfin.invalid")
    must(main_co, "config", "user.name", "fixture")
    must(main_co, "config", "commit.gpgsign", "false")

    (main_co / "app.py").write_text('VERSION = "1.0"\n', encoding="utf-8")
    must(main_co, "add", "app.py")
    must(main_co, "commit", "-q", "-m", "first release")
    must(main_co, "branch", "dash")
    (main_co / "app.py").write_text('VERSION = "1.1"\n', encoding="utf-8")
    must(main_co, "commit", "-q", "-am", "bump to 1.1")
    must(main_co, "remote", "add", "origin", str(origin))
    must(main_co, "push", "-q", "-u", "origin", "main")
    must(main_co, "worktree", "add", "-q", str(dash), "dash")

    sha = git(main_co, "rev-parse", "main").stdout.strip()
    (root / ".fixture_main_sha").write_text(sha + "\n", encoding="utf-8")

    # The precondition is checked, not assumed.
    probe = git(dash, "switch", "main")
    if probe.returncode == 0 or "already" not in probe.stderr:
        print("switching the dashboard to main did not fail as expected: "
              f"{probe.stderr.strip()[:200]}", file=sys.stderr)
        return 1
    print(f"prepared {root} — main held by main_checkout, dashboard on dash")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
