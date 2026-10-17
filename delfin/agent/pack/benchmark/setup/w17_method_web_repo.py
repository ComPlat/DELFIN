#!/usr/bin/env python3
"""Build the local web fixtures for the two w17_method web tasks.

Not a setup on its own: w17_method_web_*_setup.py call
``build(workspace, variant)``.

Both variants seed an archived page (plain .html file, no network — the
runner is offline and the task says so) plus a decoy summary file that
contradicts the page. The method under test: read the SOURCE, not the
summary.

variant="manual": docs/manual.html states the command flag is
``--baseline``. notes/summary.txt claims ``--bessel``. Task: answer in
docs/answer.txt with the exact flag name taken from the manual.

variant="changelog": docs/CHANGELOG.md states the option was REMOVED in
v2.4 and renamed to ``strict``. notes/summary.txt claims it still exists
as ``lenient``. Task: docs/answer.txt states the option's status and its
new name.

Each variant asserts its preconditions after building; a failed one exits
non-zero so the attempt is reported as unmeasured, never as a model
failure.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

_MANUAL_HTML = """<!DOCTYPE html>
<html><head><title>toolkit manual</title></head>
<body>
<h1>toolkit 1.0 manual</h1>
<p>To compare a measurement against a stored reference run, pass the
<code>--baseline</code> flag together with the path of the reference
file.</p>
</body></html>
"""

_CHANGELOG_MD = """# toolkit changelog

## v2.4 (latest)
- The `lenient` option was REMOVED. Its replacement is the `strict`
  option, which behaves the same way under a new name.

## v2.3
- Introduced `lenient`.
"""

_SUMMARY_MANUAL = ("summary: use --bessel for baseline comparison\n")
_SUMMARY_CHANGELOG = ("summary: the lenient option still exists in v2.4; "
                      "use --lenient on the command line\n")

_BRIEF_MANUAL = ("Goal: answer from the archived manual (docs/manual.html), "
                 "NOT from notes/summary.txt. Success criterion: "
                 "docs/answer.txt holds the exact command-line flag the "
                 "manual gives for a baseline comparison.\n")
_BRIEF_CHANGELOG = ("Goal: answer from the archived changelog "
                    "(docs/CHANGELOG.md), NOT from notes/summary.txt. "
                    "Success criterion: docs/answer.txt states whether the "
                    "lenient option still exists in v2.4 and, if not, the "
                    "name of its replacement.\n")


def _run(cmd: list[str], cwd: Path) -> None:
    # Sandboxed test environments deny access to ~/.gitconfig; a fixture
    # repo must not depend on a home config. GIT_CONFIG_GLOBAL redirects.
    env = dict(os.environ)
    env["GIT_CONFIG_GLOBAL"] = os.devnull
    env["GIT_CONFIG_SYSTEM"] = os.devnull
    # git also picks up ~/.config/git/ignore as the default excludes file
    # (XDG); redirect that too, into the fixture directory itself.
    env["XDG_CONFIG_HOME"] = str(cwd)
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True,
                          env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stdout[-300:]}"
                           f"{proc.stderr[-300:]}")


def build(workspace: Path, *, variant: str) -> None:
    if variant not in ("manual", "changelog"):
        raise ValueError(f"unknown variant {variant!r}")
    root = workspace / "web_notes"
    if root.exists():
        raise RuntimeError("web_notes/ already exists; refusing to seed "
                           "over an earlier attempt")
    root.mkdir(parents=True)
    _run(["git", "init", "-q"], root)
    (root / "docs").mkdir()
    (root / "notes").mkdir()
    if variant == "manual":
        (root / "docs" / "manual.html").write_text(_MANUAL_HTML,
                                                   encoding="utf-8")
        (root / "notes" / "summary.txt").write_text(_SUMMARY_MANUAL,
                                                    encoding="utf-8")
        brief = _BRIEF_MANUAL
        precondition = "--baseline" in _MANUAL_HTML
    else:
        (root / "docs" / "CHANGELOG.md").write_text(_CHANGELOG_MD,
                                                    encoding="utf-8")
        (root / "notes" / "summary.txt").write_text(_SUMMARY_CHANGELOG,
                                                    encoding="utf-8")
        brief = _BRIEF_CHANGELOG
        precondition = ("lenient" in _CHANGELOG_MD
                        and "strict" in _CHANGELOG_MD)
    (root / "m1_brief.txt").write_text(brief, encoding="utf-8")
    _run(["git", "add", "-A"], root)
    _run(["git", "-c", "user.email=fixture@example.com",
          "-c", "user.name=fixture", "commit", "-q", "-m", "seed"], root)
    if not precondition:
        raise RuntimeError(f"{variant} precondition broken")
    print(f"seeded web_notes/{variant}: preconditions hold")


if __name__ == "__main__":
    build(Path(sys.argv[1]),
          variant=sys.argv[2] if len(sys.argv) > 2 else "manual")
