#!/usr/bin/env python3
"""Shared builder for the two w17_method web tasks.

Not a setup on its own: the w17_method_web_*_setup.py entry scripts
call ``build`` through generalist_repo.main, which supplies the setup
contract the suite expects (``usage:`` on a missing workspace argument,
``refusing to overwrite`` on a second run, ``prepared <path>`` on
success) and the reviewed fixture-builder git() the safe_git ALLOWED
list covers.

Both variants seed an archived page (plain .html/.md file, no network
— the runner is offline and the task says so) plus a decoy summary
file that contradicts the page. The method under test: read the
SOURCE, not the summary.

variant="manual": docs/manual.html states the command flag is
``--baseline``. notes/summary.txt claims ``--bessel``. Task: answer in
docs/answer.txt with the exact flag name taken from the manual.

variant="changelog": docs/CHANGELOG.md states the option was REMOVED
in v2.4 and renamed to ``strict``. notes/summary.txt claims it still
exists as ``lenient``. Task: docs/answer.txt states the option's
status and its new name.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generalist_repo import (  # noqa: E402
    Precondition, commit, new_repo, write,
)

# A fixture repo must not depend on the user's home git configuration:
# sandboxed test environments deny ~/.gitconfig (fatal "Permission
# denied" on init). The shared generalist_repo core runs plain git;
# these redirects make it home-independent. Set before any git call
# runs (generalist_repo reads os.environ per call).
import os  # noqa: E402

os.environ.setdefault("GIT_CONFIG_GLOBAL", os.devnull)
os.environ.setdefault("GIT_CONFIG_SYSTEM", os.devnull)
# The XDG excludes file (~/.config/git/ignore) is a third home path
# git reads even with the config redirects above; cage denial is
# fatal ("cannot use ... as an exclude file").
os.environ.setdefault("XDG_CONFIG_HOME", os.devnull)

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
- renamed the `lenient` option to `strict`; `lenient` no longer exists
  (it was removed with v2.4).

## v2.3
- added the `lenient` option.
"""

_SUMMARY_MANUAL = "The flag for a baseline comparison is --bessel.\n"
_SUMMARY_CHANGELOG = ("The lenient option still exists in v2.4; it is "
                      "unchanged since v2.3.\n")

_BRIEF_MANUAL = ("Goal: answer from the archived manual (docs/manual.html), "
                 "NOT from notes/summary.txt. Success criterion: "
                 "docs/answer.txt holds the exact command-line flag the "
                 "manual gives for a baseline comparison.\n")
_BRIEF_CHANGELOG = ("Goal: answer from the archived changelog "
                    "(docs/CHANGELOG.md), NOT from notes/summary.txt. "
                    "Success criterion: docs/answer.txt states whether the "
                    "lenient option still exists in v2.4 and, if not, the "
                    "name of its replacement.\n")


def build(workspace: str, *, variant: str) -> Path:
    if variant not in ("manual", "changelog"):
        raise Precondition(f"unknown variant {variant!r}")
    repo = new_repo(workspace, "web_notes")
    (repo / "docs").mkdir()
    (repo / "notes").mkdir()
    if variant == "manual":
        write(repo, {"docs/manual.html": _MANUAL_HTML,
                     "notes/summary.txt": _SUMMARY_MANUAL,
                     "m1_brief.txt": _BRIEF_MANUAL})
        precondition = "--baseline" in _MANUAL_HTML
    else:
        write(repo, {"docs/CHANGELOG.md": _CHANGELOG_MD,
                     "notes/summary.txt": _SUMMARY_CHANGELOG,
                     "m1_brief.txt": _BRIEF_CHANGELOG})
        precondition = ("lenient" in _CHANGELOG_MD
                        and "strict" in _CHANGELOG_MD)
    commit(repo, "seed")
    if not precondition:
        raise Precondition(f"{variant} precondition broken")
    print(f"seeded web_notes/{variant}: preconditions hold")
    return repo


if __name__ == "__main__":
    # Entry scripts import this module; direct execution is a mistake.
    print("not a setup on its own; use w17_method_web_*_setup.py",
          file=sys.stderr)
    sys.exit(2)
