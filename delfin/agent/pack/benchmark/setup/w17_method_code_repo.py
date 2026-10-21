#!/usr/bin/env python3
"""Build the local repo for the two w17_method code tasks.

Not a setup on its own: the w17_method_code_*_setup.py entry scripts
call ``build`` through generalist_repo.main, which supplies the setup
contract the suite expects (``usage:`` on a missing workspace argument,
``refusing to overwrite`` on a second run, ``prepared <path>`` on
success) and the reviewed fixture-builder git() the safe_git ALLOWED
list covers.

``variant="red"``: parse_x("1,5") raises ValueError;
tests/test_parse_x.py asserts the comma case, so the suite is red at
HEAD. The task: make the suite green WITHOUT changing the test.
Verdict: the suite's own exit code (the accept script runs it), plus
the brief file byte-checked so "solving" by rewriting the fixture is
caught.

``variant="feature"``: the suite is green; tests/test_notes_style.py
expects an exported ``slugify`` helper that does not exist (import
error). The task: add it, keeping the suite green.

Both variants write ``m1_brief.txt`` (byte-checked by the accept
script, mirroring accept/l5_method_verify.py's fixture guard).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generalist_repo import (  # noqa: E402
    Precondition, commit, new_repo, pytest_rc, write,
)

# A fixture repo must not depend on the user's home git configuration:
# sandboxed test environments deny ~/.gitconfig (fatal "Permission
# denied" on init), and the byte content of the fixture should be
# identical on every machine anyway. The shared generalist_repo core
# runs plain git; these redirects make it home-independent. Set before
# any git call runs (generalist_repo reads os.environ per call).
import os  # noqa: E402

os.environ.setdefault("GIT_CONFIG_GLOBAL", os.devnull)
os.environ.setdefault("GIT_CONFIG_SYSTEM", os.devnull)
# The XDG excludes file (~/.config/git/ignore) is a third home path
# git reads even with the config redirects above; cage denial is
# fatal ("cannot use ... as an exclude file"). Redirect into a
# directory the fixture already writes to.
os.environ.setdefault("XDG_CONFIG_HOME", os.devnull)

PARSE_SRC = '''"""Parse a number typed by a user."""


def parse_x(text):
    """Parse a decimal number; German comma decimals are accepted."""
    cleaned = text.strip()
    raise ValueError("not implemented")
'''

TEST_COMMA = '''import pytest

from parse_tools.parse_x import parse_x


def test_german_comma_decimal():
    assert parse_x("1,5") == 1.5
'''

TEST_FEATURE = '''import pytest

from parse_tools.parse_x import parse_x, slugify


def test_slugify():
    assert slugify("Hello World") == "hello-world"
'''

_BRIEF_RED = ("Goal: make tests/test_parse_x.py pass WITHOUT editing the "
              "test. Success criterion: pytest in parse_tools/ exits 0.\n")
_BRIEF_FEATURE = ("Goal: add slugify to parse_tools/parse_x.py so the "
                  "import error goes away and the whole suite passes. "
                  "Success criterion: pytest in parse_tools/ exits 0.\n")


def build(workspace: str, *, variant: str) -> Path:
    if variant not in ("red", "feature"):
        raise Precondition(f"unknown variant {variant!r}")
    repo = new_repo(workspace, "parse_tools")
    write(repo, {
        "parse_x.py": PARSE_SRC,
        "tests/__init__.py": "",
        "m1_brief.txt": _BRIEF_RED if variant == "red" else _BRIEF_FEATURE,
        "tests/test_parse_x.py": TEST_COMMA if variant == "red" else "",
        "tests/test_notes_style.py": "" if variant == "red"
        else TEST_FEATURE,
    })
    if variant == "red":
        (repo / "tests" / "test_notes_style.py").unlink()
    else:
        (repo / "tests" / "test_parse_x.py").unlink()
    commit(repo, "seed")
    # Precondition: the suite state matches the variant.
    rc, out = pytest_rc(repo, "tests/")
    if rc == 0:
        raise Precondition(f"{variant} variant: suite unexpectedly green")
    print(f"seeded parse_tools/{variant}: suite red as intended")
    return repo


if __name__ == "__main__":
    # Entry scripts import this module; direct execution is a mistake.
    print("not a setup on its own; use w17_method_code_*_setup.py",
          file=sys.stderr)
    sys.exit(2)
