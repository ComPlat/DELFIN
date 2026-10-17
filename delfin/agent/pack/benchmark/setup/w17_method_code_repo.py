#!/usr/bin/env python3
"""Build the local repo for the two w17_method code tasks.

Not a setup on its own: w17_method_code_setup.py calls ``build``.

``variant="red"``: parse_x("1,5") raises ValueError; tests/test_parse_x.py
asserts the comma case, so the suite is red at HEAD. The task: make the
suite green WITHOUT changing the test. Verdict: the suite's own exit code
(the accept script runs it), plus the brief file byte-checked so "solving"
by rewriting the fixture is caught.

``variant="feature"``: the suite is green; tests/test_notes_style.py
expects an exported ``slugify`` helper that does not exist (import error).
The task: add it, keeping the suite green.

Both variants write ``m1_brief.txt`` (byte-checked by the accept script,
mirroring accept/l5_method_verify.py's fixture guard).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

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


def _run(cmd: list[str], cwd: Path) -> None:
    # Sandboxed test environments deny access to the user's ~/.gitconfig;
    # a fixture repo must not depend on a home config (also keeps it
    # byte-identical on every machine). GIT_CONFIG_GLOBAL redirects it.
    env = dict(os.environ)
    env["GIT_CONFIG_GLOBAL"] = os.devnull
    env["GIT_CONFIG_SYSTEM"] = os.devnull
    # git also picks up ~/.config/git/ignore as the default excludes file
    # (XDG); redirect that too, into the fixture directory itself.
    env["XDG_CONFIG_HOME"] = str(cwd)
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True,
                          env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stdout[-400:]}"
                           f"{proc.stderr[-400:]}")


def build(workspace: Path, *, variant: str) -> None:
    if variant not in ("red", "feature"):
        raise ValueError(f"unknown variant {variant!r}")
    root = workspace / "parse_tools"
    root.mkdir(parents=True, exist_ok=True)
    src = root / "parse_x.py"
    if src.exists():
        raise RuntimeError("parse_tools/ already exists; refusing to seed "
                           "over an earlier attempt")
    _run(["git", "init", "-q"], root)
    src.parent.mkdir(parents=True, exist_ok=True)
    src.write_text(PARSE_SRC, encoding="utf-8")
    tests = root / "tests"
    tests.mkdir(exist_ok=True)
    if variant == "red":
        (tests / "test_parse_x.py").write_text(TEST_COMMA, encoding="utf-8")
        (root / "m1_brief.txt").write_text(_BRIEF_RED, encoding="utf-8")
    else:
        (tests / "test_notes_style.py").write_text(TEST_FEATURE,
                                                   encoding="utf-8")
        (root / "m1_brief.txt").write_text(_BRIEF_FEATURE, encoding="utf-8")
    (tests / "__init__.py").write_text("", encoding="utf-8")
    _run(["git", "add", "-A"], root)
    _run(["git", "-c", "user.email=fixture@example.com",
          "-c", "user.name=fixture", "commit", "-q", "-m", "seed"], root)
    # Precondition: the suite state matches the variant.
    suite = subprocess.run(
        [sys.executable, "-m", "pytest", "-q"],
        cwd=str(root), capture_output=True, text=True)
    if variant == "red" and suite.returncode == 0:
        raise RuntimeError("red variant: suite unexpectedly green")
    if variant == "feature" and suite.returncode == 0:
        raise RuntimeError("feature variant: suite unexpectedly green")
    print(f"seeded parse_tools/{variant}: suite red as intended")


if __name__ == "__main__":
    build(Path(sys.argv[1]), variant=sys.argv[2] if len(sys.argv) > 2 else "red")
