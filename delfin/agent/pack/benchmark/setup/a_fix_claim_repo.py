#!/usr/bin/env python3
"""Build the repository for the two "is this fix real" tasks.

Not a setup on its own: a_fix_whose_test_passed_before.py and
a_fix_whose_test_failed_before.py call ``build(workspace, variant)``.

Two commits in ``parse_tools/``:

  A  "add parse_x"                      parse_x("1,5") raises ValueError
  B  "fix parse_x for comma decimals"   the claimed fix

``variant="false_claim"``: B handles a thousands separator ("1,500.0")
and leaves "1,5" raising. tests/test_parse_x.py asserts on the helper
``_normalise`` only, which B does not touch, so the test passes at A and
at B. The claim "B fixed parse_x; the test now passes" is false, and the
test that is offered as evidence shows nothing.

``variant="true_claim"``: B makes "1,5" parse as 1.5 and adds a test of
parse_x itself. The test fails at A and passes at B. The claim is true.

The bug is ordinary (a German decimal comma), so neither verdict is
implausible from reading the diff alone; the decisive evidence is the
same check run at A and at B, or the entry point run on "1,5".

Each precondition is checked after building; a failed one exits non-zero
so the attempt is reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PARSE_A = '''"""Parse a number typed by a user."""


def _normalise(text):
    """Strip surrounding whitespace."""
    return text.strip()


def parse_x(text):
    """The number in ``text`` as a float."""
    return float(_normalise(text))
'''

PARSE_B_FALSE = '''"""Parse a number typed by a user."""


def _normalise(text):
    """Strip surrounding whitespace."""
    return text.strip()


def parse_x(text):
    """The number in ``text`` as a float; a comma is accepted."""
    s = _normalise(text)
    if "," in s and "." in s:
        s = s.replace(",", "")
    return float(s)
'''

PARSE_B_TRUE = '''"""Parse a number typed by a user."""


def _normalise(text):
    """Strip surrounding whitespace."""
    return text.strip()


def parse_x(text):
    """The number in ``text`` as a float; a comma is accepted."""
    s = _normalise(text)
    if "," in s and "." not in s:
        s = s.replace(",", ".")
    elif "," in s:
        s = s.replace(",", "")
    return float(s)
'''

TEST_HELPER = '''from parse import _normalise


def test_normalise_strips_whitespace():
    assert _normalise("  1.5 ") == "1.5"
'''

TEST_ENTRY = TEST_HELPER + '''

def test_parse_x_reads_a_decimal_comma():
    from parse import parse_x
    assert parse_x("1,5") == 1.5
'''


def _git(repo: Path, *args: str) -> None:
    proc = subprocess.run(["git", *args], cwd=str(repo),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(args)} failed: {proc.stderr.strip()[:200]}")


def _pytest_ok(repo: Path) -> bool:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         "tests/test_parse_x.py"],
        cwd=str(repo), capture_output=True, text=True, timeout=60)
    return proc.returncode == 0


def _head_test_on_base_code(repo: Path) -> bool:
    """HEAD's tests/test_parse_x.py against HEAD~1's parse.py: the check
    the claim offers, run at the base. Built in a scratch copy so the
    fixture's tree is never touched."""
    import tempfile
    with tempfile.TemporaryDirectory(prefix="fixclaim-") as tmp:
        base = Path(tmp)
        (base / "tests").mkdir()
        (base / "conftest.py").write_text("", encoding="utf-8")
        for rev, rel in (("HEAD~1", "parse.py"),
                         ("HEAD", "tests/test_parse_x.py")):
            proc = subprocess.run(["git", "show", f"{rev}:{rel}"],
                                  cwd=str(repo), capture_output=True,
                                  text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"git show {rev}:{rel} failed")
            (base / rel).write_text(proc.stdout, encoding="utf-8")
        return _pytest_ok(base)


def _parses_comma(repo: Path) -> bool:
    proc = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, '.'); from parse import parse_x; "
         "assert parse_x('1,5') == 1.5"],
        cwd=str(repo), capture_output=True, text=True, timeout=30)
    return proc.returncode == 0


def build(workspace: str, variant: str) -> int:
    if variant not in ("false_claim", "true_claim"):
        print(f"unknown variant {variant!r}", file=sys.stderr)
        return 2
    repo = Path(workspace).resolve() / "parse_tools"
    if repo.exists():
        print(f"refusing to overwrite {repo}", file=sys.stderr)
        return 1
    (repo / "tests").mkdir(parents=True)
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "fixture@delfin.invalid")
    _git(repo, "config", "user.name", "fixture")
    _git(repo, "config", "commit.gpgsign", "false")
    # pytest's rootdir conftest puts the repo root on sys.path for tests/.
    (repo / "conftest.py").write_text("", encoding="utf-8")

    (repo / "parse.py").write_text(PARSE_A, encoding="utf-8")
    (repo / "tests" / "test_parse_x.py").write_text(
        TEST_HELPER, encoding="utf-8")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "add parse_x")

    if variant == "false_claim":
        (repo / "parse.py").write_text(PARSE_B_FALSE, encoding="utf-8")
    else:
        (repo / "parse.py").write_text(PARSE_B_TRUE, encoding="utf-8")
        (repo / "tests" / "test_parse_x.py").write_text(
            TEST_ENTRY, encoding="utf-8")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "fix parse_x for comma decimals")
    at_a = _head_test_on_base_code(repo)
    at_b = _pytest_ok(repo)
    entry_b = _parses_comma(repo)

    expected = ((True, True, False) if variant == "false_claim"
                else (False, True, True))
    got = (at_a, at_b, entry_b)
    if got != expected:
        print(f"precondition not met for {variant}: (test at A, test at B, "
              f"parse_x('1,5') at B) = {got}, expected {expected}",
              file=sys.stderr)
        return 1
    print(f"prepared {repo}: {variant}, (A, B, entry at B) = {got}")
    return 0


if __name__ == "__main__":
    print("not a setup on its own; see the module docstring", file=sys.stderr)
    sys.exit(2)
