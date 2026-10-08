"""DELFIN required git, gh and pytest, and checked for none of them.

Three prerequisites the product depends on, none of which it verified or
named a remedy for.

**pytest.** Declared only inside the ``dev`` extra, which the default
install does not select, so every site install produced an environment
whose test tool could not start. `python -m pytest` with no pytest exits 1
and writes no report, which reached the agent as "no report file produced"
-- a parse complaint about a run that never began. A field report describes
what follows: told that, the agent builds a runner of its own, in a venv
or as a wrapper script under the home directory.

**git and gh.** No check anywhere for either binary, for a GitHub login,
or for a commit identity. git's absence was noticed only as a side effect
of the remote check -- reported as "git remote", with no fix string, and
only where an origin remote exists. gh was never probed at all, while the
write gate REQUIRES the pull-request route, so the one tool the sanctioned
path needs was the one nothing checked.

Universal: every binary is shadowed through ``shutil.which`` and
``subprocess.run``, so these assertions hold on a host that has all three
and on one that has none. Nothing is skipped for a missing tool -- a check
that disappears where the tool is absent is no check.
"""

from __future__ import annotations

import subprocess
import tomllib
from pathlib import Path

import pytest

from delfin.agent import doctor as D
from delfin.agent import test_runner as TR
from delfin.agent.push_diagnosis import diagnose

_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# The default install ships a test runner
# ---------------------------------------------------------------------------

def _extras() -> dict:
    data = tomllib.loads((_ROOT / "pyproject.toml").read_text())
    return data["project"]["optional-dependencies"]


def _names(reqs) -> set[str]:
    out = set()
    for r in reqs:
        out.add(str(r).split(";")[0].strip().lower()
                .split("[")[0].split(">")[0].split("=")[0].split("<")[0]
                .strip())
    return out


def test_pytest_is_declared_in_an_extra_of_its_own():
    assert "pytest" in _names(_extras()["test"])


def test_the_dev_extra_still_brings_it():
    """Splitting it out must not take it away from a developer install.

    Written out rather than pulled in as a self-referential extra, which
    resolves against PyPI -- so the two lists have to be kept in step, and
    that is what this asserts.
    """
    assert _names(_extras()["test"]) <= _names(_extras()["dev"])


def test_the_default_install_selects_it():
    """Read out of the installer, not assumed: this default is what every
    site install and the dashboard's install button both use."""
    script = (_ROOT / "delfin" / "installers" / "install_delfin.sh").read_text()
    import re
    m = re.search(r'DELFIN_EXTRAS="\$\{DELFIN_EXTRAS-([^}"]*)\}"', script)
    assert m, "the installer's default extras could not be read"
    selected = {e.strip() for e in m.group(1).split(",") if e.strip()}
    assert "test" in selected, (
        f"default extras {sorted(selected)} bring no pytest, so the agent's "
        "test tool cannot start after a default install")
    assert selected <= set(_extras()), "an extra that pyproject does not define"


# ---------------------------------------------------------------------------
# A missing pytest is reported as the environment, not as a parse failure
# ---------------------------------------------------------------------------

def test_the_runner_names_the_interpreter_that_has_no_pytest(tmp_path,
                                                             monkeypatch):
    spawned = []

    def fake_run(argv, **kw):
        spawned.append(list(argv))
        return subprocess.CompletedProcess(argv, 1, "", "No module named pytest")

    monkeypatch.setattr(subprocess, "run", fake_run)
    out = TR.run_tests(tmp_path, python="/usr/bin/python3")
    assert out["status"] == "error"
    assert "pytest is not installed" in out["error"]
    assert "/usr/bin/python3" in out["error"]
    assert "pip install" in out["fix"]


def test_it_does_not_launch_a_suite_it_knows_cannot_run(tmp_path, monkeypatch):
    monkeypatch.setattr(
        subprocess, "run",
        lambda argv, **kw: subprocess.CompletedProcess(argv, 1, "", ""))
    called = []
    from delfin.agent import contained_run as CR
    monkeypatch.setattr(CR, "run", lambda *a, **k: called.append(a))
    TR.run_tests(tmp_path, python="/usr/bin/python3")
    assert called == [], "pytest was launched although it is absent"


def test_it_tells_the_agent_not_to_build_its_own_runner(tmp_path, monkeypatch):
    """The behaviour the field report describes. The remedy has to be in
    the answer, or the model supplies one."""
    monkeypatch.setattr(
        subprocess, "run",
        lambda argv, **kw: subprocess.CompletedProcess(argv, 1, "", ""))
    out = TR.run_tests(tmp_path, python="/usr/bin/python3")
    note = (out.get("note") or "").lower()
    assert "do not" in note
    assert "wrapper" in note or "interpreter" in note


def test_a_present_pytest_is_not_in_the_way(tmp_path, monkeypatch):
    """The preflight must not become the thing that stops every run."""
    monkeypatch.setattr(
        subprocess, "run",
        lambda argv, **kw: subprocess.CompletedProcess(argv, 0, "", ""))
    assert TR._pytest_missing_from("/usr/bin/python3") == ""


def test_an_interpreter_that_cannot_be_asked_counts_as_present(monkeypatch):
    """A preflight that cannot reach the interpreter must not stand in for
    a verdict about the suite."""
    def boom(*a, **k):
        raise OSError("no exec")

    monkeypatch.setattr(subprocess, "run", boom)
    assert TR._pytest_missing_from("/usr/bin/python3") == ""


# ---------------------------------------------------------------------------
# git and gh: four questions, four rows, every one with a remedy
# ---------------------------------------------------------------------------

def _rows(monkeypatch, *, git=None, gh=None, identity=True, logged_in=True):
    """Drive the check with the host's tools shadowed."""
    import shutil

    monkeypatch.setattr(
        shutil, "which",
        lambda name: {"git": git, "gh": gh}.get(name))

    def fake_run(argv, **kw):
        argv = list(argv)
        if argv[:2] == ["git", "config"]:
            return subprocess.CompletedProcess(
                argv, 0 if identity else 1, "me\n" if identity else "", "")
        if argv[:3] == ["gh", "auth", "status"]:
            return subprocess.CompletedProcess(
                argv, 0 if logged_in else 1, "", "")
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(D.subprocess, "run", fake_run)
    return {r["check"]: r for r in D._check_git_tooling({"workspace": "."})}


def test_the_check_is_registered():
    assert any(attr == "_check_git_tooling" for _, attr in D._CHECK_ATTRS)


def test_all_four_questions_are_asked(monkeypatch):
    rows = _rows(monkeypatch, git="/usr/bin/git", gh="/usr/bin/gh")
    assert set(rows) == {"git installed", "git identity",
                         "gh installed", "gh authenticated"}


@pytest.mark.parametrize("name", ["git installed", "git identity",
                                  "gh installed", "gh authenticated"])
def test_a_failure_always_carries_a_remedy(monkeypatch, name):
    rows = _rows(monkeypatch, git=None, gh=None, identity=False,
                 logged_in=False)
    if name not in rows:
        # git absent hides the identity question; ask it the other way.
        rows = _rows(monkeypatch, git="/usr/bin/git", gh="/usr/bin/gh",
                     identity=False, logged_in=False)
    row = rows[name]
    assert row["status"] != D.PASS
    assert row["fix"], f"{name} reports a problem with no remedy"


def test_presence_and_authentication_are_different_questions(monkeypatch):
    """The assertion the old code could not make at all: gh is there and
    still cannot open a pull request."""
    rows = _rows(monkeypatch, git="/usr/bin/git", gh="/usr/bin/gh",
                 logged_in=False)
    assert rows["gh installed"]["status"] == D.PASS
    assert rows["gh authenticated"]["status"] != D.PASS
    assert "gh auth login" in rows["gh authenticated"]["fix"]


def test_the_identity_remedy_names_the_command(monkeypatch):
    rows = _rows(monkeypatch, git="/usr/bin/git", gh="/usr/bin/gh",
                 identity=False)
    assert "git config" in rows["git identity"]["fix"]


def test_a_missing_gh_says_why_it_is_needed(monkeypatch):
    rows = _rows(monkeypatch, git="/usr/bin/git", gh=None)
    assert "pull request" in rows["gh installed"]["fix"]


def test_everything_present_passes(monkeypatch):
    rows = _rows(monkeypatch, git="/usr/bin/git", gh="/usr/bin/gh")
    assert all(r["status"] == D.PASS for r in rows.values()), rows


def test_the_check_never_raises_with_nothing_installed(monkeypatch):
    """The module's stated contract."""
    rows = _rows(monkeypatch, git=None, gh=None)
    assert rows


# ---------------------------------------------------------------------------
# And the failure text maps to an action
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text,expect", [
    ("gh: command not found", "GitHub CLI is not installed"),
    ("bash: line 1: git: command not found", "git is not installed"),
    ("To get started with GitHub CLI, please run: gh auth login",
     "no active login"),
    ("*** Please tell me who you are.", "no commit identity"),
])
def test_a_tool_failure_is_named_not_passed_through_raw(text, expect):
    d = diagnose(text)
    assert d is not None, text
    assert expect in d.cause, (text, d.cause)


@pytest.mark.parametrize("text,expect", [
    ("Could not resolve hostname github.com", "no DNS"),
    ("fatal: Authentication failed for https://github.com/x/y",
     "no usable credentials"),
    ("! [rejected] main -> main (non-fast-forward)",
     "commits this branch does not"),
])
def test_the_existing_diagnoses_still_win_where_they_should(text, expect):
    d = diagnose(text)
    assert d is not None, text
    assert expect in d.cause, (text, d.cause)
