"""Phase 2 of U2: a push/pr readiness function the gate can consult.

The doctor already *reports* git/gh state (``_check_git_tooling``,
``_check_push``) but nothing is a boolean gate another part can ask
*before* ``git push`` / ``gh pr create``. The gate that guards those two
commands must, before letting one through, know whether this host can
actually push: is git there, is an identity set, does the remote answer,
is gh there and logged in.

``ready_for_push`` answers that. It is a thin composition of the two
probes the doctor report already runs -- ``_check_git_tooling`` and
``_check_push`` -- so there is one set of facts about this checkout and a
second, divergent copy cannot creep in. A test needs no live tool and no
network: it shadows ``shutil.which`` (git/gh on PATH) and
``subprocess.run`` (config / ls-remote / gh auth status), the doctor's own
injection point (see tests/test_delfin_checks_the_tools_it_requires.py).

The gate reads the return: every row PASS means ready; any non-PASS row
names the first thing in the way and its fix. Every row carries prose
``fix`` only: installing git or gh is a system-package change and the
login is the user's (``! gh auth login``), so per the ``_row`` contract
(doctor.py:54-74) no readiness row declares a ``command`` the agent would
execute.
"""

from __future__ import annotations

import shutil
import subprocess

from delfin.agent import doctor as D


def _fake_run(*, git=True, identity=True, remote=True, has_origin=True,
              gh=True, logged_in=True):
    """A fake ``subprocess.run`` shaped like the real one, keyed on argv.

    ``_check_git_tooling`` / ``_check_push`` call it through the module
    attribute ``D.subprocess.run`` that a test monkeypatches, exactly as
    the doctor's own tests do. When git is said to be absent the git
    probes raise ``FileNotFoundError``, which is what a real shell does
    with no git on PATH and what ``_check_push`` is built to catch.
    """
    def run(argv, **kw):
        argv = list(argv)
        if not git and argv and argv[0] == "git":
            raise FileNotFoundError("no git on PATH")
        if argv[:2] == ["git", "config"] and argv[2] == "--get":
            field = argv[3]
            if field == "credential.helper":
                return subprocess.CompletedProcess(argv, 0, "", "")
            if not identity:
                return subprocess.CompletedProcess(argv, 1, "", "")
            return subprocess.CompletedProcess(argv, 0, f"{field} set\n", "")
        if argv[:3] == ["git", "remote", "get-url"]:
            if not has_origin:
                err = "error: No such remote 'origin'"
                return subprocess.CompletedProcess(argv, 2, "", err)
            return subprocess.CompletedProcess(
                argv, 0, "git@github.com:o/r.git\n", "")
        if "ls-remote" in argv:
            if not remote:
                err = "fatal: Could not read from remote repository."
                return subprocess.CompletedProcess(argv, 128, "", err)
            return subprocess.CompletedProcess(argv, 0, "", "")
        if argv[:2] == ["git", "remote"]:
            return subprocess.CompletedProcess(argv, 0, "origin\n", "")
        if argv[:3] == ["gh", "auth", "status"]:
            return subprocess.CompletedProcess(argv, 0 if logged_in else 1,
                                               "", "")
        if not git:
            raise FileNotFoundError("no git on PATH")
        return subprocess.CompletedProcess(argv, 0, "", "")

    return run


def _rows(monkeypatch, **kw):
    """Run ``ready_for_push`` with the host's tools shadowed (no network)."""
    which = {"git": None, "gh": None}
    if kw.get("git", True):
        which["git"] = "/usr/bin/git"
    if kw.get("gh", True):
        which["gh"] = "/usr/bin/gh"
    monkeypatch.setattr(shutil, "which", lambda name: which.get(name))
    monkeypatch.setattr(D.subprocess, "run", _fake_run(**kw))
    return {r["check"]: r for r in D.ready_for_push(".")}


def _by(monkeypatch, **kw):
    return _rows(monkeypatch, **kw)


# -- the five readiness questions are asked and answer per gate ------------

def test_five_readiness_questions_are_asked(monkeypatch):
    rows = _by(monkeypatch)
    assert set(rows) == {"git installed", "git identity", "git remote",
                         "gh installed", "gh authenticated"}


def test_everything_ready_means_every_row_passes(monkeypatch):
    rows = _by(monkeypatch)
    assert all(r["status"] == D.PASS for r in rows.values()), rows


def test_every_row_reads_the_standard_contract(monkeypatch):
    """U1-interoperable row format the gate and /fix both read."""
    for row in _by(monkeypatch).values():
        assert set(row) == {"check", "status", "detail", "fix"}


def test_no_readiness_row_declares_a_command(monkeypatch):
    """Installing git/gh and the login are the user's: prose fix only.

    This is the ``_row`` contract (doctor.py:54-74): a standard doctor row
    that has a machine-actionable remedy DECLARES its ``command``; these
    rows deliberately do not, because installing a system package and
    handing over the user's login are not things the agent improvises.
    """
    missing = _by(monkeypatch, git=False, identity=False, gh=False,
                  logged_in=False)
    assert not any(r.get("command") for r in missing.values())


# -- each missing piece is named with its own fix ---------------------------

def test_git_missing_is_a_warning_with_an_actionable_fix(monkeypatch):
    row = _by(monkeypatch, git=False)["git installed"]
    assert row["status"] != D.PASS
    assert row["fix"], "installing git needs an actionable remedy"
    assert "git" in row["fix"].lower()


def test_git_missing_still_asks_about_gh(monkeypatch):
    """The questions fail independently: no git must not hide a gh check."""
    rows = _by(monkeypatch, git=False, gh=False)
    assert "git installed" in rows
    assert "gh installed" in rows


def test_identity_missing_is_a_warning_with_a_prose_fix_only(monkeypatch):
    row = _by(monkeypatch, identity=False)["git identity"]
    assert row["status"] != D.PASS
    assert "user.name" in row["fix"]


def test_remote_unreachable_is_a_warning_with_a_diagnosis(monkeypatch):
    row = _by(monkeypatch, remote=False)["git remote"]
    assert row["status"] != D.PASS
    assert row["fix"], "the remote does not answer; a push needs a remedy"


def test_remote_without_origin_is_a_warning(monkeypatch):
    row = _by(monkeypatch, has_origin=False)["git remote"]
    assert row["status"] != D.PASS
    assert "origin" in (row["fix"] + row["detail"]).lower()


def test_gh_missing_is_a_warning_with_an_actionable_fix(monkeypatch):
    row = _by(monkeypatch, gh=False)["gh installed"]
    assert row["status"] != D.PASS
    assert row["fix"], "installing gh needs an actionable remedy"
    assert "gh" in row["fix"].lower()


def test_gh_logged_out_is_a_warning_with_a_login_fix(monkeypatch):
    """The login is the user's (`! gh auth login`), never an agent command."""
    row = _by(monkeypatch, logged_in=False)["gh authenticated"]
    assert row["status"] != D.PASS
    assert "login" in row["fix"]
    assert not row.get("command"), "the login is the user's, not an agent action"


# -- the runner never raises; a broken probe still answers ----------------

def test_a_runner_that_raises_still_yields_rows_never_an_exception(monkeypatch):
    """A host where nothing can even be asked must still answer, per gate.

    No exception propagates and no single broken probe hides the other
    readiness questions. ``git config`` / ``git ls-remote`` / ``gh auth
    status`` are the only probes that touch subprocess here; the runner
    raises OSError on every one, so the subprocess-touching rows degrade
    to a WARN with a fix. Presence is ``shutil.which`` -- found here --
    so ``git installed`` / ``gh installed`` stay PASS; that is the real
    shape of a host whose binaries are on PATH but whose configuration
    cannot be asked.
    """
    def boom(argv, **kw):
        raise OSError("no exec for anything")

    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/" + name)
    monkeypatch.setattr(D.subprocess, "run", boom)
    rows = {r["check"]: r for r in D.ready_for_push(".")}
    assert {c: rows[c]["status"] for c in rows} == {
        "git installed": D.PASS,
        "git identity": D.WARN,
        "git remote": D.WARN,
        "gh installed": D.PASS,
        "gh authenticated": D.WARN,
    }
    assert all(rows[c]["fix"] for c in ("git identity", "git remote",
                                        "gh authenticated")), \
        "every subprocess failure names a remedy"


def test_the_function_is_exported_for_the_gate():
    __import__("delfin.agent.doctor")
    assert hasattr(D, "ready_for_push")


def test_a_missing_remote_even_when_everything_else_is_present(monkeypatch):
    """One red row is enough to refuse, and the others still split out."""
    rows = _by(monkeypatch, remote=False)
    assert rows["git remote"]["status"] != D.PASS
    assert rows["git identity"]["status"] == D.PASS
    assert rows["gh authenticated"]["status"] == D.PASS


def test_gh_not_logged_in_but_git_fully_ready(monkeypatch):
    """gh login is independent of the git side; both rows answer."""
    rows = _by(monkeypatch, logged_in=False)
    assert rows["gh authenticated"]["status"] != D.PASS
    assert rows["git remote"]["status"] == D.PASS
    assert "login" in rows["gh authenticated"]["fix"]
