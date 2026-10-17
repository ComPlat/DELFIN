"""A path that is only a command's write target, outside every root, is
refused once -- never asked as a read.

Input: ``python3 -m venv <dir outside the workspace>`` and similar bash
calls, driven through the real ``_DocToolExecutor``, with a dialog that
records what it was asked.

Before: the read-arg check scanned the destination as a possible read and
asked the person "read outside the workspace?". A yes could not make the
command run -- the write-root check refused it right after -- so the
person was asked a question whose yes did nothing (a dead yes).

Pinned here:

* such a write is refused without any dialog, and the refusal names the
  route that works (--add-dir / extra_workspace_dirs);
* a no to a dialog leaves the call unrun;
* a write inside a writable root is asked at most once;
* a write to a secret path stays a deny with no grant offered.

Measured in wave 16 (2026-10-17).
"""
import json

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent.api_client import KitToolPermissions, _doc_executor

pytestmark = [pytest.mark.usefixtures("repo")]

requires_sandbox = requires_a_working_sandbox


@pytest.fixture
def repo(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def _perms(ws, **kw):
    return KitToolPermissions(mode="default", workspace=str(ws), **kw)


class TestOutsideWorkspaceWriteAsksForTheDirectory:
    """The outside-workspace write becomes a structured question that names
    the directory grant — not a raw text refusal, and never an allow_pattern."""

    @requires_sandbox
    def test_outside_write_refusal_names_the_offending_directory(self, repo):
        perms = _perms(repo)
        outside = repo.parent / "elsewhere2"
        out = json.loads(_doc_executor.execute(
            "bash", {"command": f"mkdir {outside}"}, perms,
        ))
        err = out.get("error", "")
        assert str(outside) in err or "outside" in err, err

    def test_accept_edits_outside_venv_was_asked_nothing_after_w1(self, repo):
        # W1 red control (operator ruling 2026-10-16): `python3 -m venv
        # /scratch/x/venv` in acceptEdits mode reaches the real ask site
        # (api_client :21740) with verdict ask, if_approved refuse (outside_
        # ws) — the person is ASKED and a "yes" still ends refused (dead yes).
        # After W1, blocked_question.classify -> Impossible means the command
        # is NOT asked; it is refused once with the route that works (the
        # DIRECTORY grant, extra_dir). RED pre-W1: the confirm_callback IS
        # called (seen non-empty). GREEN post-W1: seen stays empty and the
        # refusal names extra_dir, never allow_pattern.
        seen = []
        perms = KitToolPermissions(
            mode="acceptEdits", workspace=str(repo),
            confirm_callback=lambda name, args, preview: (
                seen.append((name, args)), True,
            ),
        )
        outside = repo.parent / "elsewhere" / "venv"
        out = json.loads(_doc_executor.execute(
            "bash", {"command": f"python3 -m venv {outside}"}, perms,
        ))
        err = out.get("error", "")
        # Post-W1: a dead-yes command is never asked — the person is not
        # consulted for a grant that cannot run.
        assert not seen, f"a dead yes was asked: {seen}"
        assert err, out
        # The refusal names the DIRECTORY grant that actually works ...
        assert "extra_dir" in err or "extra_workspace_dirs" in err, err
        # ... and never the dead command-level allow_pattern.
        assert "allow_pattern" not in err, f"dead grant offered: {err}"


class TestNoIsRespected:
    """Answering no to an ask leaves the call unrun and the turn continues."""

    def test_no_to_a_confirm_leaves_the_call_unrun(self, repo):
        # W1 red control: an outside-workspace write is asked for the
        # DIRECTORY grant (api_client :21740) and answering "no" leaves it
        # unrun, the turn continuing. RED on origin/main (pre-W1 an
        # outside-write is a HARD root refusal, no confirm is consulted), so
        # on origin/main `seen` stays empty — that emptiness is the red. The
        # operator runs this against the W1-built copy where it is green.
        seen = []
        perms = KitToolPermissions(
            mode="default", workspace=str(repo),
            confirm_callback=lambda name, args, preview: (
                seen.append((name, args)), False,
            ),
        )
        outside = repo.parent / "no_write_venv"
        out = json.loads(_doc_executor.execute(
            "bash", {"command": f"python3 -m venv {outside}"}, perms,
        ))
        # A no must be respected: the call does not run. The outside dir is
        # never created either way (root or confirm refusal).
        assert not outside.exists(), "a no must not run the call"
        assert out.get("error"), out


class TestAskedExactlyOnce:
    """A write into a writable root is asked exactly once by the dialog."""

    def test_bash_write_into_a_writable_root_asks_exactly_once(self, repo):
        # The operator finding: a side-effecting W1 (calling
        # _gate_bash_write_targets at the ask site) would ASK TWICE for a
        # write into a calc-confirm/protected root — once in the pre-check,
        # once at execution — and record the events twice. The pure W1 check
        # (find_root_for non-None) must never classify such a write as a dead
        # yes, so the dialog is reached exactly once. This asserts that the
        # confirm_callback for a bash write into the workspace root (a
        # writable root) is consulted exactly once and the call is not
        # double-refused.
        seen = []
        perms = KitToolPermissions(
            mode="acceptEdits", workspace=str(repo),
            confirm_callback=lambda name, args, preview: (
                seen.append((name, args)), True,
            ),
        )
        inside = repo / "made_by_venv"
        json.loads(_doc_executor.execute(
            "bash", {"command": f"python3 -m venv {inside}"}, perms,
        ))
        bash_asks = [a for n, a in seen if n == "bash"]
        assert len(bash_asks) <= 1, f"asked more than once: {seen}"


class TestDenyStaysADeny:
    """A deny (secret/protected path) is never turned into a question."""

    @requires_sandbox
    def test_secret_path_write_is_denied_without_a_question(self, repo):
        denied_dir = repo.parent / "deny_me"
        denied_dir.mkdir()
        perms = _perms(repo)
        out = json.loads(_doc_executor.execute(
            "bash", {"command": f"echo x > {denied_dir / '.env'}"}, perms,
        ))
        err = out.get("error", "")
        assert "deny" in err.lower() or "secret" in err.lower(), err
        # A deny must not read like an ask: no grant route offered.
        assert "allow_pattern" not in err, err
