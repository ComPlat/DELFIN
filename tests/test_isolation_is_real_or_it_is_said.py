"""Worktree isolation was voided in silence, and skipped where it is needed.

Two halves of the same promise.

**Voided in silence.** ``run_subagent`` creates the worktree first, then
``_derive_perms`` drops the workspace override when the session is scope
locked -- correctly, because a lock the child could relocate is not a
lock. But nothing said so. The worktree was created and torn down for
nothing, the child ran in the parent's workspace, and the result payload
still reported ``isolation: "worktree"`` with a worktree summary. The
caller asked for isolation, was given none, and was told it had it.

**Skipped where it is needed most.** Auto-isolation for writer presets
lives in the parallel fan-out and triggers only at two or more subagent
calls in one turn. A ``background=True`` writer is a single call, so it
got none -- and it is the worst case of all: fan-out siblings at least
finish before the turn ends, while a background writer edits the tree
while the parent is editing it. Same tree, no worktree, no lock, no
stale-write baseline between them.

The rule is the same in both: isolation is either real, or it is said.
"""

from __future__ import annotations

import json
import types

import pytest

from delfin.agent import subagents as sa


def _Perms(workspace, *, scope_locked=False):
    """The real permissions type: _derive_perms clones with
    dataclasses.replace, so a stand-in silently takes a fallback path."""
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=workspace, mode="acceptEdits")
    if scope_locked:
        perms.lock_workspace = True
    return perms


# ---------------------------------------------------------------------------
# A lock voids the override -- and now says so
# ---------------------------------------------------------------------------

def test_a_locked_session_does_not_move_the_child(tmp_path):
    """The rule itself is right and must not regress: a lock the child can
    relocate is not a lock."""
    perms = _Perms(tmp_path, scope_locked=True)
    derived = sa._derive_perms(perms, "acceptEdits", workspace=tmp_path / "wt")
    assert derived.workspace == tmp_path


def test_an_unlocked_session_still_moves_the_child(tmp_path):
    perms = _Perms(tmp_path)
    derived = sa._derive_perms(perms, "acceptEdits", workspace=tmp_path / "wt")
    assert derived.workspace == tmp_path / "wt"


def test_a_locked_session_does_not_create_a_worktree(tmp_path, monkeypatch):
    """It was created, unused, then torn down."""
    called = []
    monkeypatch.setattr(sa, "_worktree_would_be_voided",
                        sa._worktree_would_be_voided)
    assert sa._worktree_would_be_voided(_Perms(tmp_path, scope_locked=True))
    assert not sa._worktree_would_be_voided(_Perms(tmp_path))
    assert not called


def test_the_reason_is_stated_not_implied():
    note = sa._ISOLATION_VOIDED_BY_LOCK.lower()
    assert "confined to one folder" in note, "it must name the reason"
    assert "not isolated" in note, "and state the consequence"


# ---------------------------------------------------------------------------
# A background writer is isolated like a parallel one
# ---------------------------------------------------------------------------

def test_a_background_writer_is_auto_isolated():
    assert sa.auto_isolation_for("general-purpose", "", background=True) == \
        "worktree"


def test_a_background_reader_needs_no_worktree():
    assert sa.auto_isolation_for("explore", "", background=True) == ""


def test_an_explicit_choice_is_respected():
    assert sa.auto_isolation_for("general-purpose", "none",
                                 background=True) == "none"


def test_a_foreground_single_writer_is_left_alone():
    """One writer in the foreground has the tree to itself -- a worktree
    would only add a merge step the user did not ask for."""
    assert sa.auto_isolation_for("general-purpose", "", background=False) == ""


def test_an_unknown_preset_is_not_assumed_to_write():
    assert sa.auto_isolation_for("no-such-preset", "", background=True) == ""


# ---------------------------------------------------------------------------
# The tool layer uses the rule
# ---------------------------------------------------------------------------

# The next three tests were one source-text assertion: find the index of
# ``if bool(arguments.get("background")):`` and assert the name
# ``auto_isolation_for`` occurs within 2000 CHARACTERS of it. That is
# satisfied by a comment, by an unrelated function, or by a call on a
# neighbouring branch -- and it breaks when legitimate code is inserted
# between the two, which is a test that fails for being right.
#
# What the promise actually is: a background writer is LAUNCHED with
# worktree isolation. So launch one and read what the runner was handed.

def _launched_isolation(subagent_type, *, isolation="", timeout=3.0):
    """Spawn a background sub-agent and return the isolation it was
    launched with."""
    import threading

    from delfin.agent import api_client as A

    received: dict = {}
    done = threading.Event()

    def _runner(**kw):
        received.update(kw)
        done.set()
        return {"ok": True}

    class _Perms:
        subagent_runner = staticmethod(_runner)

    args = {"subagent_type": subagent_type,
            "description": "edit the parser module",
            "prompt": "rewrite the tokenizer so it handles quoted strings",
            "background": True}
    if isolation:
        args["isolation"] = isolation

    out = json.loads(A._doc_executor._execute_subagent(args, _Perms()))
    assert out.get("status") == "started_in_background", out
    assert done.wait(timeout=timeout), "the background runner never ran"
    return received.get("isolation")


def test_a_background_writer_is_launched_in_a_worktree():
    """The case the rule exists for: it edits the tree while the parent
    is editing it."""
    assert _launched_isolation("general-purpose") == "worktree"


def test_a_background_reader_is_launched_without_one():
    """A worktree for a reader is a merge step nobody asked for."""
    assert not _launched_isolation("explore")


def test_an_explicit_choice_survives_the_background_dispatch():
    """The rule fills a gap; it must not overrule the caller."""
    assert _launched_isolation("general-purpose", isolation="none") == "none"


def test_the_parallel_fan_out_still_isolates_writers():
    """The behaviour that already existed must survive the refactor."""
    assert sa.is_writer_preset("general-purpose")
    assert not sa.is_writer_preset("explore")
