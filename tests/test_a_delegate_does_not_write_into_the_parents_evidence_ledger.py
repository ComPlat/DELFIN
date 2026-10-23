"""A delegate's reads landed straight in the PARENT's evidence ledger.

A sub-agent runs on ``copy.copy(parent_client)``. Every turn the client
starts a fresh per-turn ledger and, guarded by ``hasattr``, a session
ledger::

    self._observed_files = set()
    if not hasattr(self, "_observed_files_session"):
        self._observed_files_session = set()

A shallow copy already HAS the second attribute -- bound to the PARENT's
set object. The guard therefore skipped the rebind, and the union run
after every tool call (``self._observed_files_session |= ...``) mutated
the parent's set in place.

That set is what the parent's ungrounded-citation guard reads and what
the task-completion check reads. So a child greps engine.py, the parent
never opens it, the parent writes "engine.py:412 does X", and the guard
stays silent. It survives the child's report being DISCARDED: on the
fan-out timeout path the parent records "did not finish in time" and
moves on while the abandoned thread keeps unioning into the ledger, so a
run the parent explicitly rejected still disarms the parent's honesty
guard.

``_merge_delegate_evidence`` is the intended channel -- narrow, derived
from the child's own tool trace, and applied only to a report the parent
accepted. This was a second, unfiltered one running beside it.
"""

from __future__ import annotations

import copy

import pytest

from delfin.agent import api_client as ac
from delfin.agent import subagents as sa


@pytest.fixture(autouse=True)
def _iso(monkeypatch, tmp_path):
    monkeypatch.setattr(sa, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")


class _Event:
    def __init__(self, text: str) -> None:
        self.type = "text_delta"
        self.text = text
        self.tool_name = ""
        self.tool_input = ""
        self.tool_output = ""
        self.input_tokens = 0
        self.output_tokens = 0


class _Client:
    """The client shape ``run_subagent`` copies: a real (shallow-copyable)
    object carrying both evidence ledgers."""

    def __init__(self) -> None:
        self.model = "parent-model"
        self._permissions = None
        self._observed_files: set[str] = set()
        self._observed_files_session: set[str] = {"parent/opened.py"}
        self.ran: list["_Client"] = []

    def set_permissions(self, perms) -> None:
        self._permissions = perms

    def stream_message(self, messages=None, system="", max_tokens=0):
        # What the real tool loop does after every tool result.
        self._observed_files = {"child/grepped.py"}
        self._observed_files_session |= self._observed_files
        self.ran.append(self)
        return iter([_Event("looked around")])


def _run(parent: _Client):
    return sa.run_subagent(
        subagent_type="explore",
        description="find the caller",
        prompt="Find every caller of the thing, carefully.",
        parent_client=parent,
        parent_perms=None,
    )


# ---------------------------------------------------------------------------
# The delegate cannot reach the parent's ledger
# ---------------------------------------------------------------------------

def test_a_delegates_reads_do_not_reach_the_parent_ledger():
    parent = _Client()
    _run(parent)
    assert "child/grepped.py" not in parent._observed_files_session


def test_the_parents_own_ledger_survives_the_delegate():
    parent = _Client()
    _run(parent)
    assert parent._observed_files_session == {"parent/opened.py"}


def test_the_child_ledger_is_a_different_object():
    parent = _Client()
    _run(parent)
    child = parent.ran[0]
    assert child is not parent
    assert child._observed_files_session is not parent._observed_files_session


def test_the_child_still_keeps_its_own_evidence():
    """Isolation, not amnesia: the child's own ledger is intact, which is
    what the verified merge channel is derived from."""
    parent = _Client()
    _run(parent)
    assert "child/grepped.py" in parent.ran[0]._observed_files_session


# ---------------------------------------------------------------------------
# ...and the client refuses to adopt a ledger it does not own
# ---------------------------------------------------------------------------

class _Bare:
    """Anything shallow-copyable; the ledger init only touches attributes."""


def test_a_shallow_copy_starts_its_own_session_ledger():
    parent = _Bare()
    ac._begin_observed_ledgers(parent)
    parent._observed_files_session.add("parent/opened.py")

    child = copy.copy(parent)
    ac._begin_observed_ledgers(child)
    child._observed_files_session.add("child/grepped.py")

    assert parent._observed_files_session == {"parent/opened.py"}
    assert child._observed_files_session == {"child/grepped.py"}


def test_the_owning_client_accumulates_across_turns():
    """The session ledger is cumulative for the client that owns it."""
    client = _Bare()
    ac._begin_observed_ledgers(client)
    client._observed_files_session.add("turn1.py")
    ac._begin_observed_ledgers(client)
    assert client._observed_files_session == {"turn1.py"}


def test_the_per_turn_ledger_is_always_fresh():
    client = _Bare()
    ac._begin_observed_ledgers(client)
    client._observed_files.add("turn1.py")
    ac._begin_observed_ledgers(client)
    assert client._observed_files == set()
