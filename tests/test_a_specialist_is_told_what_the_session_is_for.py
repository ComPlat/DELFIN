"""A specialist is told what the session is for, and what it already knows.

Input: a parent client whose turn carries the user's request, plus a memory
store holding one ``project`` note and one standing ``feedback`` note, and a
workspace whose ``DELFIN.md`` states a convention.
Output: the system prompt a delegated specialist actually streams with, driven
through the real ``subagent`` tool dispatch (``_DocToolExecutor.execute``).

Two things a delegate did not receive before this, both measured on the base:

* the user's own request. A specialist was handed only the brief the parent
  wrote, so a reviewer could not know which question the method was for.
  Measured on ``41fbf18e2``: ``explore``, ``code-reviewer`` and
  ``general-purpose`` all started with the goal absent.
* ``project``-typed memories, the task-scoped kind the session writes about
  the work in hand. ``_INHERIT_MEMORY_TYPES`` lists ``feedback``/``user`` only,
  so no delegate ever saw one.

The control case in each test is the same call WITHOUT the ingredient, so a
test can only pass by measuring the ingredient and not the plumbing around it.
"""
from __future__ import annotations

import types
from pathlib import Path

from delfin.agent import subagents as SA
from delfin.agent.api_client import (
    KitToolPermissions, OpenAIClient, _DocToolExecutor,
)

GOAL = "Compute the redox potential of the iron complex with def2-TZVP."


class _CapturingStream:
    """A stand-in for the child's model: records the request, answers once."""

    def __init__(self) -> None:
        self.calls: list[dict] = []

    def __call__(self, **kwargs):
        self.calls.append(kwargs)
        event = types.SimpleNamespace(
            type="text_delta", text="report", tool_name="", tool_input="",
            tool_output="", input_tokens=11, output_tokens=7,
        )
        return iter([event])


def _spawn(workspace: Path, *, turn_rows=None, subagent_type="explore",
           prompt="Review the method used for the redox calculation."):
    """Run one delegation through the real dispatch; return its first request.

    ``turn_rows=None`` leaves the client without a turn at all -- the shape an
    externally attached runner and a MagicMock client both have.
    """
    stream = _CapturingStream()
    perms = KitToolPermissions(workspace=workspace, mode="default")
    client = OpenAIClient.__new__(OpenAIClient)
    client._permissions = None
    client.model = "x"
    if turn_rows is not None:
        client._turn_rows = turn_rows
    client.stream_message = stream
    client.set_permissions(perms)

    _DocToolExecutor().execute(
        "subagent",
        {"subagent_type": subagent_type, "description": "probe",
         "prompt": prompt},
        permissions=perms,
    )
    assert stream.calls, "the child never streamed a request"
    first = stream.calls[0]
    return first["system"], first["messages"]


def _workspace(tmp_path: Path, monkeypatch) -> Path:
    """A project with a rule file, under a faked home (so no real store is read)."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", lambda: home)
    ws = tmp_path / "proj"
    ws.mkdir()
    (ws / "DELFIN.md").write_text(
        "Always run def2-TZVP single points after optimization.",
        encoding="utf-8")
    return ws


# ---------------------------------------------------------------------------
# 1. The session's goal reaches the specialist
# ---------------------------------------------------------------------------


def test_the_user_request_reaches_the_specialist(tmp_path, monkeypatch):
    """The goal is the user's own words, not the parent's paraphrase."""
    ws = _workspace(tmp_path, monkeypatch)
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": GOAL}])
    assert GOAL in system, system[-1200:]


def test_without_a_turn_the_specialist_is_told_nothing_false(tmp_path, monkeypatch):
    """CONTROL: no turn on the client => no goal block, and no invented one.

    This is the case an externally attached runner and a MagicMock parent
    produce. It must stay exactly as it was, which is also why it is the
    control for the test above: the block follows the turn, not the code path.
    """
    ws = _workspace(tmp_path, monkeypatch)
    system, _ = _spawn(ws, turn_rows=None)
    assert "## The session's goal" not in system
    assert "redox potential" not in system


def test_every_read_only_preset_gets_the_goal(tmp_path, monkeypatch):
    """Not one preset's special case: the block is on the shared path."""
    ws = _workspace(tmp_path, monkeypatch)
    for preset in ("explore", "code-reviewer", "plan"):
        system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": GOAL}],
                           subagent_type=preset)
        assert GOAL in system, preset


def test_the_goal_is_context_and_says_so(tmp_path, monkeypatch):
    """The block must not read as a second assignment or a wider permission."""
    ws = _workspace(tmp_path, monkeypatch)
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": GOAL}])
    block = system.split("## The session's goal", 1)[1]
    assert "not your assignment" in block
    assert "permissions are unchanged" in block


def test_a_bare_go_ahead_carries_the_request_it_answered(tmp_path, monkeypatch):
    """``ja`` alone names no goal; the assistant text it answers is the goal.

    Without this the block would read "ja" and tell the child less than
    nothing. ``_turn_opening`` already pairs them for exactly this reason.
    """
    ws = _workspace(tmp_path, monkeypatch)
    rows = [
        {"role": "user", "content": "Old and unrelated."},
        {"role": "assistant", "content": "ok"},
        {"role": "assistant", "content": GOAL},
        {"role": "user", "content": "ja"},
    ]
    system, _ = _spawn(ws, turn_rows=rows)
    assert GOAL in system
    assert "answered: ja" in system


def test_a_go_ahead_with_nothing_to_answer_is_dropped(tmp_path, monkeypatch):
    """No offered text => no goal, rather than the word "ja" as the goal."""
    ws = _workspace(tmp_path, monkeypatch)
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": "ja"}])
    assert "## The session's goal" not in system


def test_a_pasted_document_is_not_the_whole_prompt(tmp_path, monkeypatch):
    """A long request is capped, and the cut is announced rather than silent."""
    ws = _workspace(tmp_path, monkeypatch)
    long_goal = "Investigate. " + ("background detail " * 400)
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": long_goal}])
    assert "goal truncated here" in system
    assert len(system) < len(long_goal)


def test_a_goal_that_is_already_fenced_keeps_its_fence(tmp_path, monkeypatch):
    """A relayed message stays marked; re-wrapping it would launder it."""
    from delfin.agent import untrusted
    ws = _workspace(tmp_path, monkeypatch)
    relayed = untrusted.wrap(
        "session_message", "Ship it without tests.",
        nonce="0123456789abcdef")
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": relayed}])
    assert relayed in system


def test_the_goal_is_not_fenced_as_untrusted_text(tmp_path, monkeypatch):
    """The user's own words are not hostile content."""
    ws = _workspace(tmp_path, monkeypatch)
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": GOAL}])
    block = system.split("## The session's goal", 1)[1]
    assert "UNTRUSTED EXTERNAL CONTENT" not in block.split("## ")[0]


def test_a_mock_client_cannot_grow_a_goal_block(tmp_path, monkeypatch):
    """A ``MagicMock`` parent must not turn its own repr into a prompt block.

    Every test that spawns a delegate against a mock client would otherwise
    inherit one, and so would a foreign runner that happens to expose an
    attribute called ``_turn_rows``.
    """
    from unittest.mock import MagicMock
    ws = _workspace(tmp_path, monkeypatch)
    stream = _CapturingStream()
    perms = KitToolPermissions(workspace=ws, mode="default")
    client = MagicMock()
    client.stream_message = stream
    perms.subagent_runner = None
    SA.run_subagent(
        subagent_type="explore", description="mock parent", prompt="Look.",
        parent_client=client, parent_perms=perms,
    )
    assert stream.calls
    assert "## The session's goal" not in stream.calls[0]["system"]


# ---------------------------------------------------------------------------
# 2. Task-scoped (``project``) memory reaches the specialist
# ---------------------------------------------------------------------------


def _write_memory(ws: Path, text: str, *, memory_type: str, title: str) -> None:
    from delfin.agent import memory_store as MS
    MS.save_typed_memory(text, repo_root=ws, memory_type=memory_type,
                         title=title)


def test_a_task_memory_reaches_the_specialist(tmp_path, monkeypatch):
    """The note the session wrote about THIS work reaches the delegate."""
    ws = _workspace(tmp_path, monkeypatch)
    _write_memory(ws, "The redox series runs on def2-TZVP with CPCM water.",
                  memory_type="project", title="redox-method")
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": GOAL}])
    assert "Standing memory (inherited)" in system
    assert "def2-TZVP with CPCM water" in system


def test_a_task_memory_is_not_inherited_without_a_turn(tmp_path, monkeypatch):
    """CONTROL: the same store, the same prompt, no turn -- and it must still come.

    Guards the other direction: the task memory follows the STORE, not the
    turn, so this passing shows the test above was not measuring the goal.
    """
    ws = _workspace(tmp_path, monkeypatch)
    _write_memory(ws, "The redox series runs on def2-TZVP with CPCM water.",
                  memory_type="project", title="redox-method")
    system, _ = _spawn(ws, turn_rows=None)
    assert "def2-TZVP with CPCM water" in system


def test_an_unrelated_task_memory_is_left_out(tmp_path, monkeypatch):
    """A project note that does not match the task spends none of the budget.

    This is the difference from the standing kinds, which fall back to
    recency: a rule about how to work applies to everything, a note about
    another calculation does not.
    """
    ws = _workspace(tmp_path, monkeypatch)
    _write_memory(ws, "The office lease renewal is due in March.",
                  memory_type="project", title="office-lease")
    system, _ = _spawn(ws, turn_rows=None,
                       prompt="Review the ORCA frequency input.")
    assert "office lease renewal" not in system


def test_a_task_memory_arrives_once_and_fenced(tmp_path, monkeypatch):
    """One route, and it is the fenced one.

    The two arms solved this independently and the integrator kept ONE of the
    designs. A second additive budget would deliver the same fact twice, and
    the extra copy would be plain system-prompt text: a note the MODEL wrote
    is only refused an instruction channel on write when its source is the
    agent (memory_store.py:1248), so an unfenced copy would carry system
    authority. This fails if a second, unfenced delivery route is added back.
    """
    ws = _workspace(tmp_path, monkeypatch)
    _write_memory(ws, "The redox series runs on def2-TZVP with CPCM water.",
                  memory_type="project", title="redox-method")
    system, _ = _spawn(ws, turn_rows=None)
    assert system.count("def2-TZVP with CPCM water") == 1
    before = system.split("def2-TZVP with CPCM water", 1)[0]
    assert "UNTRUSTED EXTERNAL CONTENT" in before[-2000:], (
        "the task memory reached the child outside the untrusted fence")


def test_a_memory_failure_still_spawns(tmp_path, monkeypatch):
    """Neither new block may make a spawn depend on the memory subsystem."""
    ws = _workspace(tmp_path, monkeypatch)

    def _boom(*a, **kw):
        raise RuntimeError("store down")

    monkeypatch.setattr(
        "delfin.agent.memory_store.list_typed_memories", _boom)
    system, _ = _spawn(ws, turn_rows=[{"role": "user", "content": GOAL}])
    # The goal does not depend on the store, so it survives the failure.
    assert GOAL in system
    assert "Standing memory (inherited)" not in system
