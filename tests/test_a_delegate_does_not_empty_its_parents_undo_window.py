"""A delegate's turn must not empty the parent's undo window.

The per-session ``ToolContext`` moved the tool executor's per-turn state
off the process-wide singleton, so two sessions in one process stopped
sharing an undo window. ``run_subagent`` runs the child against a SHALLOW
copy of the parent client (``copy.copy``), which shares every attribute
the copy does not reset -- and the reset list named only the two evidence
ledgers. So ``tool_context.for_client(child)`` returned the PARENT's
context object, the child's turn start called ``begin_turn()`` on it, and
the parent's undo window was emptied mid-turn while ``live_messages`` was
rebound to the child's conversation.

What that costs: a parent that delegated anything and then called
``undo_changes(scope="turn")`` reverted the delegate's files and left its
own edits on disk -- and reported success either way. The sharing the
context was built to end, surviving it by one missing line.

Measured by probe before the fix: ``for_client(parent) is
for_client(copy.copy(parent))`` was True, and the copy's ``begin_turn()``
cleared the parent's seq list.
"""

from __future__ import annotations

import copy

from delfin.agent import api_client as A
from delfin.agent import tool_context as TC


def _client():
    """A client without __init__, which is how the suite builds one (the
    real one opens a backend)."""
    return A.OpenAIClient.__new__(A.OpenAIClient)


def _child_of(parent):
    """The child as PRODUCTION builds it: a shallow copy, then the real
    isolation function. Re-implementing that list here would make the test
    pass whatever production forgot -- which is the defect itself."""
    from delfin.agent.subagents import _isolate_sub_client

    child = copy.copy(parent)
    _isolate_sub_client(child)
    return child


def test_run_subagent_isolates_the_copy_it_makes():
    """The copy and the isolation must stay together: a shallow copy that
    reaches the child unisolated is the defect."""
    import inspect

    from delfin.agent import subagents

    src = inspect.getsource(subagents.run_subagent)
    i = src.index("_copy.copy(parent_client)")
    assert "_isolate_sub_client(sub_client)" in src[i:i + 600], (
        "the shallow copy is no longer isolated at the call site")


def test_a_childs_turn_start_leaves_the_parents_window_alone():
    parent = _client()
    parent_ctx = TC.for_client(parent)
    parent_ctx._turn_change_seqs = [11, 12]

    child = _child_of(parent)
    child_ctx = TC.for_client(child)
    assert child_ctx is not parent_ctx, "the child shares the parent's context"

    child_ctx.begin_turn()
    assert parent_ctx._turn_change_seqs == [11, 12], (
        "the child's turn start emptied the parent's undo window")


def test_the_childs_conversation_does_not_become_the_parents():
    """``live_messages`` is what history_search resolves live refs against;
    on a shared context the parent's refs resolved to the delegate's
    history."""
    parent = _client()
    TC.for_client(parent).live_messages = [{"role": "user", "content": "parent"}]
    child = _child_of(parent)
    TC.for_client(child).live_messages = [{"role": "user", "content": "child"}]
    said = [m["content"] for m in TC.for_client(parent).live_messages]
    assert said == ["parent"], said
