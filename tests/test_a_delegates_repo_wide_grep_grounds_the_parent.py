"""A delegate sent out to FIND something grounded the parent with nothing.

The parent's own harvest treats a repo-wide grep as a read: the call
declares no path, so the hits are taken from the OUTPUT
(``_observe_read_files`` → ``_paths_in_grep_output``). The child's harvest
never got that branch — it reads ``_declared_paths(args)`` and stops. A
delegate whose whole job is to locate something therefore contributed
nothing to the parent's ledger, and a finding it relayed read as
unsourced.

MEASURED before touching it, over the 99 recorded delegate sessions in
~/.delfin/subagent_sessions:

    grep_file calls by delegates      107
      with a path argument             85
      with NO path argument            22   (21%)
    sessions containing at least one   14   of 99  (14%)

For comparison, across all calls — parent and child together — only 10.6%
of greps are repo-wide. Delegates do it twice as often, which is what one
would expect of an agent dispatched to find a thing nobody has located
yet. So in roughly one delegate session in seven the parent was penalised
for relaying a true finding, which makes delegating more expensive than
doing the work — the opposite of what the delegation reminder is for.

The rule #105 established is untouched and is asserted below: paths come
from the tool's RESULT, never from arguments or from a command the
delegate typed. A grep result is the tool's own output; a bash command
line is the delegate's prose.
"""

from __future__ import annotations

import pytest

from delfin.agent.subagents import collect_report_evidence


def _call(name, args, output):
    return {"name": name, "input": args, "output": output}


_HITS = (
    "delfin/agent/engine.py:412: def stream_response(self, user_message):\n"
    "delfin/agent/api_client.py:15479: _observe_read_files(\n"
)


def _files_read(calls):
    return collect_report_evidence({}, tool_calls=calls)["files_read"]


# ---------------------------------------------------------------------------
# The hole
# ---------------------------------------------------------------------------

def test_a_repo_wide_grep_grounds_the_files_it_showed():
    read = _files_read([_call("grep_file", {"pattern": "stream_response"},
                              _HITS)])
    assert "delfin/agent/engine.py" in read
    assert "delfin/agent/api_client.py" in read


def test_a_grep_with_a_path_still_grounds_that_path():
    read = _files_read([_call("grep_file",
                              {"path": "delfin/agent/engine.py",
                               "pattern": "x"},
                              "delfin/agent/engine.py:1: x\n")])
    assert read == ["delfin/agent/engine.py"]


def test_the_declared_path_wins_and_the_output_is_not_also_harvested():
    """Mirror of the parent: a declared path is what the tool acted on.
    Harvesting both would let one call ground half the repository."""
    read = _files_read([_call("grep_file",
                              {"path": "delfin/agent/engine.py",
                               "pattern": "x"},
                              _HITS)])
    assert read == ["delfin/agent/engine.py"]


def test_the_mcp_prefixed_name_is_recognised_too():
    read = _files_read([_call("mcp__delfin-docs__grep_file",
                              {"pattern": "stream_response"}, _HITS)])
    assert "delfin/agent/engine.py" in read


# ---------------------------------------------------------------------------
# What must still ground nothing
# ---------------------------------------------------------------------------

def test_a_grep_that_matched_nothing_grounds_nothing():
    assert _files_read([_call("grep_file", {"pattern": "zzz"}, "")]) == []


def test_a_grep_that_errored_grounds_nothing():
    """The error text has to LOOK like a hit, or the test proves nothing.

    My first version used a message with no ``file:line:`` in it, so the
    extractor found nothing either way and the guard could be deleted
    without a single test noticing. A real error message quotes the place
    it failed, which is exactly the shape of a hit.
    """
    err = ('{"error": "grep failed at delfin/agent/engine.py:412: '
           'unterminated subpattern"}')
    assert _files_read([_call("grep_file", {"pattern": "f("}, err)]) == []


def test_a_call_with_no_recorded_output_grounds_nothing():
    """A wall-clock cut leaves an input with no result behind. Counting it
    would make a truncated run ground more than a clean one."""
    assert _files_read([_call("grep_file", {"pattern": "x"}, None)]) == []


def test_a_path_the_delegate_only_typed_still_grounds_nothing():
    """#105, unchanged: a bash command names files in prose, and naming is
    not reading. Only a read/write tool's own result counts."""
    read = _files_read([_call(
        "bash", {"command": "grep -rn stream_response delfin/agent/engine.py"},
        "delfin/agent/engine.py:412: def stream_response(...)\n")])
    assert read == []


def test_a_report_that_merely_mentions_a_path_grounds_nothing():
    ev = collect_report_evidence(
        {"report": "I checked delfin/agent/engine.py and it is fine."},
        tool_calls=[])
    assert ev["files_read"] == []


# ---------------------------------------------------------------------------
# It reaches the parent
# ---------------------------------------------------------------------------

def test_the_parent_ledger_receives_the_hits():
    """End to end: child evidence → verification block → parent's set."""
    import json

    from delfin.agent.api_client import _merge_delegate_evidence

    read = _files_read([_call("grep_file", {"pattern": "stream_response"},
                              _HITS)])
    payload = json.dumps({"verification": {"evidence": {
        "files_touched": read}}})
    observed: set = set()
    _merge_delegate_evidence(observed, "subagent_result", payload)
    assert "delfin/agent/engine.py" in observed


def test_nothing_is_recorded_twice():
    read = _files_read([
        _call("grep_file", {"pattern": "a"}, _HITS),
        _call("grep_file", {"pattern": "b"}, _HITS),
    ])
    assert len(read) == len(set(read))


@pytest.mark.parametrize("junk", [None, 42, {"a": 1}, [1, 2]])
def test_a_result_that_is_not_text_does_not_raise(junk):
    assert _files_read([_call("grep_file", {"pattern": "x"}, junk)]) == []
