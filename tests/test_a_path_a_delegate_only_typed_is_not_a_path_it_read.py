"""Paths a delegate merely TYPED counted as paths it had observed.

``_paths_from_call`` falls back to joining every string argument of a call
and regex-scraping path-shaped tokens out of the blob. So::

    bash: echo "checked delfin/agent/engine.py and delfin/agent/hooks.py"

harvested two paths, which became ``files_touched`` in the verdict, which
the parent unions into its own evidence ledger. Neither file was opened
by anyone, and the parent's ungrounded-citation guard was disarmed for
both.

The second half is the same fault from the other end:
``_output_saw_content(None)`` returns True, so a tool_use whose result
never arrived -- exactly what a wall-clock cut leaves behind -- had its
input paths counted. The truncated run therefore grounded MORE than the
clean one, where a refused read is dropped.

The loose scrape has a legitimate job on the LENIENT side: it suppresses
false alarms in the child's own report check, where over-collecting can
only reduce flags. It must not be a source of grounding. The two
directions are separated here: ``files`` stays generous, ``files_read``
carries only what a read or write tool was actually pointed at and
actually returned, and only ``files_read`` is published as
``files_touched`` for the parent to merge.
"""

from __future__ import annotations

from delfin.agent import subagents as sa


def _ev(calls):
    return sa.collect_report_evidence({"final_text": "x", "tool_calls": calls})


def _touched(calls) -> list[str]:
    v = sa.verify_subagent_report({"final_text": "Report.", "tool_calls": calls})
    return list(v["evidence"]["files_touched"])


# ---------------------------------------------------------------------------
# Typing a path is not reading it
# ---------------------------------------------------------------------------

_ECHOED = {
    "name": "bash",
    "input": {"command": ('echo "checked delfin/agent/engine.py and '
                          'delfin/agent/hooks.py"')},
    "output": "checked delfin/agent/engine.py and delfin/agent/hooks.py",
}


def test_a_path_named_in_a_command_is_not_grounding_evidence():
    assert _ev([_ECHOED])["files_read"] == []


def test_a_typed_path_is_not_published_to_the_parent():
    assert _touched([_ECHOED]) == []


def test_the_loose_scrape_is_still_there_for_the_lenient_side():
    """Over-collecting only suppresses false alarms in the child's own
    report check; that job is unchanged."""
    files = _ev([_ECHOED])["files"]
    assert any("engine.py" in f for f in files)


def test_a_bash_read_does_not_ground_the_parent_either():
    """A path scraped out of a command string is a guess about what the
    command did, whichever command it was."""
    ev = _ev([{"name": "bash",
               "input": {"command": "sed -n 1,40p delfin/agent/engine.py"},
               "output": "def run(): ..."}])
    assert ev["files_read"] == []


# ---------------------------------------------------------------------------
# A read tool, pointed at a file, that returned something
# ---------------------------------------------------------------------------

def test_a_successful_read_grounds_the_parent():
    ev = _ev([{"name": "read_file", "input": {"path": "a.py"},
               "output": "1  import os"}])
    assert ev["files_read"] == ["a.py"]


def test_a_write_grounds_the_parent():
    ev = _ev([{"name": "write_file", "input": {"path": "b.py"},
               "output": "written"}])
    assert "b.py" in ev["files_read"]


def test_an_edit_record_grounds_the_parent():
    ev = _ev([{"name": "multi_edit",
               "input": {"edits": [{"path": "c.py", "old": "a", "new": "b"}]},
               "output": "2 edits applied"}])
    assert "c.py" in ev["files_read"]


def test_a_refused_read_grounds_nothing():
    ev = _ev([{"name": "read_file", "input": {"path": "secret.py"},
               "output": '{"error": "denied by the folder lock"}'}])
    assert ev["files_read"] == []


# ---------------------------------------------------------------------------
# A result that never arrived is not an observation
# ---------------------------------------------------------------------------

_CUT_OFF = {"name": "read_file", "input": {"path": "delfin/agent/engine.py"}}


def test_a_call_whose_result_never_arrived_grounds_nothing():
    assert _ev([_CUT_OFF])["files_read"] == []


def test_a_truncated_run_does_not_ground_more_than_a_clean_one():
    cut = set(_touched([_CUT_OFF]))
    refused = set(_touched([{**_CUT_OFF, "output": "(file not found)"}]))
    assert cut <= refused | set()
    assert cut == set()


def test_missing_bookkeeping_is_still_not_evidence_of_fabrication():
    """The lenient set keeps its old rule: nothing recorded, nothing
    concluded."""
    assert "delfin/agent/engine.py" in _ev([_CUT_OFF])["files"]


# ---------------------------------------------------------------------------
# What the child actually changed inside an isolated worktree
# ---------------------------------------------------------------------------

def test_the_worktree_diff_still_grounds_the_parent():
    """git --stat is the framework's own record of the child's edits, not
    the child's prose."""
    ev = sa.collect_report_evidence({
        "final_text": "x",
        "tool_calls": [],
        "worktree": {"diff_summary": " delfin/agent/engine.py | 4 ++--"},
    })
    assert any("engine.py" in f for f in ev["files_read"])
