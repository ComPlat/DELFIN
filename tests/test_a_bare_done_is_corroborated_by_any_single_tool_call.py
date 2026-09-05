""""Done" was backed by whatever tool call happened to be in the trace.

The mutation family was tightened once already: ``I implemented the fix
in engine.py`` used to be corroborated by ANY write in the trace, so
``mkdir -p /tmp/scratch`` certified a change to engine.py. It is now
paired with the file the claim names.

The generic family -- "done", "erledigt", "fertig", "abgeschlossen" --
kept the old rule: the only thing that could fail it was a trace with
ZERO tool calls. One ``task_update``, one ``bash_status`` poll, one
question to the user, and the delegate's claim that the work was finished
came back supported. It is now paired with the same thing every other
family is paired with: evidence in the trace that the run wrote, ran or
read something.
"""

from __future__ import annotations

from delfin.agent import subagents as sa


def _verdict(report: str, calls):
    return sa.verify_subagent_report({"result": report}, tool_calls=calls)


def _unsupported(verdict) -> list[str]:
    return [u["claim"] for u in verdict["unsupported"]]


def _supported(verdict) -> list[str]:
    return [s["claim"] for s in verdict["supported"]]


# ---------------------------------------------------------------------------
# The pairing
# ---------------------------------------------------------------------------

def test_a_bookkeeping_call_does_not_finish_the_work():
    v = _verdict("Done.", [
        {"name": "task_update", "input": {"task_id": 3, "status": "completed"},
         "output": '{"ok": true}'},
    ])
    assert _unsupported(v), "one bookkeeping call still certified 'done'"
    assert "wrote, ran or read" in v["unsupported"][0]["why"]


def test_polling_a_job_is_not_doing_the_work():
    v = _verdict("Erledigt.", [
        {"name": "bash_status", "input": {"job_id": "ab12"},
         "output": '{"running": true}'},
    ])
    assert _unsupported(v)


def test_asking_the_user_is_not_doing_the_work():
    v = _verdict("Fertig, alles abgeschlossen.", [
        {"name": "ask_user_question", "input": {"question": "welche Datei?"},
         "output": "the second one"},
    ])
    assert _unsupported(v)


# ---------------------------------------------------------------------------
# What still passes -- the rule must not swallow honest reports
# ---------------------------------------------------------------------------

def test_a_read_that_returned_something_backs_it():
    v = _verdict("Done — the parser is in editblock.py.", [
        {"name": "read_file", "input": {"path": "editblock.py"},
         "output": "def parse(...):"},
    ])
    assert any("Done" in c for c in _supported(v))


def test_a_command_that_ran_backs_it():
    v = _verdict("Done.", [
        {"name": "bash", "input": {"command": "python -c 'print(1)'"},
         "output": "1"},
    ])
    assert any("Done" in c for c in _supported(v))


def test_a_write_backs_it():
    v = _verdict("Erledigt.", [
        {"name": "write_file", "input": {"path": "notes.md"},
         "output": '{"status": "ok"}'},
    ])
    assert any("Erledigt" in c for c in _supported(v))


def test_the_no_calls_rule_still_speaks_first():
    """A trace with nothing in it keeps its own, clearer reason."""
    v = _verdict("Done.", [])
    assert v["unsupported"]
    assert "no tool calls at all" in v["unsupported"][0]["why"]


def test_the_other_families_are_unchanged():
    """The mutation and verification rules must keep their own wording."""
    v = _verdict("I fixed the bug in engine.py.", [
        {"name": "read_file", "input": {"path": "engine.py"},
         "output": "code"},
    ])
    assert any("file-writing tool call" in u["why"] for u in v["unsupported"])
