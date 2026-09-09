"""Called with no calc_id, get_calc_info returned somebody else's run.

Swept on 2026-09-09 by calling every tool in the catalogue with no
arguments at all -- the cheapest possible probe, and one nothing had ever
run. Most tools were fine. Six answered a malformed call as though it
were a finding about the world:

* ``get_calc_info`` returned a REAL calculation, the first in the index,
  with its functional, basis set and energies. A model that spells the
  key ``id`` instead of ``calc_id`` gets a confident, complete answer
  about a run it never asked for, and the grounding rules then have it
  cite those numbers in the answer. Not a refusal and not an error:
  plausible data about the wrong subject, which is the worst thing a
  scientific tool can hand back.
* ``search_docs`` answered ``results: []`` next to the corpus statistics
  -- which reads as "1539 sections searched, nothing on your topic". A
  negative finding about the documentation, from a search never made.
* ``list_files`` listed the entire tree: a required ``pattern`` silently
  defaulting to everything.
* ``list_sections`` and ``read_section`` answered "Document '' not found"
  and printed every document id in the index.
* ``history_get`` and ``history_search`` blamed the session for a missing
  ``ref`` / ``query``, so the diagnosis pointed at the harness.

Same shape as ``cron_delete`` answering "not_found" for a missing
``entry_id``, recorded a day earlier -- that one cost a wrong diagnosis
in a live session, and it is fixed here too, by the same guard rather
than one tool at a time.

The guard reads the schema, so a tool added tomorrow is covered without
anyone remembering to cover it. Which is the point: this was never a set
of six bugs, it was one missing check at a single entry point.
"""

from __future__ import annotations

import json

import pytest

import delfin.agent.api_client as A


SKIP = {
    # Would act, or reach the network, even with no arguments.
    "bash", "bash_background", "web_fetch", "web_search", "push_notification",
    "remote_trigger", "cron_create", "publish_report", "draft_email",
    "orchestrate", "subagent", "skill", "bash_kill", "watch_job",
    "schedule_wakeup", "worktree_merge", "enter_worktree", "exit_worktree",
    "delegate",
}


def _tools_with_required():
    for tool in A._DOC_TOOLS_OPENAI:
        fn = tool.get("function") or {}
        req = (fn.get("parameters") or {}).get("required") or []
        if req and fn.get("name") not in SKIP:
            yield fn["name"], list(req)


def _perms(tmp_path):
    return A.KitToolPermissions(mode="default", workspace=str(tmp_path))


@pytest.mark.parametrize("name,required", list(_tools_with_required()),
                         ids=lambda v: v if isinstance(v, str) else "")
def test_every_tool_refuses_a_call_with_nothing_in_it(name, required, tmp_path):
    """And says WHICH argument, so the model fixes the call rather than
    concluding something about the data."""
    out = A._doc_executor.execute(name, {}, _perms(tmp_path))
    try:
        err = json.loads(out).get("error", "")
    except (json.JSONDecodeError, AttributeError):
        err = str(out)
    assert err, f"{name} answered a call with no arguments: {str(out)[:200]}"
    assert any(r in err for r in required), f"{name}: {err[:200]}"
    assert "required" in err.lower(), f"{name}: {err[:200]}"


def test_the_two_cases_that_were_reported_before_the_sweep(tmp_path):
    """Both were found in the field, written down as separate fixes, and
    never applied. Neither needed its own fix in the end.

    cron_delete answered "not_found" for a missing entry_id, so passing
    `id` instead of `entry_id` -- a real mistake, made in a real session
    -- read as "that entry is gone" and cost a wrong diagnosis.

    schedule_wakeup accepted an EMPTY prompt and scheduled a real wake-up
    with nothing to do: an agent woken at 3am to carry out an empty
    instruction.
    """
    perms = _perms(tmp_path)
    wrong_key = A._doc_executor.execute("cron_delete", {"id": 5}, perms)
    assert "entry_id is required" in str(wrong_key)
    assert "not_found" not in str(wrong_key)

    for prompt in ("", "   "):
        out = A._doc_executor.execute(
            "schedule_wakeup",
            {"delay_seconds": 600, "prompt": prompt, "reason": "r"}, perms)
        assert "prompt is required" in str(out), prompt


def test_the_message_says_nothing_was_looked_up(tmp_path):
    """The failure this replaces was a model believing a negative result.
    The message has to close that reading explicitly."""
    out = A._doc_executor.execute("get_calc_info", {}, _perms(tmp_path))
    err = json.loads(out)["error"]
    assert "calc_id" in err
    assert "nothing was looked up" in err


# ---------------------------------------------------------------------------
# ...without refusing calls that are fine
# ---------------------------------------------------------------------------

def test_an_empty_value_that_means_something_still_runs(tmp_path):
    """`new_string: ""` deletes the block it matched and `content: ""`
    creates an empty file. Both work on main; a guard that treats every
    empty required argument as missing would have broken them, which is a
    regression written in the name of a fix."""
    f = tmp_path / "f.txt"
    f.write_text("keep\nDROP\nkeep2\n", encoding="utf-8")
    perms = A.KitToolPermissions(mode="bypassPermissions",
                                 workspace=str(tmp_path))
    A._doc_executor.execute("read_file", {"path": "f.txt"}, perms)
    out = A._doc_executor.execute(
        "edit_file",
        {"path": "f.txt", "old_string": "DROP\n", "new_string": ""}, perms)
    assert "error" not in str(out).lower()[:40]
    assert f.read_text(encoding="utf-8") == "keep\nkeep2\n"

    out = A._doc_executor.execute(
        "write_file", {"path": "e.txt", "content": ""}, perms)
    assert (tmp_path / "e.txt").read_text(encoding="utf-8") == ""


def test_zero_is_a_value(tmp_path):
    """`cell_idx: 0` is the first cell, not a missing argument."""
    assert A._argument_present(0) is True
    assert A._argument_present(False) is True
    assert A._argument_present("") is False
    assert A._argument_present(None) is False
    assert A._argument_present([]) is False


def test_the_alias_a_weak_model_writes_still_counts(tmp_path):
    """`file_path` for `path` has been tolerated for as long as there have
    been weak models. A guard reading the schema literally would have
    taken that back."""
    (tmp_path / "f.txt").write_text("hello\n", encoding="utf-8")
    perms = _perms(tmp_path)
    out = A._doc_executor.execute("read_file", {"file_path": "f.txt"}, perms)
    assert "hello" in str(out)
    assert A._missing_required_argument("read_file", {"file_path": "x"}) is None
    assert A._missing_required_argument("read_file", {}) is not None


def test_a_tool_with_no_required_arguments_is_untouched():
    assert A._missing_required_argument("calc_summary", {}) is None
    assert A._missing_required_argument("list_docs", {}) is None


def test_an_unknown_tool_is_left_to_the_dispatcher(tmp_path):
    """The near-miss suggestion for a hallucinated tool name must still be
    what a model gets, rather than a complaint about arguments."""
    out = A._doc_executor.execute("read_fle", {}, _perms(tmp_path))
    assert "Unknown tool" in str(out)


def test_the_guard_never_raises():
    for args in (None, [], "not a dict", {"a": object()}):
        A._missing_required_argument("get_calc_info", args)  # type: ignore
    A._missing_required_argument("", {})
