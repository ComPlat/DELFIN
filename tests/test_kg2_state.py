"""KG2 phase 2: one state, two writers.

TaskState gains the fields the wave brief names -- success_criterion,
done-with-evidence, decisions, the person's latest instruction
(verbatim, marked as the person's) -- plus event hooks fed from tool
events (task tool updates, commits, test runs, registered waits) with
NO model call. The same store is rendered by working_state at
compaction, so the person's instruction and the open work survive
compaction and restart identically.

Red control: every case here ran RED on 019c6a31 (fields missing,
render silent about them). After the change each case goes green.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import task_state as ts


def _state(tmp_path):
    return ts.open(tmp_path / "task_state.json")


def test_new_fields_round_trip(tmp_path):
    st = _state(tmp_path)
    st.success_criterion = "gate green on all three phases"
    st.set_person_instruction("run the fast suite, skip the slow markers")
    st.decisions = [
        {"decision": "no third state mechanism", "reason": "brief forbids it"}
    ]
    st.done.append(
        {"what": "phase 1 inventory", "evidence": "sent to s15 + s29"}
    )
    st.save()
    other = ts.open(tmp_path / "task_state.json")
    other.load()
    assert other.success_criterion == "gate green on all three phases"
    assert other.person_instruction == (
        "run the fast suite, skip the slow markers")
    assert other.person_instruction_is_fresh is True
    assert other.decisions == [
        {"decision": "no third state mechanism",
         "reason": "brief forbids it"}]
    assert other.done == [
        {"what": "phase 1 inventory", "evidence": "sent to s15 + s29"}]


def test_person_instruction_replaces_task_and_is_marked(tmp_path):
    st = _state(tmp_path)
    st.commit(task="old goal from a plan")
    st.set_person_instruction("stop that, work on the importer instead")
    rendered = st.render()
    # The person's newest instruction is the goal and is marked as theirs.
    assert "stop that, work on the importer instead" in rendered
    assert "Person's latest instruction" in rendered
    assert st.person_instruction_is_fresh is True
    # A fresh person instruction replaces the stored plan as the goal.
    assert st.task == "stop that, work on the importer instead"


def test_person_instruction_marked_stale_on_commit(tmp_path):
    st = _state(tmp_path)
    st.set_person_instruction("build phase 2")
    assert st.person_instruction_is_fresh is True
    st.commit(commit="abc12345")
    # Once the agent has acted (a commit is recorded) the instruction is
    # no longer fresh: a stored plan never resumes against a newer
    # instruction, and an executed instruction is no longer presented
    # as the live steering one.
    assert st.person_instruction_is_fresh is False


def test_render_bounded_with_hostile_fields(tmp_path):
    st = _state(tmp_path)
    st.commit(task="T" * 5000)
    st.success_criterion = "S" * 5000
    st.set_person_instruction("P" * 5000)
    for _ in range(50):
        st.decisions.append({"decision": "D" * 4000, "reason": "R" * 4000})
        st.done.append({"what": "W" * 4000, "evidence": "E" * 4000})
    rendered = st.render()
    assert len(rendered) <= ts._MAX_RENDER_CHARS + 200
    assert rendered  # something came out


def test_render_scrubs_secrets(tmp_path):
    st = _state(tmp_path)
    st.set_person_instruction("use the api_key: AKIAEXAMPLE123456 and go")
    rendered = st.render()
    assert "AKIAEXAMPLE123456" not in rendered


def test_missing_fields_degrade_to_defaults(tmp_path):
    # An OLD state file (written before this change) loads without error.
    path = tmp_path / "old.json"
    path.write_text(json.dumps({
        "task": "old session", "waiting_for": "", "phases": [],
        "commits": [], "findings": [],
    }), encoding="utf-8")
    st = ts.open(path)
    st.load()
    assert st.success_criterion == ""
    assert st.person_instruction == ""
    assert st.person_instruction_is_fresh is False
    assert st.decisions == []
    assert st.done == []


def test_working_state_renders_person_instruction_from_store(tmp_path,
                                                            monkeypatch):
    # Second writer, same store: working_state picks up the person's
    # instruction from the task_state store at compaction time.
    st = ts.open(tmp_path / "task_state.json")
    st.set_person_instruction("verify against the ORCA manual first")
    st.save()

    from delfin.agent import working_state as ws
    block = ws.build_working_state_block([], workspace=tmp_path)
    assert "verify against the ORCA manual first" in block
    assert "person's latest instruction" in block
