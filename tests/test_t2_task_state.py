"""T2 phase tests: persistent per-session state file + bounded render.

A long session resumes with >900k input tokens re-sent every turn. T2 adds
``delfin/agent/task_state.py``: a small JSON state file per session (task,
phases with status, commits, open findings, "waiting for") written after each
commit/handoff, whose ``render()`` hands back ONE bounded, deterministic,
secret-scrubbed block. These tests pin that contract on the UNCHANGED code:
the module does not exist yet, so every test is red before the fix.
"""
from __future__ import annotations

import json

import pytest


# A token that must never leak into a rendered block, whatever its owner.
# Shaped as an OpenAI project key so the house redactor catches it. Built at
# runtime so no contiguous credential token sits in this source file.
_SECRET = "sk-proj-" + "nacht-s14-super-secret-value"


@pytest.fixture
def state_path(tmp_path):
    return tmp_path / "task_state.json"


def _fill(state):
    """A representative state: past a first commit, mid-phase-2, waiting."""
    state.begin()
    state.commit(task="build task_state.py", phase="phase 2",
                 phase_status="in_progress", commit="abc123")
    state.add_finding(finding="render() leaks " + _SECRET)
    state.set_waiting_for(waiting_for="reviewer approval of phase 2")
    return state


def test_state_file_written_to_disk(state_path):
    from delfin.agent import task_state
    st = task_state.open(state_path)
    _fill(st)
    st.save()
    assert state_path.exists()
    data = json.loads(state_path.read_text(encoding="utf-8"))
    assert data["task"] == "build task_state.py"
    assert data["waiting_for"] == "reviewer approval of phase 2"


def test_state_round_trips_read_back(state_path):
    from delfin.agent import task_state
    st = task_state.open(state_path)
    _fill(st)
    st.save()
    again = task_state.open(state_path)
    again.load()
    # Same task, phase, commit and finding survive a read-back.
    assert again.task == "build task_state.py"
    assert any(p["name"] == "phase 2" and p["status"] == "in_progress"
               for p in again.phases)
    assert "abc123" in again.commits
    assert again.waiting_for == "reviewer approval of phase 2"


def test_render_names_task_and_open_phase(state_path):
    from delfin.agent import task_state
    st = task_state.open(state_path)
    _fill(st)
    assert "build task_state.py" in st.render()
    assert "phase 2" in st.render()
    assert "in_progress" in st.render()
    assert "abc123" in st.render()


def test_render_is_bounded(state_path):
    from delfin.agent import task_state
    st = task_state.open(state_path)
    for i in range(200):  # far more commits than should ever be shown
        st.commit(task="bulk", phase="phase 2", phase_status="in_progress",
                  commit=f"bulk{i:06x}")
    rendered = st.render()
    assert len(rendered) <= 2200, "render() must stay within the fixed budget"
    assert "bulk000000" not in rendered  # oldest dropped by the newest-first cap
    assert "bulk0000c7" in rendered  # the newest commit stays visible (i=199 == 0xc7)


def test_render_is_deterministic(state_path):
    from delfin.agent import task_state
    st = task_state.open(state_path)
    _fill(st)
    assert st.render() == st.render()


def test_render_never_leaks_secret(state_path):
    from delfin.agent import task_state
    st = task_state.open(state_path)
    _fill(st)
    assert _SECRET not in st.render()
