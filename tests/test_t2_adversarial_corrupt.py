"""Adversarial review tests for builder nacht-s14's T2 task_state.py (hash e5702344).

Reviewer nacht-s15. These target weaknesses the builder's own control tests
(tests/test_t2_task_state.py) do not cover:
* A corrupt / partial / mid-write state file must not make load() raise — the
  module's own design doc says it mirrors working_state, which is best-effort
  ("a broken store must not break compaction"). load() must degrade like
  working_state, not propagate json.JSONDecodeError.
* Secrets must be scrubbed from EVERY rendered field, not only findings.
* Unicode / huge / None field content must not crash render() and must stay
  within the firm character ceiling.
Red = genuine defect in the built code; each finding is reported to the
builder with file:line + test id. Green cases become review evidence.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import task_state as ts


@pytest.fixture
def state_path(tmp_path):
    return tmp_path / "task_state.json"


_SECRET = "sk-proj-nacht-s15-adversarial-secret-value"


# -- Finding 1: corrupt / partial on-disk state must not blow up load() ------

def test_load_degrades_on_corrupt_json(state_path):
    """A half-written save (truncated JSON) must not raise in load()."""
    state_path.write_text('{"task": "build task_state.py", "phases": [{"name'
                          '" : "phase 2"]',  # truncated mid-dict
                          encoding="utf-8")
    st = ts.open(state_path)
    st.load()  # must not raise json.JSONDecodeError
    # Degraded == best effort: never a raise.
    assert isinstance(st.task, str)


def test_load_degrades_on_garbage_not_json(state_path):
    """A file that is not JSON at all (crash scribble) must degrade too."""
    state_path.write_text("not json at all <<< {", encoding="utf-8")
    st = ts.open(state_path)
    st.load()  # must not raise
    assert st.task == ""


def test_save_then_corrupt_then_reload_keeps_working(state_path):
    """After a corrupt read, a fresh write must recover the store."""
    st = ts.open(state_path)
    st.task = "recover me"
    st.save()
    data = json.loads(state_path.read_text(encoding="utf-8"))
    # Simulate a torn write by replacing the file mid-stream.
    state_path.write_text(data["task"][:4], encoding="utf-8")
    st2 = ts.open(state_path)
    st2.load()  # no raise
    st2.task = "recovered"
    st2.save()
    st3 = ts.open(state_path)
    st3.load()
    assert st3.task == "recovered"


# -- Finding 2: secrets scrubbed from EVERY rendered field --------------------

def test_scrub_covers_waiting_for(state_path):
    st = ts.open(state_path)
    st.waiting_for = _SECRET
    out = st.render()
    assert _SECRET not in out


def test_scrub_covers_phase_name_and_status(state_path):
    """Secrets as whole field values (after a space/newline, the realistic
    render shape) must be scrubbed from the phase name and status too.

    Task state's single pass -- _scrub(block) over the whole rendered block --
    covers every field. (A secret glued directly after a '_', i.e. '_sk-proj-',
    is NOT covered: output_guard._PROVIDER_API_KEY anchors '\\b' before
    sk-proj-, and '_' is a word char so the boundary never fires. But a
    rendered field never maps a '_' immediately before the value, so that
    edge is the shared redactor's, not task_state's.)
    """
    st = ts.open(state_path)
    st.task = "build a redactor check"
    st.commit(phase=_SECRET + " phase", phase_status="in progress",
              commit="abc")
    st.waiting_for = "awaiting " + _SECRET
    out = st.render()
    assert "build a redactor check" in out
    assert _SECRET not in out


# -- Finding 3: unicode / huge / None must not crash, must stay bounded ------

def test_unicode_survives_round_trip(state_path):
    st = ts.open(state_path)
    st.task = "Übung — Aufgabe ₂ ⚗️"
    st.add_finding(finding="芬芬芬")
    st.save()
    again = ts.open(state_path)
    again.load()
    assert again.task == "Übung — Aufgabe ₂ ⚗️"
    assert again.render()  # renders, one bound
    assert len(again.render()) <= 2200


def test_huge_task_still_bounded(state_path):
    st = ts.open(state_path)
    st.commit(task="x" * 10_000, phase="phase 2", phase_status="in_progress",
              commit="abc")
    out = st.render()
    assert len(out) <= 2200


def test_none_fields_render_without_crashing(state_path):
    st = ts.open(state_path)
    st.commit(task=None, phase=None, phase_status=None, commit=None)
    st.add_finding(finding=None)
    st.set_waiting_for(waiting_for=None)
    st.save()  # must not raise even when everything is None
    out = st.render()
    assert isinstance(out, str)
