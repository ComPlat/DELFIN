"""Adversarial pins on the multi-part task rule (phase 2, 07f5af55).

The builder's tests (tests/test_evidence_multipart.py) cover the happy
path and the declared scope. These pins check the CONTRACT boundaries an
independent reviewer would want pinned: what triggers the parts path,
what does not, and what survives re-entrancy. Every case here is correct
behaviour on the current code -- they would go RED if a future change
relaxed or over-fires the parts rule. The subject-ladder-skip gap is a
separate finding (sent to the builder), NOT pinned here.

Reference: delfin/agent/task_evidence.py --
  _PART_ITEM_RE (:275): a numbered ("1."/"1)") or lettered ("a)"/"a.")
  item at line start. Dash/bullet lists are deliberately NOT parts.
  _enumerated_parts (:279): >=2 items of the same kind form a parts list;
  a mixing of "1." with "a)" is a fragment and returns [].
  _part_verdict (:304): each part judged as a task of its own via the
  same ladder; the part's OWN description is empty, so parts never
  recurse (no double-merging of an inner list).
"""

from __future__ import annotations

from delfin.agent.task_evidence import check_completion_claim


def _check(subject, description="", **kw):
    return check_completion_claim(subject, description, **kw)


def test_a_single_numbered_line_is_not_a_parts_list():
    # One numbered line is "prose numbering", not a decomposition: the
    # subject ladder must run, not the parts path.
    r = _check("Fix the sort", "1. Do one thing", changes=[])
    assert r["verdict"] == "unmet"  # single-subject write task, no change
    # ... and crucially it is NOT an "unchecked" from an empty parts list.
    assert r["kind"] != "parts"


def test_a_dash_list_is_never_a_parts_list_even_with_work_words():
    # Dash lists are context notes ("use the gate"), not parts. So a
    # dash list must NOT turn a subject into a per-part check: the parts
    # path must not fire, and the subject must be judged by its own
    # ladder (here a write task whose change is not in the window, so
    # unmet -- NOT a vacuous unchecked from an empty parts list).
    desc = "- fix the sort\n- run the tests through the gate"
    r = _check("Fix the ordering bug", desc, changes=["delfin/report.py"])
    assert r["kind"] != "parts"
    assert r["verdict"] == "unmet"  # subject ladder ran; change not in window


def test_mixed_numbering_is_a_fragment_not_parts():
    # "1." mixing with "a)" is a fragment per the builder; the parts path
    # must NOT fire and the subject ladder must handle the subject alone.
    desc = "1. Fix the sort\n  a) also check the merge"
    r = _check("Fix the ordering bug", desc, changes=["delfin/report.py"])
    assert r["kind"] != "parts"


def test_parts_never_recurse_into_an_inner_list():
    # A part whose own text enumerates sub-items MUST NOT be decomposed a
    # second time. _PART_ITEM_RE only matches at line start, so the inner
    # "a)"/"b)" (mid-line) are not separate parts -- part "2" is one part
    # judged as a single subject. That part ("verify: run unit tests run
    # lint") has no evidence, so the aggregate is unmet and the note names
    # the unevidenced OUTER part (by its label "1"), not a deeper
    # sub-branch label -- the inner a)/b) were never split out.
    desc = ("1. fix the sort\n"
            "2. verify:  a) run unit tests  b) run lint")
    r = _check("Fix the ordering bug", desc, changes=[])
    assert r["kind"] == "parts"
    assert r["verdict"] == "unmet"  # part 1 unevidenced holds the task back
    note = str(r.get("note", ""))
    assert "part 1" in note  # names the OUTER part, not an inner a)/b) item


def test_an_unevidenced_part_holds_the_whole_task_back_and_is_named():
    # The wave-10 failure: 1 of 5 parts done called "fertig". When one
    # part is a concrete object (a written file) with no recorded write,
    # the aggregate must be unmet and the note must name that part.
    # Part 1 is evidenced by an in-window write; part 2 is not -- so the
    # note names part 2, the one without evidence.
    write1 = {"path": "delfin/report.py", "ts": 999, "created": True}
    desc = ("1. Write delfin/report.py\n"
            "2. Write delfin/report.md")
    r = _check("Deliver the report", desc,
               changes=[write1], observed=[], tests=None, window_start=0)
    assert r["verdict"] == "unmet"
    assert r["kind"] == "parts"
    note = str(r.get("note", ""))
    assert "part 2" in note and "report.md" in note
