"""Tests for the trace-derived behaviour classifier (behavioural-parity eval).

Three of the four target behaviours are ORDER-derived ("plan BEFORE act",
"read BEFORE edit", "run AFTER edit"), which the presence-only signal rubric
cannot express.  ``behavior_flags`` derives them from the ordered
``Trajectory.tool_calls`` using the agent's REAL snake_case tool vocabulary.

These tests run against hand-built synthetic trajectories — zero tokens, zero
live model — and exist to falsify classifier bugs (risk R2) BEFORE any live
run spends money.  They are the CI contract for the classifier's semantics.
"""

from __future__ import annotations

from delfin.agent import benchmark as bm
from delfin.agent.benchmark import Task, Trajectory


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _task(behavior: str, tid: str = "t") -> Task:
    return Task(id=tid, task_class="behavior", mode="solo", prompt="",
                behavior=behavior)


def _call(name: str, inp: str = "") -> dict:
    return {"name": name, "input": inp}


def _fp(path: str) -> str:
    """A tool_call input string carrying a file_path (as the runner emits)."""
    return '{"file_path": "%s"}' % path


def _traj(tool_calls: list[dict], text: str = "") -> Trajectory:
    return Trajectory(text=text, tool_calls=tool_calls)


# ---------------------------------------------------------------------------
# planned — a PLAN call before the first mutate AND first acting exec
# ---------------------------------------------------------------------------


def test_planned_true_when_task_create_precedes_edit():
    traj = _traj([_call("task_create"), _call("read_file", _fp("a.py")),
                  _call("edit_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("plan"), traj) == {"planned": 1}


def test_planned_true_via_exit_plan_mode():
    traj = _traj([_call("exit_plan_mode"), _call("write_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("plan"), traj) == {"planned": 1}


def test_planned_false_when_edit_precedes_plan():
    traj = _traj([_call("edit_file", _fp("a.py")), _call("task_create")])
    assert bm.behavior_flags(_task("plan"), traj) == {"planned": 0}


def test_planned_false_when_no_plan_call():
    traj = _traj([_call("read_file", _fp("a.py")),
                  _call("edit_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("plan"), traj) == {"planned": 0}


def test_planned_false_when_acting_exec_precedes_plan():
    # Running pytest (acting exec) before planning = did not plan first.
    traj = _traj([_call("bash", "pytest tests/ -q"), _call("task_create")])
    assert bm.behavior_flags(_task("plan"), traj) == {"planned": 0}


def test_planned_true_when_readonly_bash_precedes_plan():
    # A read-only scout (ls) before task_create must NOT disqualify planning.
    traj = _traj([_call("bash", "ls 01_xtb_opt/"), _call("task_create"),
                  _call("write_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("plan"), traj) == {"planned": 1}


# ---------------------------------------------------------------------------
# scouted — a read/scout before the first mutate
# ---------------------------------------------------------------------------


def test_scouted_true_read_before_edit():
    traj = _traj([_call("read_file", _fp("a.py")),
                  _call("edit_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("scout"), traj) == {"scouted": 1}


def test_scouted_false_edit_with_no_read():
    traj = _traj([_call("edit_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("scout"), traj) == {"scouted": 0}


def test_scouted_false_edit_before_read():
    traj = _traj([_call("edit_file", _fp("a.py")),
                  _call("read_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("scout"), traj) == {"scouted": 0}


def test_scouted_true_via_readonly_bash():
    traj = _traj([_call("bash", "cat 02_xtb_thermo/thermo.py"),
                  _call("edit_file", _fp("thermo.py"))])
    assert bm.behavior_flags(_task("scout"), traj) == {"scouted": 1}


def test_scouted_true_via_explore_subagent():
    traj = _traj([_call("subagent", '{"subagent_type": "explore"}'),
                  _call("edit_file", _fp("a.py"))])
    assert bm.behavior_flags(_task("scout"), traj) == {"scouted": 1}


def test_scouted_true_when_read_and_no_mutate():
    # Read-only investigation with no edit still counts as scouting.
    traj = _traj([_call("read_file", _fp("a.py")), _call("grep_file", "{}")])
    assert bm.behavior_flags(_task("scout"), traj) == {"scouted": 1}


# ---------------------------------------------------------------------------
# verified — after the LAST mutate, run something OR read the edited file back
# ---------------------------------------------------------------------------


def test_verified_true_run_after_write():
    traj = _traj([_call("write_file", _fp("s.py")),
                  _call("bash", "python s.py")])
    assert bm.behavior_flags(_task("verify"), traj) == {"verified": 1}


def test_verified_true_run_tests_after_edit():
    traj = _traj([_call("edit_file", _fp("m.py")), _call("run_tests", "{}")])
    assert bm.behavior_flags(_task("verify"), traj) == {"verified": 1}


def test_verified_false_edit_with_no_run():
    traj = _traj([_call("write_file", _fp("s.py"))])
    assert bm.behavior_flags(_task("verify"), traj) == {"verified": 0}


def test_verified_true_readback_of_edited_file():
    traj = _traj([_call("edit_file", _fp("thermo.py")),
                  _call("bash", "cat thermo.py")])
    assert bm.behavior_flags(_task("verify"), traj) == {"verified": 1}


def test_verified_false_read_back_wrong_file():
    # Reading a DIFFERENT file after editing does not verify the edit.
    traj = _traj([_call("edit_file", _fp("thermo.py")),
                  _call("bash", "cat other.py")])
    assert bm.behavior_flags(_task("verify"), traj) == {"verified": 0}


def test_verified_false_run_before_edit_only():
    # Running BEFORE the last mutate, then editing, is not verification.
    traj = _traj([_call("bash", "python s.py"),
                  _call("write_file", _fp("s.py"))])
    assert bm.behavior_flags(_task("verify"), traj) == {"verified": 0}


def test_verified_true_readfile_tool_readback():
    traj = _traj([_call("edit_file", _fp("thermo.py")),
                  _call("read_file", _fp("thermo.py"))])
    assert bm.behavior_flags(_task("verify"), traj) == {"verified": 1}


# ---------------------------------------------------------------------------
# asked — clarifying question (tool or prose) with NO guessing
# ---------------------------------------------------------------------------


def test_asked_true_via_tool():
    traj = _traj([_call("ask_user_question", "{}")])
    assert bm.behavior_flags(_task("ask"), traj) == {"asked": 1}


def test_asked_true_via_prose_question():
    traj = _traj([], text="Which molecule should I optimize? Please specify.")
    assert bm.behavior_flags(_task("ask"), traj) == {"asked": 1}


def test_asked_false_when_guessed_by_mutating():
    # Asked but then guessed (edited a file) → guessing overrides asking.
    traj = _traj([_call("ask_user_question", "{}"),
                  _call("edit_file", _fp("x.py"))])
    assert bm.behavior_flags(_task("ask"), traj) == {"asked": 0}


def test_asked_false_when_guessed_by_running():
    traj = _traj([_call("bash", "xtb mol.xyz --opt")],
                 text="I'll optimize the geometry.")
    assert bm.behavior_flags(_task("ask"), traj) == {"asked": 0}


def test_asked_false_no_question_no_action():
    traj = _traj([], text="Here is some general information about geometry.")
    assert bm.behavior_flags(_task("ask"), traj) == {"asked": 0}


# ---------------------------------------------------------------------------
# tagging / aggregation
# ---------------------------------------------------------------------------


def test_untagged_task_yields_no_flags():
    traj = _traj([_call("read_file", _fp("a.py"))])
    assert bm.behavior_flags(Task(id="u", task_class="x", mode="solo",
                                  prompt=""), traj) == {}


def test_only_tagged_behaviour_is_scored():
    # A scout task must not emit planned/verified/asked keys.
    traj = _traj([_call("read_file", _fp("a.py")),
                  _call("edit_file", _fp("a.py"))])
    flags = bm.behavior_flags(_task("scout"), traj)
    assert set(flags) == {"scouted"}


def test_score_outcome_populates_behavior_field():
    t = _task("scout")
    traj = _traj([_call("read_file", _fp("a.py")),
                  _call("edit_file", _fp("a.py"))])
    res = bm.score_outcome(t, traj, model="m")
    assert res.behavior == {"scouted": 1}
    # Behaviour flags must NOT leak into the quality components.
    assert "scouted" not in res.components
    assert res.quality_0_100 == sum(res.components.values())


def test_behavior_rates_aggregates_per_flag():
    t_scout = _task("scout", tid="s")
    good = bm.score_outcome(
        t_scout, _traj([_call("read_file", _fp("a.py")),
                        _call("edit_file", _fp("a.py"))]), model="m")
    bad = bm.score_outcome(
        t_scout, _traj([_call("edit_file", _fp("a.py"))]), model="m")
    t_ver = _task("verify", tid="v")
    ver = bm.score_outcome(
        t_ver, _traj([_call("write_file", _fp("s.py")),
                      _call("bash", "python s.py")]), model="m")
    rates = bm.behavior_rates([good, bad, ver])
    assert rates["scouted"] == {"rate": 0.5, "n": 2}
    assert rates["verified"] == {"rate": 1.0, "n": 1}


def test_behavior_rates_reads_dict_records():
    # behavior_rates must also work on JSONL-loaded dicts, not just objects.
    records = [
        {"behavior": {"planned": 1}},
        {"behavior": {"planned": 0}},
        {"behavior": {}},          # untagged contributes nothing
        {},                         # no behavior key at all
    ]
    rates = bm.behavior_rates(records)
    assert rates == {"planned": {"rate": 0.5, "n": 2}}


def test_aggregate_replicates_means_behavior_flags():
    t = _task("scout")
    r_hit = bm.score_outcome(
        t, _traj([_call("read_file", _fp("a.py")),
                  _call("edit_file", _fp("a.py"))]), model="m")
    r_miss = bm.score_outcome(
        t, _traj([_call("edit_file", _fp("a.py"))]), model="m")
    agg = bm.aggregate_replicates([r_hit, r_miss, r_hit])
    assert agg.behavior["scouted"] == 2 / 3


def test_mcp_namespaced_tool_names_are_normalized():
    # Live backends (KIT/MCP) emit ``mcp__<server>__<tool>`` names — the
    # classifier must strip the namespace or a real read/edit is invisible.
    # (Regression: KIT emitted mcp__delfin-docs__read_file → scouted=0.)
    scout = Trajectory(tool_calls=[
        {"name": "mcp__delfin-docs__read_file", "input": _fp("xtb_sample.out")},
    ])
    assert bm.behavior_flags(_task("scout"), scout) == {"scouted": 1}

    verify = Trajectory(tool_calls=[
        {"name": "mcp__delfin-fs__write_file", "input": _fp("s.py")},
        {"name": "mcp__delfin-sh__bash", "input": "python s.py"},
    ])
    assert bm.behavior_flags(_task("verify"), verify) == {"verified": 1}

    plan = Trajectory(tool_calls=[
        {"name": "mcp__delfin__task_create", "input": "{}"},
        {"name": "mcp__delfin-fs__edit_file", "input": _fp("a.py")},
    ])
    assert bm.behavior_flags(_task("plan"), plan) == {"planned": 1}


def test_dict_shaped_tool_input_is_handled():
    # The runner may hand inputs as dicts rather than JSON strings.
    traj = Trajectory(tool_calls=[
        {"name": "read_file", "input": {"file_path": "a.py"}},
        {"name": "edit_file", "input": {"file_path": "a.py"}},
    ])
    assert bm.behavior_flags(_task("scout"), traj) == {"scouted": 1}
