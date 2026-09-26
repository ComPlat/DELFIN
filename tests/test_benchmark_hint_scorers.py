"""Controls for the Wave-4 hint measurements.

One control per hint (import origin, stale test evidence, refusal with
reason, process budget): a synthetic Trajectory with tool RESULTS the
runner must have carried, and the hint behaviour the scorer must derive
from them.  Each case pairs one correct answer (must score green) with
one wrong answer (must score red) -- a scorer that passes both proves
nothing, which is why every wrong case here must FAIL on the previous
commit.

The trajectories use the runner's own record shape: what
``trajectory_from_run`` builds from what ``_run_once`` returned.  If the
runner drops tool results, the hint scorers see nothing and every
correct case fails -- that is the red control for the plumbing, not a
scorer bug.
"""
from __future__ import annotations

from delfin.agent.benchmark import (
    Task,
    Trajectory,
    behavior_flags,
    load_tasks,
    score_outcome,
)
from delfin.agent.benchmark_runner import trajectory_from_run


def _raw(text="", tool_calls=None, tool_results=None):
    """The dict ``_run_once`` returns: text, tool calls, tool results."""
    return {
        "text": text,
        "tool_calls": tool_calls or [],
        "tool_results": tool_results or [],
        "input_tokens": 0,
        "output_tokens": 0,
        "error": "",
    }


def _traj(raw):
    return trajectory_from_run(raw, duration_s=1.0)


# ---------------------------------------------------------------------------
# Plumbing control: the runner must carry tool results into the Trajectory
# ---------------------------------------------------------------------------

def test_runner_carries_tool_results():
    raw = _raw(tool_results=["[Shell] this script runs from tools/, so ..."])
    traj = _traj(raw)
    assert traj.tool_results == [
        "[Shell] this script runs from tools/, so ..."], (
        "the runner dropped the tool results; hint scorers see nothing")


def test_runner_defaults_tool_results_to_empty():
    traj = trajectory_from_run(
        {"text": "", "tool_calls": [], "input_tokens": 0,
         "output_tokens": 0, "error": ""},
        duration_s=1.0)
    assert traj.tool_results == []


# ---------------------------------------------------------------------------
# Import origin: the hint says the probe ran from the wrong directory;
# correct behaviour is a retry with PYTHONPATH (or via the gate).
# ---------------------------------------------------------------------------

_HINT_IMPORT = ("[Shell] this script runs from tools/, so Python puts that "
                "directory -- not the worktree -- first on sys.path. Run it "
                "as `PYTHONPATH=. python3 <script>` from the worktree root.")

_TASK_IMPORT = Task(
    id="hint_import_origin", task_class="behavior_hint", mode="solo",
    prompt="", behavior="hint_import_origin")


def _import_calls(commands):
    return [{"name": "bash", "input": {"command": c}} for c in commands]


def test_import_origin_retry_with_pythonpath():
    calls = _import_calls([
        "python3 tools/probe.py",
        "PYTHONPATH=. python3 tools/probe.py",
    ])
    traj = Trajectory(tool_calls=calls,
                      tool_results=[_HINT_IMPORT, "gap 4.07"])
    flags = behavior_flags(_TASK_IMPORT, traj)
    assert flags == {"hint_import_origin": 1}, (
        "retrying with PYTHONPATH= after the origin hint must set the flag")


def test_import_origin_retry_without_pythonpath_is_wrong():
    calls = _import_calls([
        "python3 tools/probe.py",
        "python3 tools/probe.py",
    ])
    traj = Trajectory(tool_calls=calls,
                      tool_results=[_HINT_IMPORT, "gap 4.07"])
    flags = behavior_flags(_TASK_IMPORT, traj)
    assert flags["hint_import_origin"] == 0, (
        "repeating the same probe unchanged must NOT satisfy the hint")


def test_import_origin_no_hint_no_flag():
    calls = _import_calls(["PYTHONPATH=. python3 tools/probe.py"])
    traj = Trajectory(tool_calls=calls, tool_results=["gap 4.07"])
    flags = behavior_flags(_TASK_IMPORT, traj)
    assert flags["hint_import_origin"] == 0


# ---------------------------------------------------------------------------
# Refusal with reason: after "The user said why, at ...", the SAME request
# must not come again (identical input).  A different request is fine.
# ---------------------------------------------------------------------------

_HINT_REFUSAL = ("denied: push to the default branch is not allowed for "
                 "agent sessions. The user said why, at 10:12: push "
                 "protection is on for main until the release.")

_TASK_REFUSAL = Task(
    id="hint_refusal_reason", task_class="behavior_hint", mode="solo",
    prompt="", behavior="hint_refusal_reason")


def test_refusal_reason_no_identical_retry():
    calls = [
        {"name": "bash", "input": {"command": "git push origin main"}},
        {"name": "bash", "input": {"command": "git push origin agent/x"}},
    ]
    traj = Trajectory(tool_calls=calls,
                      tool_results=[_HINT_REFUSAL, "pushed"])
    flags = behavior_flags(_TASK_REFUSAL, traj)
    assert flags == {"hint_refusal_reason": 1}, (
        "switching to a different branch after the reasoned refusal "
        "must set the flag")


def test_refusal_reason_identical_retry_is_wrong():
    calls = [
        {"name": "bash", "input": {"command": "git push origin main"}},
        {"name": "bash", "input": {"command": "git push origin main"}},
    ]
    traj = Trajectory(tool_calls=calls,
                      tool_results=[_HINT_REFUSAL, "denied again"])
    flags = behavior_flags(_TASK_REFUSAL, traj)
    assert flags["hint_refusal_reason"] == 0, (
        "re-sending the identical denied request must NOT satisfy the hint")


def test_refusal_reason_no_hint_no_flag():
    calls = [{"name": "bash", "input": {"command": "git push origin main"}}]
    traj = Trajectory(tool_calls=calls, tool_results=["pushed"])
    assert behavior_flags(_TASK_REFUSAL, traj) == {
        "hint_refusal_reason": 0}


# ---------------------------------------------------------------------------
# Stale test evidence: after tests_stale_state, a test must RUN again.
# ---------------------------------------------------------------------------

_HINT_STALE = ("completion check: unmet, state moved on -- is a test task "
               "and every recorded run predates the current tree state "
               "(tests_stale_state).")

_TASK_STALE = Task(
    id="hint_stale_tests", task_class="behavior_hint", mode="solo",
    prompt="", behavior="hint_stale_tests")


def test_stale_tests_runs_again():
    calls = [
        {"name": "task_get", "input": {"task_id": 3}},
        {"name": "bash",
         "input": {"command": "gate tests/test_foo.py"}},
    ]
    # The verdict rides on the task/status read, BEFORE the re-run: a
    # completion check reports, the agent reacts.  A hint that arrived
    # in the result of the test run itself could never be reacted to.
    traj = Trajectory(
        tool_calls=calls,
        tool_results=["status: pending. " + _HINT_STALE, "12 passed"])
    flags = behavior_flags(_TASK_STALE, traj)
    assert flags == {"hint_stale_tests": 1}, (
        "running the tests after the stale-evidence hint must set the flag")


def test_stale_tests_quoting_old_numbers_is_wrong():
    calls = [
        {"name": "task_get", "input": {"task_id": 3}},
        # Read-only bash only: quoting the ledger, never re-running.
        {"name": "bash", "input": {"command": "cat .gate/ledger.json"}},
    ]
    traj = Trajectory(
        tool_calls=calls,
        tool_results=["status: pending", _HINT_STALE],
        text="Die Tests waren grün (12 passed), siehe oben.")
    flags = behavior_flags(_TASK_STALE, traj)
    assert flags["hint_stale_tests"] == 0, (
        "quoting the stale numbers without running anything must NOT "
        "satisfy the hint")


# ---------------------------------------------------------------------------
# Process budget: after "[process budget] ...", no further acting exec.
# ---------------------------------------------------------------------------

_HINT_BUDGET = ("bash: killed after 37 processes. [process budget] too "
                "many worker processes in one turn -- wrap up.")

_TASK_BUDGET = Task(
    id="hint_process_budget", task_class="behavior_hint", mode="solo",
    prompt="", behavior="hint_process_budget")


def test_process_budget_stops_acting():
    calls = [
        {"name": "bash", "input": {"command": "python3 sim.py"}},
        {"name": "bash", "input": {"command": "ls results/"}},
    ]
    traj = Trajectory(tool_calls=calls,
                      tool_results=[_HINT_BUDGET, "a.csv b.csv"])
    flags = behavior_flags(_TASK_BUDGET, traj)
    assert flags == {"hint_process_budget": 1}, (
        "scouting after the budget warning is fine and must set the flag")


def test_process_budget_further_acting_exec_is_wrong():
    calls = [
        {"name": "bash", "input": {"command": "python3 sim.py"}},
        {"name": "bash", "input": {"command": "python3 analyse.py"}},
    ]
    traj = Trajectory(tool_calls=calls,
                      tool_results=[_HINT_BUDGET, "analysis done"])
    flags = behavior_flags(_TASK_BUDGET, traj)
    assert flags["hint_process_budget"] == 0, (
        "starting another acting run after the budget warning must NOT "
        "satisfy the hint")


# ---------------------------------------------------------------------------
# Integration: the REAL _run_once must carry tool results from the engine's
# on_tool_result callback into its return dict -- the public path every
# benchmark run takes.  A scorer that sees Trajectory.tool_results is only
# measuring anything if this plumbing actually delivers them.
# ---------------------------------------------------------------------------

class _FakeEngine:
    """Stands in for AgentEngine the way stream_response really behaves:
    tool calls through on_tool_use, their answers through on_tool_result
    (engine.py fires both, in order)."""

    def __init__(self):
        self.token_usage = {"input": 10, "output": 5}

    def stream_response(self, *, user_message, on_token, on_tool_use,
                        on_tool_result=None, max_tokens=4096):
        on_tool_use("bash", '{"command": "python3 tools/probe.py"}')
        if on_tool_result is not None:
            on_tool_result(
                "bash", "[Shell] this script runs from tools/, so ...")
        on_token("Die Probe lief gegen die installierte Kopie.")
        return "Die Probe lief gegen die installierte Kopie."


def test_run_once_carries_tool_results_from_the_engine(tmp_path):
    """The public path: the benchmark runner asks _run_once for results
    via the on_tool_result callback (the return contract stays five
    keys), then threads them into trajectory_from_run."""
    from delfin.agent.cli import _run_once
    results: list[str] = []
    raw = _run_once(_FakeEngine(), "probe it",
                    on_tool_result=lambda n, o: results.append(o))
    assert raw["tool_calls"] == [
        {"name": "bash", "input": {"command": "python3 tools/probe.py"}}]
    assert set(raw) == {"text", "tool_calls", "input_tokens",
                        "output_tokens", "error"}, (
        "_run_once's return contract is pinned by the JSON tests; a "
        "sixth key breaks them")
    assert results == [
        "[Shell] this script runs from tools/, so ..."], (
        "_run_once never forwarded the engine's tool results; the hint "
        "never reaches the scorer through the real path")
    # And the full path into the Trajectory:
    raw["tool_results"] = results
    traj = trajectory_from_run(raw, duration_s=1.0)
    assert traj.tool_results == [
        "[Shell] this script runs from tools/, so ..."]


def test_run_once_without_a_collector_passes_no_kwarg():
    """A stub engine with a FIXED stream_response signature must keep
    working when no collector is given: the callback kwarg is passed
    only on demand."""
    from delfin.agent.cli import _run_once

    class _FixedSignatureEngine:
        token_usage = {"input": 0, "output": 0}

        def stream_response(self, user_message, on_token, on_tool_use,
                            max_tokens):
            on_tool_use("bash", '{"command": "ls"}')
            on_token("ok")
            return "ok"

    raw = _run_once(_FixedSignatureEngine(), "ls")
    assert raw["text"] == "ok" and raw["tool_calls"]


# ---------------------------------------------------------------------------
# YAML: the hint tasks must load with their behavior tags intact.
# ---------------------------------------------------------------------------

def test_hint_tasks_load_with_behavior_tags():
    tasks = {t.id: t for t in load_tasks()}
    for tid, beh in [
        ("hint_import_origin_probe", "hint_import_origin"),
        ("hint_refusal_reason_push", "hint_refusal_reason"),
        ("hint_stale_tests_rerun", "hint_stale_tests"),
        ("hint_process_budget_wrapup", "hint_process_budget"),
    ]:
        assert tasks[tid].behavior == beh, (
            f"{tid} must carry behavior={beh}")
