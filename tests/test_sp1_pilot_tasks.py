"""SP1 phase 2: pilot_tasks_sp1.yaml loads through the real loader and
its tasks are answerable under the citation bar.

Real records rule: the YAML is loaded through benchmark.load_tasks with
an explicit path (the way the CLI will load it), not through a local
parser. The citation-requirement check reads the prompts themselves;
the answerability check calls verify_citations against THIS repository
with a hand-built answer whose facts are the grep-verified locations
each task's comment names.
"""
from pathlib import Path

from delfin.agent.benchmark import load_tasks
from delfin.agent.explore_pilot import PilotTask, verify_citations, run_pilot

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / \
    "benchmark" / "pilot_tasks_sp1.yaml"


def _tasks():
    tasks = load_tasks(_PACK)
    assert tasks, "the pilot task file must load through load_tasks"
    return tasks


def test_six_tasks_load_with_required_fields():
    tasks = _tasks()
    assert len(tasks) == 6
    for t in tasks:
        assert t.id.startswith("sp1_")
        assert t.task_class == "pilot"
        assert t.mode == "solo"
        assert t.prompt
        assert t.max_tool_calls > 0


def test_every_prompt_requires_citations():
    """SA1 reviewer carry: without the requirement in the prompt the
    arms are compared at an unequal bar."""
    for t in _tasks():
        p = t.prompt.lower()
        assert "path:line" in p, f"{t.id}: prompt must demand path:line"
        assert "file.py:line" in p, f"{t.id}: must state the citation form"


def test_not_in_default_glob():
    """The file must stay OUTSIDE benchmark's tasks_*.yaml glob so the
    pilot suite never joins the production benchmark."""
    parent = _PACK.parent
    assert _PACK.name not in {p.name for p in parent.glob("tasks_*.yaml")}


# Facts per task id: (file, line, snippet that must be on that line).
# Each verified by grep at the time the task was written; the test
# re-checks them against the CURRENT tree so a drifted task fails here,
# not in a live measurement.
_FACTS = {
    "sp1_pilot_guard": [
        ("delfin/agent/explore_pilot.py", 37, "PILOT_TYPE"),
        ("delfin/agent/explore_pilot.py", 49, "check_delegate"),
    ],
    "sp1_control_parsing": [
        ("delfin/config.py", 569, "parse_control_text"),
        ("delfin/config.py", 632, "validate_control_text"),
    ],
    "sp1_benchmark_runner_entry": [
        ("delfin/agent/benchmark_runner.py", 1450, "run_task"),
        ("delfin/agent/benchmark.py", 323, "class Task"),
    ],
    "sp1_pilot_harness": [
        ("delfin/agent/explore_pilot.py", 255, "run_pilot"),
        ("delfin/agent/explore_pilot.py", 225, "verify_citations"),
    ],
    "sp1_scoring_result": [
        ("delfin/agent/benchmark.py", 1938, "score_outcome"),
        ("delfin/agent/benchmark.py", 520, "BenchmarkResult"),
    ],
    "sp1_pilot_wiring": [
        ("delfin/agent/explore_pilot.py", 167, "child_runner"),
        ("delfin/agent/explore_pilot.py", 352, "make_pilot_runner"),
    ],
}


def test_every_task_is_answerable_against_the_real_repo():
    """For each task: build the reference answer from the verified facts,
    run it through the full harness path (run_pilot -> verify_citations),
    and require citations_ok — the task set must be answerable."""
    tasks = {t.id: t for t in _tasks()}
    assert set(tasks) == set(_FACTS)
    for tid, facts in _FACTS.items():
        answer = " ".join(
            f"see {f}:{ln}" for f, ln, _ in facts)
        row = run_pilot(
            [PilotTask(task_id=tid, prompt=tasks[tid].prompt)],
            runner=lambda task, arm="plain", _a=answer: {
                "success": True, "final_text": _a})[0]
        assert row.arms["plain"].citations_ok, \
            f"{tid}: reference answer did not verify: " \
            f"{row.arms['plain'].citation_report.bad}"
        # And the fact still sits where the comment says: the line must
        # contain the named snippet.
        for f, ln, snippet in facts:
            line = (Path(".") / f).read_text(
                encoding="utf-8", errors="replace").splitlines()[ln - 1]
            assert snippet in line, \
                f"{tid}: fact drifted — {f}:{ln} no longer holds {snippet!r}"


def test_fakes_only_harness_run():
    """A full run_pilot over all 6 tasks with a fake runner: two arms,
    verdict computed, no model call."""
    tasks = [PilotTask(task_id=t.id, prompt=t.prompt) for t in _tasks()]

    def fake_runner(task, arm="plain"):
        # explore arm "runs" with fewer tool calls, same citations.
        tc = 1 if arm == "explore" else 4
        answer = "see delfin/agent/explore_pilot.py:37"
        return {"success": True, "final_text": answer,
                "tool_calls": tc, "duration_s": 1.0}

    rows = run_pilot(tasks, runner=fake_runner)
    assert len(rows) == 6
    for r in rows:
        assert set(r.arms) == {"plain", "explore"}
