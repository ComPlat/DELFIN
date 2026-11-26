"""A scored run's acceptance verdict is about the task that was scored.

The ten r3 long tasks share ONE acceptance script and ONE ``verify:`` name,
so before this was fixed every run of any of them was graded against all
ten fixtures at once.  A model that finished its own task was rejected
because nine fixtures it had never been shown were still unsolved: the
measured completion rate for the suite could not rise above zero, whatever
the model did.

Both directions are asserted here, through the REAL entry point
(``benchmark_runner.run_acceptance``, which launches the script the way a
scored run does):

* with the task id exported, the workspace that solved exactly that task
  is ACCEPTED -- this is the half that used to fail;
* with the task id exported, a workspace that has NOT solved it is still
  REJECTED, and the reason names the task -- so the fix does not become a
  gate that accepts anything;
* with no id exported the script keeps its documented whole-suite meaning
  (``check(task_id=None)``), which the gate test relies on.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "delfin" / "agent" / "pack" / "benchmark"

MEDIAN = "r3long-code-median"


def _setup_mod():
    sys.path.insert(0, str(BENCH / "setup"))
    try:
        import r3long_setup as setup_mod
    finally:
        sys.path.remove(str(BENCH / "setup"))
    return setup_mod


def _accept_mod():
    path = BENCH / "accept" / "r3long_accept.py"
    spec = importlib.util.spec_from_file_location("r3long_accept_scope", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def ws(tmp_path):
    setup = _setup_mod()
    assert setup.main([setup.__name__, str(tmp_path)]) == 0
    return tmp_path


def _solve_median(ws: Path) -> None:
    """Do exactly the work r3long-code-median asks for, and nothing else."""
    d = ws / "r3long" / "code" / "median"
    core = (d / "stats_core.py").read_text()
    fixed = core.replace(
        "    return s[n // 2]  # DEFECT: must be "
        "(s[n // 2 - 1] + s[n // 2]) / 2",
        "    return (s[n // 2 - 1] + s[n // 2]) / 2",
    )
    assert fixed != core, "the fixture no longer carries the defect to fix"
    (d / "stats_core.py").write_text(fixed)
    # The note the task asks for.  Its required shape is checked by the
    # acceptance script (module name plus both expected values); writing it
    # by hand rather than through a helper keeps this test independent of
    # that script's internals.
    (d / "fix_note.md").write_text(
        "changed module: stats_core.py\n"
        "corrected values: 2.5, 7.0\n",
        encoding="utf-8",
    )


def test_an_exported_task_id_accepts_the_workspace_that_solved_it(ws, monkeypatch):
    """The half that was broken: one task done, nine untouched, ACCEPTED."""
    from delfin.agent.benchmark_runner import run_acceptance

    _solve_median(ws)

    # Without an id this is the whole-suite question, and it must still say
    # no -- nine fixtures are unsolved.  This is the control: it proves the
    # workspace really is only partly done, so the next assertion is about
    # scope and not about the fixture being trivially green.
    monkeypatch.delenv("DELFIN_BENCH_TASK_ID", raising=False)
    whole_ok, whole_out = run_acceptance("r3long_accept.py", ws, root=ROOT)
    assert whole_ok is False, whole_out
    assert MEDIAN not in whole_out.split(":")[0], whole_out

    monkeypatch.setenv("DELFIN_BENCH_TASK_ID", MEDIAN)
    ok, out = run_acceptance("r3long_accept.py", ws, root=ROOT)
    assert ok is True, out
    assert MEDIAN in out, out


def test_an_exported_task_id_still_rejects_the_unsolved_task(ws, monkeypatch):
    """The never-pass direction: scoping must not become accepting."""
    from delfin.agent.benchmark_runner import run_acceptance

    monkeypatch.setenv("DELFIN_BENCH_TASK_ID", MEDIAN)
    ok, out = run_acceptance("r3long_accept.py", ws, root=ROOT)
    assert ok is False, out
    assert MEDIAN in out, out


def test_an_unknown_exported_task_id_is_rejected_loudly(ws, monkeypatch):
    """A typo must not silently fall back to grading all ten."""
    from delfin.agent.benchmark_runner import run_acceptance

    monkeypatch.setenv("DELFIN_BENCH_TASK_ID", "r3long-no-such-task")
    ok, out = run_acceptance("r3long_accept.py", ws, root=ROOT)
    assert ok is False, out
    assert "r3long-no-such-task" in out, out


def test_no_exported_id_keeps_the_documented_whole_suite_meaning(ws):
    """``check(task_id=None)`` grades all ten -- the gate test relies on it."""
    acc = _accept_mod()
    ok, out = acc.check(ws)
    assert ok is False
    assert out.startswith(MEDIAN + ":"), out


def test_the_runner_exports_the_id_it_is_scoring(monkeypatch):
    """The wiring: a scored run must hand the id to the acceptance script.

    Asserted on ``_child_env`` because that is the dict the subprocess is
    given; a value set only in the parent would never reach the script.
    """
    from delfin.agent import benchmark_runner as br

    env = br._child_env(ROOT, task_id=MEDIAN)
    assert env.get(br.BENCH_TASK_ID_ENV) == MEDIAN

    # No id asked for and none in the parent: the child gets none either.
    monkeypatch.delenv(br.BENCH_TASK_ID_ENV, raising=False)
    assert br.BENCH_TASK_ID_ENV not in br._child_env(ROOT)

    # A caller that sets the variable itself is not overruled.
    monkeypatch.setenv(br.BENCH_TASK_ID_ENV, MEDIAN)
    assert br._child_env(ROOT)[br.BENCH_TASK_ID_ENV] == MEDIAN
