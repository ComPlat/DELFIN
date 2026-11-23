"""Two w36-r2 mounts, each a call on a path the product runs.

1. ``grade_run`` accepts ``max_cores`` (w21_run_verify.py:126 passes the
   task's 2-core budget) but its body never consulted it: a run whose
   CONTROL requests PAL 8 passed the whole-run acceptance. The helper
   that does the job, ``resources_within``, was test-only. A green run
   with the budget kept is the control on both sides.

2. The five D1 data tasks name ``fixtures/trap_*.csv`` in their prompts
   and carry ``fixture:`` keys, but no task declared a ``setup:``, so
   the CSVs were built by nothing: the question was never put. The
   fixture builder existed (setup/w17_data_fixtures.py) and was
   test-only. The mount is the setup key the original commit described
   (1307d8ac3: "fixtures built by setup/w17_data_fixtures.py").
"""

from __future__ import annotations

from pathlib import Path

import yaml

from delfin.agent.benchmark import load_tasks
from delfin.agent.pack.benchmark.accept import w20_common
from delfin.agent.run_evidence import Evidence

_ROOT = Path(__file__).resolve().parents[1]


def _ev(**control) -> Evidence:
    """One terminated, input-carrying step + CONTROL. The step clears
    grade_run's other checks so the resource verdict is the only thing
    that can move the outcome.

    The step names its OUTPUT as well as its input. A step that claims
    to have terminated normally while naming no output is a shape
    ``run_evidence.collect`` never produces -- it reads the termination
    state OUT of the output file -- and the execution contract refuses
    it, correctly: an input alone does not show a calculation. Without
    the output this fixture made the control pass for the wrong reason
    (measured 2026-11-22, when the contract was mounted in grade_run).
    """
    return Evidence(
        workspace="ws",
        control=dict(control),
        steps=[{"program": "orca", "input": "a.inp", "output": "a.out",
                "terminated_normally": True,
                "charge": control.get("charge"),
                "multiplicity": control.get("multiplicity")}],
    )


# ------------------------------------------------------------- mount 1
def test_grade_run_flags_a_run_over_the_core_budget():
    """PAL 8 against a 2-core budget must fail the whole-run verdict."""
    ev = _ev(**{"charge": 0, "PAL": 8})
    verdict = w20_common.grade_run(ev, max_cores=2)
    assert verdict["verdict"] == "fail", verdict
    assert any("resource" in p for p in verdict["problems"])


def test_grade_run_flags_slurm_cores_over_the_budget():
    ev = _ev(charge=0)
    ev.slurm_records.append({"requested_cores": 4})
    verdict = w20_common.grade_run(ev, max_cores=2)
    assert verdict["verdict"] == "fail", verdict


def test_grade_run_control_within_budget_still_passes():
    """The control: a run inside the budget is unaffected by the mount."""
    ev = _ev(**{"charge": 0, "PAL": 2})
    verdict = w20_common.grade_run(ev, max_cores=2)
    assert verdict["verdict"] == "pass", verdict["problems"]


def test_grade_run_default_budget_is_the_practice_two_cores():
    """The signature default is max_cores=2 -- the practice budget the
    runner's call relies on; an omitted argument does not disable the
    check, it applies the declared budget. A run inside it passes."""
    ev = _ev(**{"charge": 0, "PAL": 2})
    verdict = w20_common.grade_run(ev)
    assert verdict["verdict"] == "pass", verdict["problems"]
    over = _ev(**{"charge": 0, "PAL": 8})
    assert w20_common.grade_run(over)["verdict"] == "fail"


# ------------------------------------------------------------- mount 2
def _d1_tasks():
    return [t for t in load_tasks(_ROOT / "delfin" / "agent" /
                                   "pack" / "benchmark" /
                                   "tasks_w17_data.yaml")]


def test_every_d1_task_declares_the_fixture_builder_as_setup():
    tasks = _d1_tasks()
    assert len(tasks) == 5, [t.id for t in tasks]
    for t in tasks:
        assert t.setup == "w17_data_fixtures.py", (t.id, t.setup)


def test_the_setup_key_runs_the_fixture_builder(tmp_path):
    """The declared setup is a real script that builds the CSVs the
    prompt names — the second answer is reached by the product."""
    from delfin.agent.benchmark_runner import run_setup

    ok, output = run_setup("w17_data_fixtures.py", tmp_path,
                           root=_ROOT)
    assert ok, output
    assert (tmp_path / "fixtures" / "trap_mixed_units.csv").is_file()
    assert (tmp_path / "fixtures" / "trap_hostile_cell.csv").is_file()


def test_the_setup_refuses_an_existing_fixture_dir(tmp_path):
    """Mirror of every other setup script: a second run fails instead of
    half-rebuilding, so the attempt reads as unmeasured."""
    from delfin.agent.benchmark_runner import run_setup

    ok, _ = run_setup("w17_data_fixtures.py", tmp_path, root=_ROOT)
    assert ok
    ok2, output2 = run_setup("w17_data_fixtures.py", tmp_path, root=_ROOT)
    assert not ok2, output2


def test_yaml_setup_key_names_the_module_the_parser_loads():
    """The pack names its setup modules BY FILENAME from tasks.yaml;
    this pins that the filename in the yaml is the module on disk."""
    raw = yaml.safe_load(
        (_ROOT / "delfin" / "agent" / "pack" / "benchmark" /
         "tasks_w17_data.yaml").read_text(encoding="utf-8"))
    setups = {t.get("setup") for t in raw["tasks"]}
    assert setups == {"w17_data_fixtures.py"}
