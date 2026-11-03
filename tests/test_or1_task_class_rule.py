"""OR1 tasks must not violate the chemistry cheap-physics rule.

The rule tests/test_chemistry_is_graded_on_the_physics.py declares
task_class chemistry as the cheap xtb-physics suite: no task in that
class may expect an ORCA run. OR1's five TS tasks DO expect real ORCA
(adapters orca_ts_*, graded on real freq outputs) — so they must not
carry task_class chemistry. Red test: fails while the yaml still says
chemistry, passes once reclassed (F-CI-1, operator finding (1) on the
PR-ready check at c7b6f13f).
"""
from pathlib import Path

import yaml

_YAML = Path(
    "delfin/agent/pack/benchmark/tasks_w20_or1.yaml")


def test_or1_tasks_do_not_use_the_chemistry_class():
    data = yaml.safe_load(_YAML.read_text(encoding="utf-8"))
    for task in data["tasks"]:
        assert task.get("task_class") != "chemistry", (
            f"{task['id']}: chemistry is the cheap xtb-physics class "
            "(test_chemistry_is_graded_on_the_physics forbids ORCA "
            "expectations there); OR1 tasks run real ORCA through "
            "DELFIN's adapters and must use a class the rule does not "
            "cover")


def test_or1_tasks_carry_the_reclassed_class():
    data = yaml.safe_load(_YAML.read_text(encoding="utf-8"))
    classes = {t["task_class"] for t in data["tasks"]}
    assert classes == {"science_analysis"}, classes
